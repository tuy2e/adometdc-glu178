"""Confounder check: is the site potential just the protein's overall charge?

M1 measures the electrostatic potential at the putrescine position. A protein that is more
basic overall will have a less negative potential *everywhere*, so a difference between two
species at the site could be nothing more than a difference in bulk amino-acid composition.
If that is what is happening, M1 separates proteomes, not pockets, and the gate would pass
for the wrong reason.

This script re-reads the APBS maps that stage1_measure.py already wrote and computes, per
structure:

  net_charge        formal charge of the whole protein, from the PQR
  charge_15A        summed formal charge within 15 A of the putrescine position
  shell_median      median potential over a solvent shell around the WHOLE protein
                    (points 3.0-5.0 A outside any heavy atom, sampled evenly)
  site_minus_shell  M1 at the site minus that shell median  <- the confound-corrected value
  site_z            how far the site sits from the shell distribution, in its own SDs

If the activated / not-stimulated separation survives in site_minus_shell, it is a property
of the pocket. If it collapses, M1 is reporting bulk charge and must not be used in the gate.

Usage:  python stage1_charge_control.py <aligned_dir> [names...]
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage1_measure import PUT_CENTRE, interp, read_dx  # noqa: E402

RNG = np.random.default_rng(20261111)


def pqr_atoms(path):
    xyz, q, r = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                q.append(float(parts[-2]))
                r.append(float(parts[-1]))
    return np.array(xyz), np.array(q), np.array(r)


def shell_points(xyz, n=4000, lo=3.0, hi=5.0):
    """Points in a thin solvent shell outside the protein, sampled from its bounding box."""
    lo_c, hi_c = xyz.min(0) - hi - 1, xyz.max(0) + hi + 1
    keep = []
    for _ in range(40):
        cand = RNG.uniform(lo_c, hi_c, size=(20000, 3))
        # distance to nearest heavy atom, in chunks to stay in memory
        d = np.full(len(cand), np.inf)
        for i in range(0, len(xyz), 2000):
            block = xyz[i:i + 2000]
            d = np.minimum(d, np.linalg.norm(cand[:, None, :] - block[None, :, :], axis=2).min(1))
        sel = cand[(d >= lo) & (d <= hi)]
        keep.append(sel)
        if sum(len(k) for k in keep) >= n:
            break
    pts = np.vstack(keep)
    if len(pts) > n:
        pts = pts[RNG.choice(len(pts), n, replace=False)]
    return pts


def main(aligned_dir, names):
    work = os.path.join(aligned_dir, "_work")
    meas_path = os.path.join(aligned_dir, "stage1_measurements.json")
    meas = json.load(open(meas_path)) if os.path.exists(meas_path) else {}
    if not names:
        names = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(work, "*.pqr")))

    out = {}
    for name in names:
        pqr = os.path.join(work, "%s.pqr" % name)
        dx = next((c for c in (os.path.join(work, "%s_pot.dx" % name),
                               os.path.join(work, "%s_pot-PE0.dx" % name)) if os.path.exists(c)), None)
        if not (os.path.exists(pqr) and dx):
            print("  %-14s missing pqr or dx - skipped" % name)
            continue
        xyz, q, _ = pqr_atoms(pqr)
        origin, delta, grid = read_dx(dx)

        site = interp(origin, delta, grid, PUT_CENTRE)
        d_site = np.linalg.norm(xyz - PUT_CENTRE, axis=1)
        charge_15 = float(q[d_site <= 15.0].sum())

        pts = shell_points(xyz)
        vals = np.array([interp(origin, delta, grid, p) for p in pts])
        vals = vals[~np.isnan(vals)]
        med = float(np.median(vals))
        sd = float(np.std(vals))

        out[name] = dict(
            net_charge=round(float(q.sum()), 1),
            charge_15A=round(charge_15, 1),
            n_atoms=len(xyz),
            shell_median=round(med, 2),
            shell_sd=round(sd, 2),
            site=round(site, 2),
            site_minus_shell=round(site - med, 2),
            site_z=round((site - med) / sd, 2) if sd else None,
            m5_cavity=(meas.get(name) or {}).get("M5_cavity"),
        )
        r = out[name]
        print("  %-14s net %+6.1f e | charge within 15 A %+5.1f | shell %+6.2f +/-%4.2f | "
              "site %+6.2f | site-shell %+6.2f | z %+5.2f"
              % (name, r["net_charge"], r["charge_15A"], r["shell_median"], r["shell_sd"],
                 r["site"], r["site_minus_shell"], r["site_z"]))

    with open(os.path.join(aligned_dir, "stage1_charge_control.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nWrote stage1_charge_control.json")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
