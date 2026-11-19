"""Stage 1: the six pre-registered pocket measurements.

Measures exactly what 1_reports/PREREG_pocket_chemistry_measurements.md fixes in advance:

  M1 electrostatic potential at the putrescine position   PDB2PQR (PARSE, pH 7) + APBS, kT/e
  M2 pKa at human 15 / 174 / 178 / 256                    PROPKA3
  M3 hydrogen-bond partners of the 178 side chain         geometry on the protonated structure
  M4 basic residues filling the site                      Arg/Lys/His N within 8 A
  M5 free cavity volume                                   3 A probe within 4 A (as 10 Nov)
  M6 docking in-site fraction                             Vina 1.2.5, 12 A box, 9 modes

Usage:  python stage1_measure.py <aligned_dir> <ligand.pdbqt> [names...]
Writes <aligned_dir>/stage1_measurements.json and .tsv

Needs: pdb2pqr30, apbs, propka3, obabel, vina on PATH.
"""
import glob
import json
import os
import re
import subprocess
import sys
import warnings

import numpy as np
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

PUT_CENTRE = np.array([61.324, 11.189, 7.513])
GATE_POS = [15, 174, 178, 256]
parser = PDBParser(QUIET=True)

BASIC_N = {"ARG": ["NE", "NH1", "NH2"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# ---------------------------------------------------------------- M1: APBS ----
def read_dx(path):
    """Minimal OpenDX reader: returns (origin, delta, shape, values)."""
    origin = delta = None
    shape = None
    vals = []
    with open(path) as fh:
        deltas = []
        for line in fh:
            if line.startswith("object 1"):
                shape = tuple(int(x) for x in line.split()[-3:])
            elif line.startswith("origin"):
                origin = np.array([float(x) for x in line.split()[1:4]])
            elif line.startswith("delta"):
                deltas.append([float(x) for x in line.split()[1:4]])
            elif line and line[0].isdigit() or line.startswith("-"):
                vals.extend(float(x) for x in line.split())
        delta = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    grid = np.array(vals[: shape[0] * shape[1] * shape[2]]).reshape(shape)
    return origin, delta, grid


def interp(origin, delta, grid, pt):
    f = (pt - origin) / delta
    i = np.floor(f).astype(int)
    if np.any(i < 0) or np.any(i + 1 >= np.array(grid.shape)):
        return float("nan")
    d = f - i
    out = 0.0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((1 - d[0]) if dx == 0 else d[0]) * ((1 - d[1]) if dy == 0 else d[1]) \
                    * ((1 - d[2]) if dz == 0 else d[2])
                out += w * grid[i[0] + dx, i[1] + dy, i[2] + dz]
    return float(out)


def apbs_potential(pdb, workdir, name):
    pqr = os.path.join(workdir, "%s.pqr" % name)
    r = run(["pdb2pqr30", "--ff=PARSE", "--titration-state-method=propka",
             "--with-ph=7.0", "--keep-chain", pdb, pqr])
    if not os.path.exists(pqr):
        return None, None, None, "pdb2pqr failed: %s" % r.stderr.strip().splitlines()[-1:]
    xyz = []
    with open(pqr) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz = np.array(xyz)
    ext = xyz.max(0) - xyz.min(0)
    cg = ext + 40.0
    inp = os.path.join(workdir, "%s.in" % name)
    dx = os.path.join(workdir, "%s_pot" % name)
    with open(inp, "w") as fh:
        fh.write("""read
    mol pqr %s
end
elec name solv
    mg-auto
    dime 97 97 97
    cglen %.1f %.1f %.1f
    fglen 40.0 40.0 40.0
    cgcent mol 1
    fgcent %.3f %.3f %.3f
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    ion charge  1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 2.0
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy no
    calcforce no
    write pot dx %s
end
quit
""" % (os.path.basename(pqr), cg[0], cg[1], cg[2],
            PUT_CENTRE[0], PUT_CENTRE[1], PUT_CENTRE[2], os.path.basename(dx)))
    r = run(["apbs", os.path.basename(inp)], cwd=workdir)
    # APBS 3.x writes <stem>.dx for a serial run and <stem>-PE0.dx when built with MPI
    dxf = next((c for c in (dx + ".dx", dx + "-PE0.dx") if os.path.exists(c)), dx + ".dx")
    if not os.path.exists(dxf):
        return None, None, pqr, "apbs failed: %s" % (r.stdout or r.stderr)[-300:]
    origin, delta, grid = read_dx(dxf)
    at_centre = interp(origin, delta, grid, PUT_CENTRE)
    # mean over a 2 A sphere, sampled on a 0.5 A grid
    g = np.arange(-2.0, 2.01, 0.5)
    pts = np.array([[x, y, z] for x in g for y in g for z in g])
    pts = pts[np.linalg.norm(pts, axis=1) <= 2.0] + PUT_CENTRE
    sphere = [interp(origin, delta, grid, p) for p in pts]
    sphere = [v for v in sphere if not np.isnan(v)]
    return at_centre, (float(np.mean(sphere)) if sphere else None), pqr, None


# -------------------------------------------------------------- M2: PROPKA ----
def propka_pkas(pdb, workdir, resmap, chain):
    r = run(["propka3", os.path.abspath(pdb)], cwd=workdir)
    pka = os.path.join(workdir, os.path.basename(pdb)[:-4] + ".pka")
    if not os.path.exists(pka):
        return {}, "propka failed: %s" % (r.stderr or r.stdout)[-200:]
    want = {}
    for human in GATE_POS:
        tag = resmap.get(str(human), "-")
        if tag != "-" and re.match(r"^[A-Z]\d+$", tag):
            want[int(tag[1:])] = human
    out, insummary = {}, False
    for line in open(pka):
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            insummary = True
            continue
        if insummary:
            parts = line.split()
            if len(parts) >= 4 and parts[1].isdigit():
                num = int(parts[1])
                if num in want and (chain is None or parts[2] == chain):
                    try:
                        out[want[num]] = dict(res=parts[0], pka=float(parts[3]))
                    except ValueError:
                        pass
            elif line.strip().startswith("-") and out:
                break
    return out, None


# ------------------------------------------------- M3/M4/M5: geometry ---------
def site_geometry(pdb, pqr, resmap, chain):
    st = parser.get_structure("s", pdb)[0]
    atoms = [a for a in st.get_atoms() if a.element != "H"]
    xyz = np.array([a.coord for a in atoms])

    # M4 basic residues within 8 A
    basic, named = [], {}
    for a, c in zip(atoms, xyz):
        r = a.get_parent()
        rn = r.get_resname()
        if rn in BASIC_N and a.get_id() in BASIC_N[rn] and np.linalg.norm(c - PUT_CENTRE) <= 8.0:
            tag = "%s%d%s" % (rn.capitalize(), r.id[1], r.get_parent().id)
            if tag not in named:
                named[tag] = round(float(np.linalg.norm(c - PUT_CENTRE)), 1)
    basic = sorted(named.items(), key=lambda kv: kv[1])

    # M5 cavity, identical to the 10 Nov docking round
    g = np.arange(-4.0, 4.01, 0.5)
    pts = np.array([[x, y, z] for x in g for y in g for z in g])
    pts = pts[np.linalg.norm(pts, axis=1) <= 4.0] + PUT_CENTRE
    near = xyz[np.linalg.norm(xyz - PUT_CENTRE, axis=1) < 8.0]
    free = sum(1 for q in pts if np.min(np.linalg.norm(near - q, axis=1)) >= 3.0)
    cavity = round(free * 0.125, 1)

    # M3 H-bonds of the 178 side chain, on the protonated PQR
    hb = []
    tag178 = resmap.get("178", "-")
    if tag178 != "-" and os.path.exists(pqr):
        target = int(tag178[1:])
        prot = []
        with open(pqr) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    prot.append(dict(
                        name=line[12:16].strip(), resn=line[17:20].strip(),
                        ch=line[21].strip(), resi=int(line[22:26]),
                        xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
        sidechain = [a for a in prot if a["resi"] == target and (chain is None or a["ch"] == chain)
                     and a["name"] not in BACKBONE and a["name"][0] in "NO"]
        polar = [a for a in prot if a["name"][0] in "NO"
                 and not (a["resi"] == target and (chain is None or a["ch"] == chain))]
        hyd = [a for a in prot if a["name"][0] == "H"]
        seen = set()
        for s in sidechain:
            for p in polar:
                d = np.linalg.norm(s["xyz"] - p["xyz"])
                if d > 3.5:
                    continue
                ok = False
                for donor, acceptor in ((s, p), (p, s)):
                    for h in hyd:
                        if np.linalg.norm(h["xyz"] - donor["xyz"]) > 1.3:
                            continue
                        v1 = donor["xyz"] - h["xyz"]
                        v2 = acceptor["xyz"] - h["xyz"]
                        ang = np.degrees(np.arccos(np.clip(
                            np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)))
                        if ang >= 120:
                            ok = True
                            break
                    if ok:
                        break
                if ok:
                    tag = "%s%d%s" % (p["resn"].capitalize(), p["resi"], p["ch"])
                    if tag not in seen:
                        seen.add(tag)
                        hb.append((tag, round(float(d), 2)))
    return dict(basic=basic, cavity=cavity, hbonds=sorted(hb, key=lambda x: x[1]))


# ------------------------------------------------------------- M6: docking ----
def dock(pdb, ligand, workdir, name):
    rec = os.path.join(workdir, "%s_rec.pdbqt" % name)
    if not os.path.exists(rec):
        run(["obabel", pdb, "-O", rec, "-xr", "-p", "7.4", "--partialcharge", "gasteiger"])
    if not os.path.exists(rec):
        return None, "obabel failed to make a receptor"
    out = os.path.join(workdir, "%s_site.pdbqt" % name)
    log = run(["vina", "--receptor", rec, "--ligand", ligand,
               "--center_x", "%.3f" % PUT_CENTRE[0], "--center_y", "%.3f" % PUT_CENTRE[1],
               "--center_z", "%.3f" % PUT_CENTRE[2],
               "--size_x", "12", "--size_y", "12", "--size_z", "12",
               "--exhaustiveness", "32", "--num_modes", "9", "--seed", "20261110",
               "--out", out])
    if not os.path.exists(out):
        return None, "vina failed: %s" % (log.stdout or log.stderr)[-300:]
    poses, cur, score = [], [], None
    for line in open(out):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("REMARK VINA RESULT"):
            score = float(line.split()[3])
        elif line.startswith(("ATOM", "HETATM")):
            el = (line[77:79].strip() or line[12:14].strip())
            if el[0] in "CN":
                cur.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("ENDMDL"):
            poses.append((score, np.array(cur).mean(0)))
    if not poses:
        return None, "no poses"
    dists = [float(np.linalg.norm(c - PUT_CENTRE)) for _, c in poses]
    insite = [i for i, d in enumerate(dists) if d <= 2.5]
    best = int(np.argmin([s for s, _ in poses]))
    return dict(n_poses=len(poses), n_in_site=len(insite),
                best_score=poses[best][0], best_dist=round(dists[best], 1),
                best_in_site_score=(min(poses[i][0] for i in insite) if insite else None)), None


# ---------------------------------------------------------------------- main --
def main(aligned_dir, ligand, names):
    work = os.path.join(aligned_dir, "_work")
    os.makedirs(work, exist_ok=True)
    prep = {}
    for f in ("prep_stage1.json", "prep.json"):
        p = os.path.join(aligned_dir, f)
        if os.path.exists(p):
            prep.update(json.load(open(p)))

    if not names:
        names = sorted(os.path.basename(p)[:-12]
                       for p in glob.glob(os.path.join(aligned_dir, "*_aligned.pdb")))

    rows = {}
    for name in names:
        pdb = os.path.join(aligned_dir, "%s_aligned.pdb" % name)
        if not os.path.exists(pdb):
            print("  %-14s missing %s" % (name, pdb))
            continue
        info = prep.get(name, {})
        resmap = info.get("residues") or info.get("res") or {}
        chain = info.get("chain_used")
        notes = []

        at_c, mean_2a, pqr, err = apbs_potential(pdb, work, name)
        if err:
            notes.append(err)
        pkas, err = propka_pkas(pdb, work, resmap, chain)
        if err:
            notes.append(err)
        geom = site_geometry(pdb, pqr or "", resmap, chain)
        dk, err = dock(pdb, ligand, work, name)
        if err:
            notes.append(err)

        rows[name] = dict(
            label=info.get("label", name), residues=resmap,
            M1_centre=None if at_c is None else round(at_c, 2),
            M1_mean2A=None if mean_2a is None else round(mean_2a, 2),
            M2_pka={str(k): v["pka"] for k, v in pkas.items()},
            M3_hbonds_178=geom["hbonds"],
            M4_basic_8A=geom["basic"],
            M5_cavity=geom["cavity"],
            M6=dk, notes=notes)
        print("  %-14s pot %8s kT/e (2A mean %8s) | cavity %6.1f A3 | basic %d | dock %s | %s"
              % (name, rows[name]["M1_centre"], rows[name]["M1_mean2A"], geom["cavity"],
                 len(geom["basic"]),
                 ("%d/9 @ %.2f" % (dk["n_in_site"], dk["best_score"])) if dk else "failed",
                 "; ".join(notes) if notes else "ok"))

    with open(os.path.join(aligned_dir, "stage1_measurements.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    cols = ["name", "M1_centre", "M1_mean2A", "M5_cavity", "M6_in_site", "M6_best",
            "M6_best_dist", "M4_n_basic", "M2_178", "M2_174", "M3_n_hbonds"]
    with open(os.path.join(aligned_dir, "stage1_measurements.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for n, r in rows.items():
            d = r["M6"] or {}
            fh.write("\t".join(str(x) for x in [
                n, r["M1_centre"], r["M1_mean2A"], r["M5_cavity"],
                d.get("n_in_site"), d.get("best_score"), d.get("best_dist"),
                len(r["M4_basic_8A"]), r["M2_pka"].get("178"), r["M2_pka"].get("174"),
                len(r["M3_hbonds_178"])]) + "\n")
    print("\nWrote stage1_measurements.json / .tsv in %s" % aligned_dir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
