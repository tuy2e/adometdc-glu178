"""Track B p20: firm up or break the Galdieria result of p18 (ODC-half dN/dS 0.20 with 7/13 catalytic sites lost).

A. Bootstrap over codons (N_BOOT replicates, seed 1) for every pair of p18_dnds.tsv whose three regions were all
   interpretable. NG86 (Nei-Gojobori 1986, Jukes-Cantor), computed from per-codon synonymous/non-synonymous site
   and difference counts (Biopython's own NG86 counting, summed over the resampled codons - identical to
   cal_dn_ds on the full set). 95 % CI = 2.5 / 97.5 percentiles. Replicates with dS = 0 or a saturated
   (not computable) dS or dN give no ratio and are counted and dropped.
   Tests, within a pair: ODC below 1 -> one-sided bootstrap p = share of replicates with ratio >= 1;
   ODC vs linker -> each replicate resamples both regions; p = 2 x the smaller share of (ODC - linker) <= 0 or >= 0.
   PAML/codeml is not installed, so there is no likelihood-ratio test.
B. Where is the constraint, for the same Galdieria pairs. Sliding windows over the ODC-half codons (WINDOW codons,
   step STEP, in the order of the first sequence). Then the ODC-half codons are split by the human-ODC position
   their first-sequence residue maps to (TM-align of the ODC half onto the AlphaFold model of human ODC, CA within
   4 A, as p9/p12):
     (i)   1D7K dimer-interface positions (contact 4.5 A, results/p12_interface_residues.tsv), not catalytic
     (ii)  the 13 catalytic positions
     (iii) everything else in the ODC half (also given restricted to mapped residues)
   and dN/dS with bootstrap CIs is computed for each; (i) vs (iii) tested as in A.
C. Comparison pairs with an intact ODC half. All 311 rows of p10_all_fusions.tsv, paired within the same family
   (lineage rank ending -aceae / -idae, from p10_fusion_search.tsv) and the same domain order. For every pair:
   protein identity, then (if both CDS verify, p18 method) whole-protein dS. Kept: 0.05 <= dS <= 1.0. ODC half
   intact = p10 structural status "intact", or, with no model, all 13 sites identical to human by p10's sequence
   read (labelled). A and B are then run on every kept pair whose two ODC halves are both intact; B's position map
   then uses a model in data/p10_models or Marc's models if one exists, otherwise a local sequence alignment
   to P11926 (labelled).

Writes results/p20_dnds_ci.tsv, results/p20_windows.tsv, results/p20_pairs_C.tsv, results/p20_report.md
Caches CDS in data/p20_cds/ (reads data/p18_cds/ first).
Run:  py -3.12 p20_dnds_ci.py            (A and B on Galdieria, then the C search; C fetches CDS for ~270 sequences)
      py -3.12 p20_dnds_ci.py --no-c     (A and B only)
"""
import csv
import itertools
import json
import math
import sys
import warnings
from collections import Counter, defaultdict

import numpy as np
from Bio import BiopythonExperimentalWarning
from Bio.Data import CodonTable

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.codonalign.codonseq import _count_diff_NG86, _count_site_NG86
from tmtools import tm_align

import p7_fusion_validation as p7
import p8_putrescine_site as p8
import p10_red_algal_pattern as p10
import p18_dnds as p18
from foldtrace.io import load_structure
from p2_fetch_annotations import DATA, RESULTS
from p3_residue_grid import RUN
from p6_decarboxylase_partner import references
from p9_self_processing import ODC_SITE_LABELS

sys.stdout.reconfigure(encoding="utf8")
N_BOOT, WINDOW, STEP, SEED = 1000, 30, 10, 1
DS_MIN, DS_MAX = 0.05, 1.0
CAT = {int(l[1:]) for l in ODC_SITE_LABELS}
TABLE = CodonTable.generic_by_id[1]
P20 = DATA / "p20_cds"
P20.mkdir(parents=True, exist_ok=True)
RUN_C = "--no-c" not in sys.argv


def cached(name, url, accept="application/json"):
    """p18's cache, extended: look in data/p18_cds first, write new records to data/p20_cds."""
    for d in (DATA / "p18_cds", P20):
        if (d / name).exists():
            return (d / name).read_text(encoding="utf8")
    body = p18.http(url, accept)
    (P20 / name).write_text(body, encoding="utf8")
    return body


p18.cached = cached
p18.CDS = P20


# ---------------------------------------------------------------- per-codon NG86
_site, _diff = {}, {}


def site(c):
    if c not in _site:
        _site[c] = _count_site_NG86([c], codon_table=TABLE)
    return _site[c]


def diff(a, b):
    if (a, b) not in _diff:
        _diff[(a, b)] = _count_diff_NG86(a, b, codon_table=TABLE)
    return _diff[(a, b)]


def codon_arrays(c1, c2):
    """Per codon pair: S sites, N sites (mean of the two codons), synonymous and non-synonymous differences."""
    S = np.array([(site(a)[0] + site(b)[0]) / 2 for a, b in zip(c1, c2)])
    N = np.array([(site(a)[1] + site(b)[1]) / 2 for a, b in zip(c1, c2)])
    D = np.array([diff(a, b) for a, b in zip(c1, c2)], float).reshape(-1, 2)
    return np.column_stack([S, N, D[:, 0], D[:, 1]])


def jc(p):
    return -0.75 * math.log(1 - 4 * p / 3) if p < 0.75 else None


def ng86(arr):
    """(dN, dS) from summed per-codon counts; None where not computable."""
    if arr.shape[0] == 0:
        return None, None
    S, N, Sd, Nd = arr.sum(0)
    return (jc(Nd / N) if N else None), (jc(Sd / S) if S else None)


def ratio(arr):
    dn, ds = ng86(arr)
    return dn / ds if dn is not None and ds else None


def boot(arr, rng, n=N_BOOT):
    """Bootstrap ratios (None where undefined)."""
    k = arr.shape[0]
    out = []
    for _ in range(n):
        out.append(ratio(arr[rng.integers(0, k, k)]) if k else None)
    return out


def ci(reps):
    v = np.array([x for x in reps if x is not None])
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (None, None)


def summarize(arr, reps):
    dn, ds = ng86(arr)
    lo, hi = ci(reps)
    good = [x for x in reps if x is not None]
    return {"codons": int(arr.shape[0]), "differing_codons": int(((arr[:, 2] + arr[:, 3]) > 0).sum()),
            "dN": r3(dn), "dS": r3(ds), "dNdS": r3(dn / ds) if dn is not None and ds else "",
            "ci_low": r3(lo), "ci_high": r3(hi), "boot_defined": len(good),
            "p_ratio_ge_1": r3(sum(x >= 1 for x in good) / len(good)) if good else "",
            "saturated": ds is None or ds > 1, "too_little_divergence": ds is not None and ds < DS_MIN}


def paired_p(r1, r2):
    d = [a - b for a, b in zip(r1, r2) if a is not None and b is not None]
    if not d:
        return None, 0
    lo = sum(x <= 0 for x in d) / len(d)
    hi = sum(x >= 0 for x in d) / len(d)
    return min(1.0, 2 * min(lo, hi)), len(d)


def r3(x):
    return "" if x is None else round(x, 3)


# ---------------------------------------------------------------- pairs, codons, regions
def regions_any(arch, length):
    doms = [(d.split()[0], *map(int, d.split()[1].split("-"))) for d in arch.split(" > ")]
    ado = [(s, e) for p, s, e in doms if p in {"PF01536", "PF02675"}]
    plp = [(s, e) for p, s, e in doms if p in {"PF02784", "PF00278"}]
    a0, a1 = min(s for s, _ in ado), max(e for _, e in ado)
    if min(s for s, _ in plp) > a1:        # AdoMetDC first: ODC = first PLP start .. C-terminus (as p18)
        p0 = min(s for s, _ in plp if s > a1)
        return {"AdoMetDC": (a0, a1), "linker": (a1 + 1, p0 - 1), "ODC": (p0, length)}
    p1 = max(e for _, e in plp if e < a0)  # PLP first: ODC = N-terminus .. last PLP end before AdoMetDC
    return {"ODC": (1, p1), "linker": (p1 + 1, a0 - 1), "AdoMetDC": (a0, a1)}


def pair_codons(sx, sy, cx, cy):
    """[(pos1, pos2, codon1, codon2)] for gap-free, stop-free, unambiguous aligned codon pairs."""
    a1, a2 = p18.align(sx, sy)
    out = []
    for p1, p2, _, _ in p18.columns(a1, a2):
        if p1 is None or p2 is None:
            continue
        k1, k2 = cx[3 * (p1 - 1):3 * p1], cy[3 * (p2 - 1):3 * p2]
        if len(k1) == len(k2) == 3 and not ({k1, k2} & p18.STOPS) and set(k1 + k2) <= set("ACGT"):
            out.append((p1, p2, k1, k2))
    return out, a1, a2


def protein_identity(a1, a2):
    cols = [(x, y) for x, y in zip(a1, a2) if x != "-" and y != "-"]
    return round(100 * sum(x == y for x, y in cols) / len(cols), 1) if cols else 0.0


def region_arr(cod, lo, hi, keep=None):
    sel = [(c1, c2) for p1, _, c1, c2 in cod if lo <= p1 <= hi and (keep is None or keep(p1))]
    return codon_arrays([a for a, _ in sel], [b for _, b in sel]) if sel else np.zeros((0, 4))


# ---------------------------------------------------------------- mapping to human ODC numbering
def human_map(acc, seq, lo, hi, odc, odc_seq):
    """{residue number in acc: human ODC position}, and the method used."""
    for path in (RUN / "models" / f"{acc}.pdb", p10.MODELS / f"AF-{acc}.pdb"):
        if path.exists():
            ca, s, nums = p7.sub_structure(load_structure(str(path)), lo, hi)
            oca, oseq, onums = odc
            r = tm_align(oca, ca, oseq, s)
            on = oca @ np.asarray(r.u).T + np.asarray(r.t)
            m, i, j = {}, 0, 0
            for x, y in zip(r.seqxA, r.seqyA):
                if x != "-" and y != "-" and np.linalg.norm(on[i] - ca[j]) <= p8.OFFSET:
                    m[nums[j]] = onums[i]
                i += x != "-"
                j += y != "-"
            return m, f"structure ({path.name}, TM-align, CA <= {p8.OFFSET} A)"
    a = p10.LOCAL.align(seq[lo - 1:hi], odc_seq)[0]
    m = {lo + i: j + 1 for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))}
    return m, "sequence (local alignment to P11926; no model)"


def interface_set():
    with open(RESULTS / "p12_interface_residues.tsv", encoding="utf8") as f:
        return {int(r["pos"]) for r in csv.DictReader(f, delimiter="\t")
                if r["enzyme"] == "ODC" and r["contact_4.5"] == "True"}


# ---------------------------------------------------------------- A + B for one pair
def analyse(x, y, seqs, cds, reg, rng, odc, odc_seq, iface, label):
    cod, _, _ = pair_codons(seqs[x], seqs[y], cds[x], cds[y])
    rows, wins = [], []
    arrs = {name: region_arr(cod, lo, hi) for name, (lo, hi) in reg.items()}
    reps = {name: boot(a, rng) for name, a in arrs.items()}
    for name in ("AdoMetDC", "linker", "ODC"):
        rows.append({"set": label, "seq1": x, "seq2": y, "part": "A", "region": name, **summarize(arrs[name], reps[name])})
    p, n = paired_p(reps["ODC"], reps["linker"])
    rows[-1]["p_vs_linker"] = r3(p)
    rows[-1]["paired_reps"] = n
    lo, hi = reg["ODC"]
    m, method = human_map(x, seqs[x], lo, hi, odc, odc_seq)
    groups = {"(i) interface, not catalytic": lambda q: m.get(q) in iface - CAT,
              "(ii) catalytic (13)": lambda q: m.get(q) in CAT,
              "(iii) everything else": lambda q: not (m.get(q) in iface or m.get(q) in CAT),
              "(iii) everything else, mapped only": lambda q: q in m and not (m[q] in iface or m[q] in CAT)}
    garr = {g: region_arr(cod, lo, hi, f) for g, f in groups.items()}
    grep = {g: boot(a, rng) for g, a in garr.items()}
    for g in groups:
        rows.append({"set": label, "seq1": x, "seq2": y, "part": "B", "region": f"ODC {g}", "map_method": method,
                     **summarize(garr[g], grep[g])})
    p, n = paired_p(grep["(i) interface, not catalytic"], grep["(iii) everything else"])
    rows[-4]["p_vs_iii"] = r3(p)
    rows[-4]["paired_reps"] = n
    odc_cod = [c for c in cod if lo <= c[0] <= hi]
    for k in range(0, max(1, len(odc_cod) - WINDOW + 1), STEP):
        w = odc_cod[k:k + WINDOW]
        arr = codon_arrays([c[2] for c in w], [c[3] for c in w])
        dn, ds = ng86(arr)
        hum = [m.get(c[0]) for c in w if m.get(c[0])]
        wins.append({"set": label, "seq1": x, "seq2": y, "window": k // STEP + 1,
                     "seq1_residues": f"{w[0][0]}-{w[-1][0]}", "human_odc_span": f"{min(hum)}-{max(hum)}" if hum else "",
                     "codons": len(w), "dN": r3(dn), "dS": r3(ds), "dNdS": r3(dn / ds) if dn is not None and ds else "",
                     "status": "saturated" if ds is None or ds > 1 else "too little divergence" if ds < DS_MIN else "ok",
                     "interface_codons": sum(m.get(c[0]) in iface - CAT for c in w),
                     "catalytic_codons": sum(m.get(c[0]) in CAT for c in w)})
    return rows, wins


# ---------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(SEED)
    odc_seq = references()["P11926"]
    odc = p7.sub_structure(load_structure(str(p7.af_model("P11926"))), 1, p7.ODC_ORDERED)
    iface = interface_set()
    cand = p18.candidates()
    ents = json.loads(next((DATA / "p10").glob("uniparc__feature_id_PF01536*.json")).read_text(encoding="utf8"))
    seq_of = {e["uniParcId"]: e["sequence"]["value"] for e in ents}
    seqs, cds, notes = {}, {}, {}

    def load(acc):
        if acc in seqs:
            return
        r = cand[acc]
        seqs[acc] = seq_of[r["upi"]]
        cds[acc], notes[acc] = p18.cds_for(r["upi"], seqs[acc])

    # A + B on p18's interpretable pairs
    p18rows = list(csv.DictReader(open(RESULTS / "p18_dnds.tsv", encoding="utf8"), delimiter="\t"))
    good = [r for r in p18rows if all(r.get(f"{g}_interpretable") == "True" for g in ("AdoMetDC", "linker", "ODC"))]
    out, wins = [], []
    for r in good:
        x, y = r["seq1"], r["seq2"]
        load(x)
        load(y)
        reg = p18.regions(cand[x]["pfam_architecture"], len(seqs[x]))
        rr, ww = analyse(x, y, seqs, cds, reg, rng, odc, odc_seq, iface, "Galdieria (p18)")
        out += rr
        wins += ww
        print(f"A/B {x} - {y}: ODC {rr[2]['dNdS']} [{rr[2]['ci_low']}, {rr[2]['ci_high']}]", flush=True)

    pairs_c = []
    if RUN_C:
        allf = list(csv.DictReader(open(RESULTS / "p10_all_fusions.tsv", encoding="utf8"), delimiter="\t"))
        info = {}
        for f in allf:
            s = cand[f["accession"].split(".")[0]]
            lin = s["lineage"].split("; ")
            fam = next((x for x in reversed(lin) if x.endswith(("aceae", "idae"))), "")
            sites = f["sequence_sites"]
            n_eq = sum(f"{l}=" in sites for l in ODC_SITE_LABELS)
            intact = ("structure" if f["odc_status"] == "intact" else
                      "sequence (13/13 identical, no model)" if not f["odc_status"] and n_eq == 13 else "")
            info[f["accession"].split(".")[0]] = {"family": fam, "order": s["adometdc_then_plp"], "group": f["group"],
                                                  "organism": f["organism"], "odc_status": f["odc_status"] or f"sequence {n_eq}/13",
                                                  "intact": intact, "arch": s["pfam_architecture"]}
        byfam = defaultdict(list)
        for a, i in info.items():
            if i["family"]:
                byfam[(i["family"], i["order"])].append(a)
        cand_pairs = [(a, b) for k, v in byfam.items() for a, b in itertools.combinations(sorted(v), 2)]
        print(f"C: {len(cand_pairs)} candidate pairs in {sum(len(v) > 1 for v in byfam.values())} family groups", flush=True)
        for n, (a, b) in enumerate(cand_pairs, 1):
            load(a)
            load(b)
            row = {"seq1": a, "organism1": info[a]["organism"], "seq2": b, "organism2": info[b]["organism"],
                   "family": info[a]["family"], "group": info[a]["group"],
                   "odc_half_1": info[a]["odc_status"], "odc_half_2": info[b]["odc_status"],
                   "both_intact": bool(info[a]["intact"] and info[b]["intact"]),
                   "intact_basis": f"{info[a]['intact'] or '-'} / {info[b]['intact'] or '-'}"}
            a1, a2 = p18.align(seqs[a], seqs[b])
            row["protein_identity"] = protein_identity(a1, a2)
            if cds[a] and cds[b]:
                cod, _, _ = pair_codons(seqs[a], seqs[b], cds[a], cds[b])
                dn, ds = ng86(codon_arrays([c[2] for c in cod], [c[3] for c in cod]))
                row.update({"codons": len(cod), "dN_whole": r3(dn), "dS_whole": r3(ds) if ds is not None else "not computable",
                            "in_window": ds is not None and DS_MIN <= ds <= DS_MAX})
            else:
                row.update({"dS_whole": "no verified CDS", "in_window": False})
            pairs_c.append(row)
            if n % 100 == 0:
                print(f"  C pairs {n}/{len(cand_pairs)}", flush=True)
        for row in pairs_c:
            if row["in_window"] and row["both_intact"]:
                a, b = row["seq1"], row["seq2"]
                try:
                    reg = regions_any(info[a]["arch"], len(seqs[a]))
                except ValueError:
                    row["note"] = "regions not definable from Pfam"
                    continue
                rr, ww = analyse(a, b, seqs, cds, reg, rng, odc, odc_seq, iface, f"C: {row['family']}")
                out += rr
                wins += ww
        with open(RESULTS / "p20_pairs_C.tsv", "w", encoding="utf8", newline="") as f:
            cols = list(dict.fromkeys(k for r in pairs_c for k in r))
            w = csv.DictWriter(f, cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(pairs_c)
    for name, rows in (("p20_dnds_ci.tsv", out), ("p20_windows.tsv", wins)):
        with open(RESULTS / name, "w", encoding="utf8", newline="") as f:
            cols = list(dict.fromkeys(k for r in rows for k in r))
            w = csv.DictWriter(f, cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(rows)
    json.dump({"cds_notes": notes}, open(P20 / "_p20_cds_notes.json", "w", encoding="utf8"), indent=1)


if __name__ == "__main__":
    main()
