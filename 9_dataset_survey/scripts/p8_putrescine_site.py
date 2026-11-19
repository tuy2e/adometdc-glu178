"""Track B step 1: does the AdoMetDC half of the fusions keep the putrescine regulatory site?

Positions (human P17707 numbering, 1MSV chain A): 13, 15, 111, 174, 176, 178, 256, 285
(wild-type L13, E15, F111, D174, T176, E178, E256, F285; read from P17707, not assumed).

Mapping: exactly the p3 method (tmtools TM-align of the whole model onto 1MSV chain A, CA traces,
foldtrace's loader, alignment cached in data/p3_alignments/ and shared with p3).
Per protein and position:
    call     identical / substituted / displaced (>4.0 A) / deleted (internal gap) / not_covered
    state    kept        identical, residue pLDDT >= 70
             lost        substituted, residue pLDDT >= 70
             unresolved  pLDDT < 70, or displaced, or no residue aligned (deleted / not_covered)
             (same rule as FOLDTRACE: a displaced or unaligned position is unresolved, not lost;
              the call column keeps the detail so either reading can be checked)
    in_adometdc_pfam   whether the mapped residue lies inside an AdoMetDC Pfam match (PF01536/PF02675),
                       i.e. that the site was read from the AdoMetDC half and not from the partner

Groups (from p5_loss_vs_fusion_proteins.tsv; fold-passers, A0A8C5BDJ7 excluded):
    standalone   590
    odc_fusion    23   (the PLP-decarboxylase fusions of p6/p7)
    fused_all    101   (every fused architecture, the 23 included)
    fused_other   78   (fused_all minus the 23)
Tests, each position separately, kept vs lost (unresolved left out of the denominator):
    Fisher exact (two-sided), and within-phylum permutation exactly as p5 (labels shuffled within
    phylum, 10,000 permutations, seed 1, statistic = difference in loss rate), plus MH odds ratio.
Every test repeated for the InterPro-only set (annotation source = interpro).

Reads (never writes) the FOLDTRACE run folder.
Writes results/p8_putrescine_site.tsv (one row per protein x position),
       results/p8_putrescine_site_tests.tsv, results/p8_putrescine_site_report.txt
Run:  py -3.12 p8_putrescine_site.py          (first run aligns ~690 models, several minutes;
                                               re-runs use the cache)
      py -3.12 p8_putrescine_site.py --only-23  (quick: the 23 fusions only, no tests)
"""
import csv
import random
import sys
from collections import Counter, defaultdict

import numpy as np
from foldtrace.io import load_structure

from p3_residue_grid import RUN, REF_PDB, REF_CHAIN, OFFSET, PLDDT_OK, align, human_sequence, plddt_by_resnum
from p5_loss_vs_fusion import RESULTS, fisher, mantel_haenszel, read

POSITIONS = [13, 15, 111, 174, 176, 178, 256, 285]
ADOMETDC = ("PF01536", "PF02675")
N_PERM = 10000
STATES = ("kept", "lost", "unresolved")


def adometdc_spans(row):
    out = []
    for part in row["domains_with_positions"].split("; "):
        bits = part.split()
        if len(bits) == 2 and bits[0] in ADOMETDC and "-" in bits[1]:
            a, b = bits[1].split("-")
            out.append((int(a), int(b)))
    return out


def map_sites(acc, ref, ref_idx, human):
    mpath = RUN / "models" / f"{acc}.pdb"
    cand = load_structure(str(mpath))
    pl = plddt_by_resnum(mpath)
    a = align(acc, ref, cand)
    u, t = np.asarray(a["u"]), np.asarray(a["t"])
    ref_on_cand = ref.ca @ u.T + t
    pairs, first, last = {}, None, None
    i = j = 0
    for x, y in zip(a["seqxA"], a["seqyA"]):
        if x != "-" and y != "-":
            pairs[i] = j
            first = i if first is None else first
            last = i
        i += x != "-"
        j += y != "-"
    out = {}
    for pos in POSITIONS:
        ri = ref_idx[pos]
        h = human[pos - 1]
        if ri in pairs:
            ci = pairs[ri]
            off = float(np.linalg.norm(ref_on_cand[ri] - cand.ca[ci]))
            aa, rn, p = cand.seq[ci], cand.resnums[ci], pl.get(cand.resnums[ci])
            call = "displaced" if off > OFFSET else ("identical" if aa == h else "substituted")
        else:
            aa = rn = off = p = None
            call = "deleted" if first is not None and first < ri < last else "not_covered"
        if call in ("identical", "substituted") and p is not None and p >= PLDDT_OK:
            state = "kept" if call == "identical" else "lost"
        else:
            reason = {"displaced": "displaced", "deleted": "no_residue_aligned", "not_covered": "no_residue_aligned"}
            state = "unresolved"
            why = reason.get(call, "plddt_below_70")
            if call in ("displaced",) and p is not None and p < PLDDT_OK:
                why = "displaced;plddt_below_70"
        out[pos] = {"human_aa": h, "call": call, "cand_resnum": rn, "cand_aa": aa,
                    "ca_offset": None if off is None else round(off, 2), "plddt": p, "state": state,
                    "unresolved_reason": "" if state != "unresolved" else why}
    return out, a["tm_norm_ref"]


def perm_test(rows, g1, g2, rank="phylum"):
    """p5 method: shuffle labels within taxon; statistic = loss rate g1 - loss rate g2 (resolved only)."""
    strata = defaultdict(list)
    for r in rows:
        strata[r["tax_" + rank]].append((r["state"] == "lost", r["group"]))
    groups = list(strata.values())

    def stat(gs):
        n1 = sum(1 for g in gs for _, c in g if c == g1)
        n2 = sum(1 for g in gs for _, c in g if c == g2)
        l1 = sum(1 for g in gs for l, c in g if c == g1 and l)
        l2 = sum(1 for g in gs for l, c in g if c == g2 and l)
        return l1 / n1 - l2 / n2
    obs = stat(groups)
    rng = random.Random(1)
    hits = 0
    for _ in range(N_PERM):
        perm = []
        for g in groups:
            labels = [c for _, c in g]
            rng.shuffle(labels)
            perm.append([(l, c) for (l, _), c in zip(g, labels)])
        hits += abs(stat(perm)) >= abs(obs) - 1e-12
    tabs = []
    informative = 0
    for g in groups:
        a = sum(1 for l, c in g if c == g1 and l)
        b = sum(1 for l, c in g if c == g1 and not l)
        c_ = sum(1 for l, c in g if c == g2 and l)
        d = sum(1 for l, c in g if c == g2 and not l)
        tabs.append((c_, d, a, b))   # MH in p5 is standalone-first; g2 is standalone here
        informative += bool(a + b and c_ + d)
    return obs, (hits + 1) / (N_PERM + 1), mantel_haenszel(tabs), informative


def main():
    only23 = "--only-23" in sys.argv
    human = human_sequence()
    ref = load_structure(str(REF_PDB), chain=REF_CHAIN)
    ref_idx = {n: k for k, n in enumerate(ref.resnums)}
    missing = [p for p in POSITIONS if p not in ref_idx]
    if missing:
        sys.exit(f"positions without coordinates in 1MSV: {missing}")

    prot = {r["accession"]: r for r in read("p5_loss_vs_fusion_proteins.tsv")}
    arch = {r["accession"]: r for r in read("architecture_per_protein.tsv")}
    odc = [r["accession"] for r in read("p6_decarboxylase_partner.tsv")]
    odc_set = set(odc)
    accs = [a for a, r in prot.items() if r["class"] in ("standalone", "fused")]
    if only23:
        accs = odc
    group = {a: ("odc_fusion" if a in odc_set else prot[a]["class"]) for a in accs}

    rows = []
    for k, acc in enumerate(accs, 1):
        sites, tm = map_sites(acc, ref, ref_idx, human)
        spans = adometdc_spans(arch[acc])
        for pos in POSITIONS:
            s = sites[pos]
            rn = s["cand_resnum"]
            inside = "" if rn is None or not spans else any(a <= rn <= b for a, b in spans)
            rows.append({"accession": acc, "organism": prot[acc]["organism"], "group": group[acc],
                         "architecture": prot[acc]["architecture"], "source": prot[acc]["source"],
                         "tax_kingdom": prot[acc]["tax_kingdom"], "tax_phylum": prot[acc]["tax_phylum"] or "(none)",
                         "tm_norm_ref": round(tm, 4), "human_pos": pos, **s, "in_adometdc_pfam": inside})
        if k % 50 == 0 or k == len(accs):
            print(f"  {k}/{len(accs)}", flush=True)

    cols = ["accession", "organism", "group", "architecture", "source", "tax_kingdom", "tax_phylum", "tm_norm_ref",
            "human_pos", "human_aa", "call", "cand_resnum", "cand_aa", "ca_offset", "plddt", "state",
            "unresolved_reason", "in_adometdc_pfam"]
    out_tsv = RESULTS / ("p8_putrescine_site_23only.tsv" if only23 else "p8_putrescine_site.tsv")
    with open(out_tsv, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r[c] is None else r[c]) for c in cols})
    print(f"wrote {out_tsv}")

    # ---- per-protein view of the 23
    rep = ["Putrescine-site positions in the AdoMetDC half (p8_putrescine_site.py)",
           f"Human residues: {', '.join(f'{human[p - 1]}{p}' for p in POSITIONS)}. kept = identical & pLDDT>=70; "
           "lost = substituted & pLDDT>=70; unresolved = pLDDT<70, displaced >4 A, or nothing aligned.",
           "", "[1] The 23 AdoMetDC-PLP-decarboxylase fusions (residue / pLDDT; lower case = unresolved)"]
    by = defaultdict(dict)
    for r in rows:
        by[r["accession"]][r["human_pos"]] = r
    rep.append(f"  {'accession':<12} {'organism':<34} " + " ".join(f"{human[p - 1]}{p:<7}" for p in POSITIONS))
    for acc in odc:
        cells = []
        for p in POSITIONS:
            s = by[acc][p]
            aa = s["cand_aa"] or "-"
            aa = aa if s["state"] != "unresolved" else aa.lower()
            pl = "" if s["plddt"] is None else f"{s['plddt']:.0f}"
            cells.append(f"{aa}/{pl:<6}")
        rep.append(f"  {acc:<12} {prot[acc]['organism'][:34]:<34} " + " ".join(cells))
    outside = [r for r in rows if r["group"] == "odc_fusion" and r["in_adometdc_pfam"] is False]
    rep.append(f"  mapped residues falling outside the AdoMetDC Pfam match: {len(outside)}"
               + (" (" + ", ".join(f"{r['accession']}:{r['human_pos']}" for r in outside) + ")" if outside else ""))

    if only23:
        print("\n".join(rep))
        return

    # ---- counts and tests
    comparisons = [("odc_fusion", "standalone"), ("fused_all", "standalone"), ("fused_other", "standalone")]

    def members(g, rs):
        if g == "fused_all":
            return [dict(r, group="fused_all") for r in rs if r["group"] in ("odc_fusion", "fused")]
        if g == "fused_other":
            return [dict(r, group="fused_other") for r in rs if r["group"] == "fused"]
        return [r for r in rs if r["group"] == g]

    tests = []
    for setname, pred in (("all", lambda r: True), ("interpro_only", lambda r: r["source"] == "interpro")):
        rep.append(f"\n[2] Counts and tests, set = {setname} (kept / lost / unresolved; tests on kept vs lost)")
        for pos in [178] + [p for p in POSITIONS if p != 178]:
            rs = [r for r in rows if r["human_pos"] == pos and pred(r)]
            rep.append(f"  {human[pos - 1]}{pos}")
            cnt = {}
            for g in ("standalone", "odc_fusion", "fused_all", "fused_other"):
                c = Counter(r["state"] for r in members(g, rs))
                cnt[g] = c
                lost_aa = Counter(r["cand_aa"] for r in members(g, rs) if r["state"] == "lost")
                rep.append(f"    {g:<11} n {sum(c.values()):>3}: kept {c['kept']:>3}  lost {c['lost']:>3}  "
                           f"unresolved {c['unresolved']:>3}" + (f"   lost as {dict(lost_aa)}" if lost_aa else ""))
            for g1, g2 in comparisons:
                a, b = cnt[g1]["lost"], cnt[g1]["kept"]
                c, d = cnt[g2]["lost"], cnt[g2]["kept"]
                fp = fisher(a, b, c, d) if (a + b) and (c + d) else float("nan")
                resolved = [r for r in members(g1, rs) + members(g2, rs) if r["state"] != "unresolved"]
                if a + c == 0 or not (a + b) or not (c + d):
                    obs, pp, mh, inf = (a / (a + b) if a + b else float("nan")) - (c / (c + d) if c + d else float("nan")), float("nan"), float("nan"), 0
                    note = "no losses in either group: nothing to test" if a + c == 0 else "empty group"
                else:
                    obs, pp, mh, inf = perm_test(resolved, g1, g2)
                    note = ""
                tests.append({"set": setname, "position": f"{human[pos - 1]}{pos}", "group": g1, "vs": g2,
                              "group_kept": b, "group_lost": a, "group_unresolved": cnt[g1]["unresolved"],
                              "vs_kept": d, "vs_lost": c, "vs_unresolved": cnt[g2]["unresolved"],
                              "loss_rate_diff_pts": "" if obs != obs else round(100 * obs, 2),
                              "fisher_p": "" if fp != fp else f"{fp:.3g}",
                              "perm_within_phylum_p": "" if pp != pp else f"{pp:.3g}",
                              "mh_odds_ratio_standalone_vs_group": "" if mh != mh else f"{mh:.3g}",
                              "phyla_with_both_groups": inf, "note": note})
                rep.append(f"    {g1} vs {g2}: Fisher p = {tests[-1]['fisher_p'] or 'NA'}; within-phylum permutation "
                           f"p = {tests[-1]['perm_within_phylum_p'] or 'NA'} ({inf} phyla hold both groups)"
                           + (f"; {note}" if note else ""))
    tcols = list(tests[0].keys())
    with open(RESULTS / "p8_putrescine_site_tests.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, tcols, delimiter="\t")
        w.writeheader()
        w.writerows(tests)
    (RESULTS / "p8_putrescine_site_report.txt").write_text("\n".join(rep) + "\n", encoding="utf8")
    print("\n".join(rep))
    print(f"wrote {RESULTS / 'p8_putrescine_site_tests.tsv'}\n      {RESULTS / 'p8_putrescine_site_report.txt'}")


if __name__ == "__main__":
    sys.exit(main())
