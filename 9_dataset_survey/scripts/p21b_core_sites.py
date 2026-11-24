"""p21 B: prune the 13-site ODC intactness list to the sites an active enzyme never substitutes.

Evidence per site:
  - experimentally active ODCs in data/references/ (p6 / p10 reference set, annotation as recorded there):
    human P11926, T. brucei P07805, yeast P08432, tobacco P93351, Plasmodium falciparum Q8IJ77 (bifunctional,
    ODC activity Muller 2000), Selenomonas O50657 (LDC/ODC, EC 4.1.1.17 + 4.1.1.18). Residue at each site by local
    alignment to P11926 (BLOSUM62, 11/1, as p10).
  - Q8IJ77 also by structure (p10_all_fusions.tsv ODC_* columns, p10 gate).
  - all 65 Plasmodium rows of p10_all_fusions.tsv (sequence read, column sequence_sites; structure where present).
  - context only, not used for the core: active fold-III decarboxylases of other substrates in the same reference
    set (PBCV-1 ADC, Arabidopsis ADC1/2, E. coli SpeA, DAPDC LysA / DAPDC1, CANSDC x2, LDC2).
CORE = sites where every active ODC above carries the human residue.
Re-score: the 23 structure-verified fusions from results/p9_odc_sites_gated.tsv (read only), with p9's rule applied
to each list: degraded = >= 1 site "lost"; intact = every listed site "kept"; otherwise undetermined.

Writes results/p21b_core_sites.tsv (one row per site), results/p21b_rescore.tsv (one row per fusion),
results/p21b_report.md
"""
import csv
import sys
from collections import Counter, defaultdict

import p6_decarboxylase_partner as p6
import p10_red_algal_pattern as p10
from p2_fetch_annotations import RESULTS
from p9_self_processing import ODC_SITE_LABELS

sys.stdout.reconfigure(encoding="utf8")
ACTIVE_ODC = ["P11926", "P07805", "P08432", "P93351", "Q8IJ77", "O50657"]
OTHER_ACTIVE = ["Q84527", "Q9SI64", "O23141", "P21170", "P00861", "Q949X7", "Q56575", "A8FNH9", "A0A8F1NN86"]


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sites_by_seq(seq, odc_seq):
    a = p10.LOCAL.align(seq, odc_seq)[0]
    res = {j + 1: seq[i] for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))}
    return {l: res.get(int(l[1:]), "-") for l in ODC_SITE_LABELS}


def parse_seq_sites(s):
    out = {}
    for t in s.split():
        for l in ODC_SITE_LABELS:
            if t.startswith(l):
                out[l] = l[0] if t[len(l):] == "=" else t[len(l) + 1:]
    return out


def status(states):
    c = Counter(states)
    return "degraded" if c["lost"] else "intact" if c["kept"] == len(states) else "undetermined"


def main():
    refs = p10.references()
    labels = {**p6.REFS, **p10.EXTRA_REFS}
    odc_seq = refs["P11926"]
    ref_sites = {a: sites_by_seq(refs[a], odc_seq) for a in ACTIVE_ODC + OTHER_ACTIVE}
    allf = read("p10_all_fusions.tsv")
    q = next(r for r in allf if r["accession"] == "Q8IJ77")
    q_struct = {l: q[f"ODC_{l}"] for l in ODC_SITE_LABELS}
    plas = [r for r in allf if r["organism"].startswith("Plasmodium")]
    plas_species = {" ".join(r["organism"].split()[:2]) for r in plas}
    plas_seq = [parse_seq_sites(r["sequence_sites"]) for r in plas]

    rows, core = [], []
    for l in ODC_SITE_LABELS:
        h = l[0]
        act = {a: ref_sites[a][l] for a in ACTIVE_ODC}
        kept_all = all(v == h for v in act.values())
        pk = Counter(s.get(l, "-") for s in plas_seq)
        oth = {a: ref_sites[a][l] for a in OTHER_ACTIVE}
        if kept_all:
            core.append(l)
            why = "human residue in all 6 active ODCs"
        else:
            bad = ", ".join(f"{labels[a][0]} {v}" for a, v in act.items() if v != h)
            why = f"substituted in an active ODC: {bad}"
        rows.append({"site": l, "core": kept_all, "reason": why,
                     **{f"active:{labels[a][0]}": v for a, v in act.items()},
                     "Q8IJ77_structure": q_struct[l],
                     "plasmodium_65_human_residue": f"{pk.get(h, 0)}/{len(plas)}",
                     "plasmodium_65_other": "; ".join(f"{k} {v}" for k, v in pk.most_common() if k != h),
                     "other_active_decarboxylases_human_residue": f"{sum(v == h for v in oth.values())}/{len(oth)}",
                     "other_active_detail": " ".join(f"{labels[a][0].split()[0]}:{v}" for a, v in oth.items())})
    with open(RESULTS / "p21b_core_sites.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    st = defaultdict(dict)
    org = {}
    for r in read("p9_odc_sites_gated.tsv"):
        st[r["accession"]][r["odc_site"]] = (r["state"], r["cand_aa"])
        org[r["accession"]] = r["organism"]
    res = []
    for a, d in st.items():
        s13 = status([d[l][0] for l in ODC_SITE_LABELS])
        sc = status([d[l][0] for l in core])
        res.append({"accession": a, "organism": org[a],
                    "kept_lost_unres_13": "/".join(str(Counter(d[l][0] for l in ODC_SITE_LABELS)[k]) for k in ("kept", "lost", "unresolved")),
                    "call_13": s13,
                    "kept_lost_unres_core": "/".join(str(Counter(d[l][0] for l in core)[k]) for k in ("kept", "lost", "unresolved")),
                    "call_core": sc, "changed": s13 != sc,
                    "lost_sites": " ".join(f"{l}{d[l][1]}" for l in ODC_SITE_LABELS if d[l][0] == "lost"),
                    "unresolved_sites": " ".join(l for l in ODC_SITE_LABELS if d[l][0] == "unresolved")})
    with open(RESULTS / "p21b_rescore.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(res[0]), delimiter="\t")
        w.writeheader()
        w.writerows(res)

    drop = [r for r in rows if not r["core"]]
    L = ["# p21 B: a pruned core of the 13 ODC sites", "",
         "Generated by `scripts/p21b_core_sites.py`.", "",
         f"**Core list ({len(core)} sites): {' '.join(core)}.** Dropped: "
         + "; ".join(f"{r['site']} ({r['reason'].split(': ', 1)[1]})" for r in drop) + ".", "",
         "Rule: a site is core only if all six experimentally active ODCs in `data/references/` carry the human "
         "residue there (human, *T. brucei*, yeast, tobacco, *Plasmodium falciparum* Q8IJ77, *Selenomonas* LDC/ODC). "
         "Residues are read by sequence alignment to human ODC. The 65 *Plasmodium* fusions and nine active "
         "decarboxylases of other substrates are shown for context and do not define the core.", "",
         f"*Plasmodium* in p10: {len(plas)} rows from {len(plas_species)} species names.", "",
         "| site | core | " + " | ".join(labels[a][0] for a in ACTIVE_ODC) + " | Q8IJ77 structure | Plasmodium 65: human residue | other residues | other active decarboxylases: human residue |",
         "|---|---|" + "---|" * len(ACTIVE_ODC) + "---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['site']} | {'yes' if r['core'] else '**no**'} | "
                 + " | ".join(r[f"active:{labels[a][0]}"] for a in ACTIVE_ODC)
                 + f" | {r['Q8IJ77_structure']} | {r['plasmodium_65_human_residue']} | {r['plasmodium_65_other'] or '-'} | "
                   f"{r['other_active_decarboxylases_human_residue']} |")
    ch = [r for r in res if r["changed"]]
    L += ["", "## Re-scoring the 23 fusions (p9 gated states, p9 rule applied to each list)", "",
          f"**Calls that change: {len(ch)}.**", ""]
    L += [f"- {r['accession']} *{r['organism']}*: {r['call_13']} -> **{r['call_core']}** (13 sites k/l/u "
          f"{r['kept_lost_unres_13']}, core {r['kept_lost_unres_core']}; lost {r['lost_sites'] or 'none'}; "
          f"unresolved {r['unresolved_sites'] or 'none'})" for r in ch] or ["- none"]
    L += ["", "| accession | organism | 13 sites k/l/u | call | core k/l/u | call | lost sites |", "|---|---|---|---|---|---|---|"]
    for r in sorted(res, key=lambda r: (r["call_core"], r["organism"])):
        L.append(f"| {r['accession']} | *{r['organism']}* | {r['kept_lost_unres_13']} | {r['call_13']} | "
                 f"{r['kept_lost_unres_core']} | {r['call_core']} | {r['lost_sites'] or '-'} |")
    # every p10 fusion with a structural read (includes Q8IJ77 and the other non-23 proteins)
    p10ch, p10rows = [], []
    for r in allf:
        if not r["odc_status"] or r["set"] == "our 23":
            continue
        s = {l: r[f"ODC_{l}"].split()[-1] for l in ODC_SITE_LABELS}
        c13, cc = status(list(s.values())), status([s[l] for l in core])
        p10rows.append({"accession": r["accession"], "organism": r["organism"], "group": r["group"],
                        "p10_odc_status": r["odc_status"], "call_13": c13, "call_core": cc,
                        "lost": " ".join(f"{l}{r[f'ODC_{l}'].split()[0]}" for l in ODC_SITE_LABELS if s[l] == "lost")})
        if c13 != cc:
            p10ch.append(p10rows[-1])
    with open(RESULTS / "p21b_rescore_p10.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(p10rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(p10rows)
    L += ["", "## The p10 fusions outside the 23 (structural reads from p10_all_fusions.tsv)", "",
          f"{len(p10rows)} have a structural read. Calls that change between the 13-site and the core list: "
          f"**{len(p10ch)}**.", ""]
    L += [f"- {r['accession']} *{r['organism']}* ({r['group']}): {r['call_13']} -> **{r['call_core']}** "
          f"(lost {r['lost'] or 'none'})" for r in p10ch] or ["- none"]
    L += ["", "Note: p10's own `odc_status` also requires a partner of >= 300 aa and no unaligned stretch of >= 100 "
          "residues (\"incomplete\"); the calls here use the site rule only. All rows are in p21b_rescore_p10.tsv.",
          "", "Correction to p20: the 65 *Plasmodium* rows carry Phe at Y331 in 62, not 60 (p20 counted only rows "
          "whose sole non-identical site was Y331).", ""]
    (RESULTS / "p21b_report.md").write_text("\n".join(L) + "\n", encoding="utf8")
    print("core:", core)
    print("changed:", [(r["accession"], r["call_13"], r["call_core"]) for r in ch])


if __name__ == "__main__":
    main()
