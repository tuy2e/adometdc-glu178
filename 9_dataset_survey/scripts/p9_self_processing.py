"""Track B step 2: does the AdoMetDC half of the 23 fusions still self-process, and is the ODC half intact?

Sites: Marc's five prespecified AdoMetDC sites, read from his sites/adometdc.tsv (never written):
    Ser68  pyruvoyl self-cleavage        Ser229, His243  processing proton acceptors
    Cys82  catalytic proton donor        Glu247          substrate positioning
1MSV is the S68A proenzyme mutant, so position 68 is compared with Ser (P17707), not with the Ala in 1MSV.

Mapping and states: exactly p8 (p8_putrescine_site.map_sites: p3's TM-align onto 1MSV chain A, cached in
data/p3_alignments/), plus one state from Marc's sites file:
    kept        identical, pLDDT >= 70
    altered     substituted by a residue listed in the file's `altered` column (S229T, E247Q), pLDDT >= 70
    lost        any other substitution, pLDDT >= 70
    unresolved  pLDDT < 70, displaced > 4 A, or nothing aligned
Tests per site, 23 ODC fusions vs 590 standalone, lost vs not lost (kept + altered), unresolved left out:
Fisher exact, within-phylum permutation (p8/p5 method), MH odds ratio; repeated for the InterPro-only set.

Both halves, one row per fusion (the paper-figure table):
    AdoMetDC half  intact        all five kept or altered
                   not intact    at least one lost
                   undetermined  none lost, at least one unresolved
    ODC half       intact        13/13 in p7_partner_integrity.tsv (functional_residues_identical_in_place)
                   not intact    fewer than 13
    p7 has no per-residue pLDDT gate (identical residue within 4 A of human ODC's), so the two halves are
    not judged by quite the same rule; the table carries both counts.

Writes results/p9_self_processing.tsv (protein x site), p9_self_processing_tests.tsv,
       p9_both_halves.tsv, p9_self_processing.md
Run:  py -3.12 p9_self_processing.py            (reuses the p8 alignment cache)
      py -3.12 p9_self_processing.py --only-23  (the 23 and the both-halves table only, no tests)
"""
import csv
import re
import sys
from collections import Counter

from foldtrace.io import load_structure

import p8_putrescine_site as p8
from p3_residue_grid import RUN, REF_PDB, REF_CHAIN, PLDDT_OK, human_sequence, plddt_by_resnum
from p5_loss_vs_fusion import RESULTS, fisher, read

sys.stdout.reconfigure(encoding="utf8")
ODC_SITE_LABELS = ["K69", "H197", "S200", "G237", "E274", "P275", "G276", "R277", "Y331", "D332", "C360", "D361",
                   "Y389"]   # p7's columns (human ODC P11926 numbering)
BORDERLINE = 1.5   # pLDDT points either side of the 70 gate
ODC_GROUPS = [[69], [197, 200], [237], [274, 275, 276, 277], [331, 332], [360, 361], [389]]
SITES_FILE = RUN.parent.parent / "sites" / "adometdc.tsv"


def marc_sites():
    with open(SITES_FILE, encoding="utf8") as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    return [(r["label"], int(r["ref_resnum"]), r["expected"], r["role"], set(r["altered"].split(",")) - {""})
            for r in rows]


def main():
    only23 = "--only-23" in sys.argv
    sites = marc_sites()
    p8.POSITIONS = [s[1] for s in sites]
    human = human_sequence()
    for lab, pos, exp, _, _ in sites:
        if human[pos - 1] != exp:
            sys.exit(f"{lab}: P17707 has {human[pos - 1]} at {pos}, sites file expects {exp}")
    ref = load_structure(str(REF_PDB), chain=REF_CHAIN)
    ref_idx = {n: k for k, n in enumerate(ref.resnums)}

    prot = {r["accession"]: r for r in read("p5_loss_vs_fusion_proteins.tsv")}
    arch = {r["accession"]: r for r in read("architecture_per_protein.tsv")}
    odc = [r["accession"] for r in read("p6_decarboxylase_partner.tsv")]
    odc_set = set(odc)
    accs = odc if only23 else [a for a, r in prot.items() if a in odc_set or r["class"] == "standalone"]

    rows = []
    for k, acc in enumerate(accs, 1):
        mapped, tm = p8.map_sites(acc, ref, ref_idx, human)
        spans = p8.adometdc_spans(arch[acc])
        for lab, pos, exp, role, altered in sites:
            s = dict(mapped[pos])
            if s["state"] == "lost" and s["cand_aa"] in altered:
                s["state"] = "altered"
            rn = s["cand_resnum"]
            rows.append({"accession": acc, "organism": prot[acc]["organism"],
                         "group": "odc_fusion" if acc in odc_set else "standalone",
                         "source": prot[acc]["source"], "tax_kingdom": prot[acc]["tax_kingdom"],
                         "tax_phylum": prot[acc]["tax_phylum"] or "(none)", "tm_norm_ref": round(tm, 4),
                         "site": lab, "role": role, "human_pos": pos, **s,
                         "in_adometdc_pfam": "" if rn is None or not spans else any(a <= rn <= b for a, b in spans)})
        if k % 50 == 0 or k == len(accs):
            print(f"  {k}/{len(accs)}", flush=True)

    cols = ["accession", "organism", "group", "source", "tax_kingdom", "tax_phylum", "tm_norm_ref", "site", "role",
            "human_pos", "human_aa", "call", "cand_resnum", "cand_aa", "ca_offset", "plddt", "state",
            "unresolved_reason", "in_adometdc_pfam"]
    with open(RESULTS / "p9_self_processing.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r[c] is None else r[c]) for c in cols})

    # ---- ODC half with the same pLDDT gate (p7's mapping, read from p7_partner_integrity.tsv)
    integ = {r["accession"]: r for r in read("p7_partner_integrity.tsv")}
    odc_rows = []
    for acc in odc:
        pl = plddt_by_resnum(RUN / "models" / f"{acc}.pdb")
        for lab in ODC_SITE_LABELS:
            cellv = integ[acc][lab]
            m = re.match(r"([A-Z])(\d+) \(([\d.]+)A\)", cellv)
            aa, rn, off = (m.group(1), int(m.group(2)), float(m.group(3))) if m else (None, None, None)
            p = pl.get(rn) if rn is not None else None
            if aa is None:
                call, state, why = "no_residue_aligned", "unresolved", "no_residue_aligned"
            else:
                call = "displaced" if off > p8.OFFSET else ("identical" if aa == lab[0] else "substituted")
                if call == "displaced":
                    state, why = "unresolved", "displaced" + (";plddt_below_70" if p < PLDDT_OK else "")
                elif p < PLDDT_OK:
                    state, why = "unresolved", "plddt_below_70"
                else:
                    state, why = ("kept" if call == "identical" else "lost"), ""
            odc_rows.append({"accession": acc, "organism": prot[acc]["organism"], "odc_site": lab,
                             "human_odc_pos": int(lab[1:]), "human_aa": lab[0], "call": call, "cand_resnum": rn,
                             "cand_aa": aa, "ca_offset": off, "plddt": p, "state": state, "unresolved_reason": why})
    with open(RESULTS / "p9_odc_sites_gated.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(odc_rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in odc_rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})

    def half(st):
        return "not intact" if st["lost"] else ("undetermined" if st["unresolved"] else "intact")

    def combine(ado, odc_ok):
        return {("intact", "intact"): "both intact", ("intact", "not intact"): "AdoMetDC half only",
                ("not intact", "intact"): "ODC half only", ("not intact", "not intact"): "neither",
                }.get((ado, odc_ok), f"AdoMetDC {ado}; ODC {odc_ok}")

    def is_borderline(r):
        """The gate alone decides the call: residue in place (identical/substituted, <= 4 A) with pLDDT
        within BORDERLINE of 70."""
        return r["call"] in ("identical", "substituted") and r["plddt"] is not None \
            and abs(r["plddt"] - PLDDT_OK) < BORDERLINE

    def flipped(r, altered=()):
        """State the site would get on the other side of the gate."""
        if r["state"] != "unresolved":
            return "unresolved"
        return "kept" if r["call"] == "identical" else ("altered" if r["cand_aa"] in altered else "lost")

    # ---- both halves, one row per fusion
    labels = [s[0] for s in sites]
    alt_of = {s[0]: s[4] for s in sites}
    both = []
    for acc in odc:
        rs = {r["site"]: r for r in rows if r["accession"] == acc}
        st = Counter(r["state"] for r in rs.values())
        ado = half(st)
        n13 = integ[acc]["functional_residues_identical_in_place"]
        odc_p7 = "intact" if n13.split("/")[0] == n13.split("/")[1] else "not intact"
        ors = [r for r in odc_rows if r["accession"] == acc]
        ost = Counter(r["state"] for r in ors)
        odc_ok = half(ost)
        verdict = combine(ado, odc_ok)

        # borderline sites, each flipped on its own: does the per-fusion verdict change?
        border = []
        for side, r, lab in [("AdoMetDC", r, r["site"]) for r in rs.values()] + [("ODC", r, r["odc_site"]) for r in ors]:
            if not is_borderline(r):
                continue
            new = flipped(r, alt_of.get(lab, ()))
            if side == "AdoMetDC":
                alt = combine(half(st - Counter([r["state"]]) + Counter([new])), odc_ok)
            else:
                alt = combine(ado, half(ost - Counter([r["state"]]) + Counter([new])))
            border.append(f"{side} {lab} {r['cand_aa']} pLDDT {r['plddt']:.1f}: {r['state']} -> {new}; "
                          + (f"verdict would be '{alt}'" if alt != verdict else "verdict unchanged"))

        def cell(r):
            aa = r["cand_aa"] or "-"
            pl = "" if r["plddt"] is None else f" ({r['plddt']:.0f}{'†' if is_borderline(r) else ''})"
            return f"{aa}{pl}" if r["state"] in ("kept", "altered") else f"{aa}{pl} {r['state']}"
        both.append({"accession": acc, "organism": prot[acc]["organism"],
                     **{lab: cell(rs[lab]) for lab in labels},
                     "adometdc_sites_kept_or_altered": f"{st['kept'] + st['altered']}/{len(sites)}",
                     "adometdc_half": ado, "odc_p7_identical_in_place": n13, "odc_half_p7_ungated": odc_p7,
                     "odc_kept": ost["kept"], "odc_lost": ost["lost"], "odc_unresolved": ost["unresolved"],
                     "odc_TM_to_human_ODC": integ[acc]["TM_to_human_ODC"],
                     "odc_partner_len": integ[acc]["partner_len"], "odc_half": odc_ok, "both_halves": verdict,
                     "borderline_sites": " | ".join(border),
                     "verdict_depends_on_borderline": "yes" if any("would be" in b for b in border) else ""})
    with open(RESULTS / "p9_both_halves.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(both[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(both)

    md = ["# Track B step 2: AdoMetDC self-processing sites and the ODC half (p9)", "",
          "Generated by `scripts/p9_self_processing.py`. Sites are 📋 Marc's five prespecified AdoMetDC sites "
          "(`aim1_workflow/sites/adometdc.tsv`), mapped by the p8 method (TM-align onto 1MSV chain A, human P17707 "
          "numbering). kept = identical, pLDDT >= 70; altered = the substitution Marc's file lists as tolerated "
          "(S229T, E247Q), pLDDT >= 70; lost = other substitution, pLDDT >= 70; unresolved = pLDDT < 70, displaced "
          "> 4 A, or nothing aligned.", "",
          "ODC half: the 13 human-ODC functional residues as mapped by p7 (TM-align of the partner region onto "
          "the AlphaFold model of human ODC; residue and CA offset read from `p7_partner_integrity.tsv`), now "
          "with the same gate as the AdoMetDC half: kept = identical, <= 4 A, pLDDT >= 70; lost = substituted, "
          "<= 4 A, pLDDT >= 70; unresolved = pLDDT < 70, displaced > 4 A, or nothing aligned. Per residue in "
          "`p9_odc_sites_gated.tsv`. A half is intact when every site is kept (or altered), not intact when "
          "any site is lost, undetermined when none is lost but some are unresolved. p7's ungated count "
          "(identical within 4 A, any pLDDT) is kept alongside.", "",
          "## Both halves, per fusion", "",
          "Cells: residue (pLDDT); the state is written out when it is not kept/altered. ODC column: "
          "kept / lost / unresolved of 13 (p7 ungated count in brackets).", "",
          f"Borderline: a site whose residue is in place with pLDDT within {BORDERLINE} of 70, so the cutoff alone "
          "decides its call. † marks it (AdoMetDC cells: after the pLDDT; ODC column: after the counts, sites "
          "listed below the table). ‡ marks a verdict that changes if any one borderline site falls on the other "
          "side of 70.", "",
          "| accession | organism | " + " | ".join(labels) + " | AdoMetDC half | ODC k/l/u (p7) | ODC half | both halves |",
          "|---|---|" + "---|" * len(labels) + "---|---|---|---|"]
    for b in both:
        odc_border = "†" if "ODC " in b["borderline_sites"] else ""
        md.append(f"| {b['accession']} | *{b['organism']}* | " + " | ".join(b[l] for l in labels)
                  + f" | {b['adometdc_half']} | {b['odc_kept']}/{b['odc_lost']}/{b['odc_unresolved']}{odc_border} "
                  f"({b['odc_p7_identical_in_place']}) | {b['odc_half']} | {b['both_halves']}"
                  + ("‡" if b["verdict_depends_on_borderline"] else "") + " |")
    vc = Counter(b["both_halves"] for b in both)
    md += ["", "Borderline sites:", ""]
    md += [f"- {b['accession']} *{b['organism']}*: " + s for b in both if b["borderline_sites"]
           for s in b["borderline_sites"].split(" | ")]
    md += ["", "Totals: " + "; ".join(f"{k} {v}" for k, v in vc.most_common()) + ".",
           "ODC half: gated " + "; ".join(f"{k} {v}" for k, v in Counter(b["odc_half"] for b in both).most_common())
           + " | p7 ungated " + "; ".join(f"{k} {v}" for k, v in
                                          Counter(b["odc_half_p7_ungated"] for b in both).most_common()) + ".", ""]

    # ---- the gutted ODC halves (p7 count <= 5/13): where the losses sit and how confident the model is there
    md += ["## Gutted ODC halves (p7: at most 5 of 13 identical in place)", "",
           "Site groups along human ODC: " + " | ".join(",".join(str(p) for p in g) for g in ODC_GROUPS)
           + f" ({len(ODC_GROUPS)} groups; sites within 5 residues of each other share a group). Partner-domain "
           "pLDDT is over every model residue from p7's partner start; low stretches are runs of >= 5 "
           "consecutive residues below 70.", "",
           "| accession | organism | lost (residue, pLDDT) | unresolved (reason) | groups with a loss | "
           "lost-residue pLDDT min-max | partner pLDDT mean, % < 70 | low stretches (model numbering) |", "|---|---|---|---|---|---|---|---|"]
    for acc in odc:
        n7 = int(integ[acc]["functional_residues_identical_in_place"].split("/")[0])
        if n7 > 5:
            continue
        rs = [r for r in odc_rows if r["accession"] == acc]
        lost = [r for r in rs if r["state"] == "lost"]
        unr = [r for r in rs if r["state"] == "unresolved"]
        pl = plddt_by_resnum(RUN / "models" / f"{acc}.pdb")
        lo = int(integ[acc]["partner_from"])
        part = sorted((n, v) for n, v in pl.items() if n >= lo)
        runs, cur = [], []
        for n, v in part:
            if v < PLDDT_OK and (not cur or n == cur[-1] + 1):
                cur.append(n)
            else:
                if len(cur) >= 5:
                    runs.append((cur[0], cur[-1]))
                cur = [n] if v < PLDDT_OK else []
        if len(cur) >= 5:
            runs.append((cur[0], cur[-1]))
        grp = {k for k, g in enumerate(ODC_GROUPS) for r in lost if r["human_odc_pos"] in g}
        vals = [v for _, v in part]
        lp = [r["plddt"] for r in lost]
        md.append(f"| {acc} | *{prot[acc]['organism']}* | "
                  + (", ".join(f"{r['odc_site']}→{r['cand_aa']} ({r['plddt']:.0f})" for r in lost) or "none")
                  + " | " + (", ".join(f"{r['odc_site']} ({r['unresolved_reason']}"
                                       + ("" if r["plddt"] is None else f", {r['plddt']:.0f}") + ")" for r in unr) or "none")
                  + f" | {len(grp)}/{len(ODC_GROUPS)} | " + (f"{min(lp):.0f}-{max(lp):.0f}" if lp else "")
                  + f" | {sum(vals) / len(vals):.1f}, {100 * sum(v < PLDDT_OK for v in vals) / len(vals):.0f}% | "
                  + (", ".join(f"{a}-{b}" for a, b in runs) or "none") + " |")
    md.append("")
    outside = [r for r in rows if r["group"] == "odc_fusion" and r["in_adometdc_pfam"] is False]
    md.append(f"Mapped AdoMetDC-site residues falling outside the AdoMetDC Pfam match: {len(outside)}"
              + (" (" + ", ".join(f"{r['accession']} {r['site']}" for r in outside) + ")." if outside else "."))

    if only23:
        (RESULTS / "p9_self_processing.md").write_text("\n".join(md) + "\n", encoding="utf8")
        print("\n".join(md))
        return

    # ---- counts and tests, 23 vs standalone
    tests = []
    for setname, pred in (("all", lambda r: True), ("interpro_only", lambda r: r["source"] == "interpro")):
        md += ["", f"## Per site, 23 ODC fusions vs standalone (set: {setname})", "",
               "Tests compare lost with not lost (kept + altered); unresolved is left out.", "",
               "| site | group | n | kept | altered | lost | unresolved | Fisher p | within-phylum perm p | phyla with both |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for lab, pos, *_ in sites:
            rs = [r for r in rows if r["site"] == lab and pred(r)]
            cnt = {g: Counter(r["state"] for r in rs if r["group"] == g) for g in ("odc_fusion", "standalone")}
            a, b = cnt["odc_fusion"]["lost"], cnt["odc_fusion"]["kept"] + cnt["odc_fusion"]["altered"]
            c, d = cnt["standalone"]["lost"], cnt["standalone"]["kept"] + cnt["standalone"]["altered"]
            fp = fisher(a, b, c, d) if (a + b) and (c + d) else float("nan")
            resolved = [dict(r, state="lost" if r["state"] == "lost" else "kept")
                        for r in rs if r["state"] != "unresolved"]
            if a + c == 0 or not (a + b) or not (c + d):
                obs, pp, mh, inf = float("nan"), float("nan"), float("nan"), 0
                note = "no losses in either group: nothing to test" if a + c == 0 else "empty group"
            else:
                obs, pp, mh, inf = p8.perm_test(resolved, "odc_fusion", "standalone")
                note = ""
            lost_as = lambda g: dict(Counter(r["cand_aa"] for r in rs if r["group"] == g and r["state"] == "lost"))
            t = {"set": setname, "site": lab, "group": "odc_fusion", "vs": "standalone",
                 **{f"group_{s}": cnt["odc_fusion"][s] for s in ("kept", "altered", "lost", "unresolved")},
                 **{f"vs_{s}": cnt["standalone"][s] for s in ("kept", "altered", "lost", "unresolved")},
                 "group_lost_as": lost_as("odc_fusion") or "", "vs_lost_as": lost_as("standalone") or "",
                 "loss_rate_diff_pts": "" if obs != obs else round(100 * obs, 2),
                 "fisher_p": "" if fp != fp else f"{fp:.3g}", "perm_within_phylum_p": "" if pp != pp else f"{pp:.3g}",
                 "mh_odds_ratio_standalone_vs_group": "" if mh != mh else f"{mh:.3g}",
                 "phyla_with_both_groups": inf, "note": note}
            tests.append(t)
            for g in ("odc_fusion", "standalone"):
                k_ = cnt[g]
                first = g == "odc_fusion"
                md.append(f"| {lab if first else ''} | {g} | {sum(k_.values())} | {k_['kept']} | {k_['altered']} | "
                          f"{k_['lost']} | {k_['unresolved']} | " + (f"{t['fisher_p'] or 'NA'} | "
                          f"{t['perm_within_phylum_p'] or 'NA'} | {inf} |" if first else " | | |"))
        md.append("")
        md.append("Notes: " + "; ".join(f"{t['site']}: {t['note']}" for t in tests if t["set"] == setname and t["note"])
                  if any(t["note"] for t in tests if t["set"] == setname) else "")
    with open(RESULTS / "p9_self_processing_tests.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(tests[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(tests)
    (RESULTS / "p9_self_processing.md").write_text("\n".join(md) + "\n", encoding="utf8")
    print("\n".join(md))
    print(f"wrote {RESULTS / 'p9_self_processing.tsv'}, p9_self_processing_tests.tsv, p9_both_halves.tsv, "
          "p9_self_processing.md")


if __name__ == "__main__":
    sys.exit(main())
