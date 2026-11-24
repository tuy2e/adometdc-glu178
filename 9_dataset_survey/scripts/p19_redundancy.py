"""Track B p19: does the fused ODC half stay intact only where it is the organism's ONLY ODC?

Species = the organisms of the 23 structure-verified fusions (results/p9_both_halves.tsv). Taxids from
results/p10_fusion_search.tsv, except Closterium sp. NIES-68 (2019903; p10 recorded the parent complex 2997407).
Check: each species' census must contain its own fusion(s), otherwise the taxid is wrong.

1. Fusion ODC half, from results/p9_odc_sites_gated.tsv (read only): degraded = >= 1 site "lost";
   intact = 13/13 "kept"; otherwise undetermined. Species with several fusions: intact if any fusion is intact,
   degraded if any is degraded and none intact, else undetermined.
2. Census, the UniParc route of p10 (not active UniProt): every UniParc sequence of the taxon with a
   PLP-decarboxylase (PF02784/PF00278) or AdoMetDC (PF01536/PF02675) feature; cached in data/p19/.
   Standalone = PLP feature, no AdoMetDC feature.
3. ODC family = identity to human ODC P11926 >= ODC_ID (local alignment, BLOSUM62, gap 11/1, as p10/p17).
   The threshold sits in the gap p17 saw (25-31 % DAPDC-like vs 42-48 % ODC-like); every identity is listed.
4. 13 human-ODC sites per ODC-family protein: structural (p10 odc_structural: TM-align to human ODC, CA within
   4 A, pLDDT >= 70) when an AlphaFold model already exists in data/p10_models/ (nothing new is downloaded),
   otherwise the sequence alignment. Intact standalone ODC = ODC family and >= 11 of 13 kept.
5. Proteome (the confounder): UniProt proteomes of the taxon (cached in data/p19/). Genome-derived = proteome
   built on a GCA_/GCF_ genome assembly. Complete = BUSCO complete >= BUSCO_MIN of the lineage set.
   A species with no intact standalone ODC is "cannot tell" unless it has a complete genome-derived proteome:
   absence there cannot be told from missing data, and it is excluded from the tests.
6. Fisher exact (two-sided, p5): fused ODC half intact / degraded x intact standalone ODC yes / no.
   Test A: all species with a decided fusion status and a decided standalone status.
   Test B: only species with a complete genome-derived proteome.

Writes results/p19_redundancy.tsv (one row per species), results/p19_proteins.tsv (every census protein),
results/p19_report.md
Run:  py -3.12 p19_redundancy.py
"""
import csv
import json
import sys
from collections import Counter, defaultdict

import p10_red_algal_pattern as p10
import p7_fusion_validation as p7
from foldtrace.io import load_structure
from p2_fetch_annotations import DATA, RESULTS, get
from p5_loss_vs_fusion import fisher
from p6_decarboxylase_partner import references
from p9_self_processing import ODC_SITE_LABELS

sys.stdout.reconfigure(encoding="utf8")
P19 = DATA / "p19"
P19.mkdir(parents=True, exist_ok=True)
p10.CACHE = P19            # p10's cached() / uniparc_search() write here, not into data/p10/
ODC_ID, BUSCO_MIN, INTACT_MIN = 0.365, 0.80, 11
TAXID_FIX = {"Closterium sp. NIES-68": 2019903}
ADO, PLP = {"PF01536", "PF02675"}, {"PF02784", "PF00278"}


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fusion_status():
    st = defaultdict(Counter)
    for r in read("p9_odc_sites_gated.tsv"):
        st[r["accession"]][r["state"]] += 1
    return {a: ("degraded" if c["lost"] else "intact" if c["kept"] == 13 else "undetermined") for a, c in st.items()}


def proteomes(taxid):
    p = P19 / f"proteomes_{taxid}.json"
    if not p.exists():
        st, _, body = get(f"https://rest.uniprot.org/proteomes/search?query=organism_id:{taxid}&format=json&size=50")
        p.write_text(body, encoding="utf8")
    out = []
    for r in json.loads(p.read_text(encoding="utf8")).get("results", []):
        ga = r.get("genomeAssembly") or {}
        b = (r.get("proteomeCompletenessReport") or {}).get("buscoReport") or {}
        frac = b["complete"] / b["total"] if b.get("total") else None
        out.append({"id": r.get("id"), "type": r.get("proteomeType"), "proteins": r.get("proteinCount"),
                    "assembly": ga.get("assemblyId") or "", "genome": (ga.get("assemblyId") or "").startswith(("GCA_", "GCF_")),
                    "busco": frac, "busco_txt": f"{b['complete']}/{b['total']} ({b.get('lineageDb', '').split()[0]})" if frac else "n/a"})
    return out


def best_proteome(ps):
    if not ps:
        return None
    return max(ps, key=lambda p: (p["genome"], p["busco"] or 0, p["proteins"] or 0))


def seq_identity_and_sites(seq, odc_seq):
    a = p10.LOCAL.align(seq, odc_seq)[0]
    pairs = [(i, j) for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))]
    ident = sum(seq[i] == odc_seq[j] for i, j in pairs) / len(pairs) if pairs else 0.0
    res = {j + 1: seq[i] for i, j in pairs}
    kept = sum(res.get(int(l[1:])) == l[0] for l in ODC_SITE_LABELS)
    return ident, len(pairs), kept, " ".join(f"{l}{'=' if res.get(int(l[1:])) == l[0] else ':' + res.get(int(l[1:]), '-')}"
                                             for l in ODC_SITE_LABELS)


def main():
    odc_seq = references()["P11926"]
    odc = p7.sub_structure(load_structure(str(p7.af_model("P11926"))), 1, p7.ODC_ORDERED)
    fstat = fusion_status()
    tax = {}
    for r in read("p10_fusion_search.tsv"):
        if r["in_our_23"]:
            tax[r["in_our_23"]] = int(r["taxid"])
    fus_by_sp = defaultdict(list)
    for r in read("p9_both_halves.tsv"):
        fus_by_sp[r["organism"]].append(r["accession"])
    q = " OR ".join(f"feature_id:{x}" for x in sorted(PLP | ADO))
    prot_rows, sp_rows = [], []
    for sp, accs in fus_by_sp.items():
        taxid = TAXID_FIX.get(sp, tax[accs[0]])
        ents = p10.uniparc_search(f"organism_id:{taxid} AND ({q})")
        found = set()
        standalone = []
        for e in ents:
            ids = {s[0] for s in p10.pfam_spans(e)}
            ups = [a.split(".")[0] for a in (e.get("uniProtKBAccessions") or [])]
            if set(ups) & set(accs):
                found |= set(ups) & set(accs)
            if not ids & PLP or ids & ADO:
                continue
            seq = e["sequence"]["value"]
            ident, n, kept_seq, sites_seq = seq_identity_and_sites(seq, odc_seq)
            fam = ident >= ODC_ID
            row = {"species": sp, "taxid": taxid, "upi": e["uniParcId"], "uniprot": ",".join(ups), "length": len(seq),
                   "pfam": " > ".join(f"{a} {b}-{c}" for a, b, c in p10.pfam_spans(e)),
                   "identity_to_P11926": round(ident, 3), "aligned_residues": n,
                   "odc_family": fam, "sites_method": "", "sites_kept": "", "sites": "", "intact_odc": False}
            if fam:
                model = next((p10.MODELS / f"AF-{u}.pdb" for u in ups if (p10.MODELS / f"AF-{u}.pdb").exists()), None)
                if model:
                    o = p10.odc_structural(odc, model, 1)
                    c = Counter(v["state"] for v in o["sites"].values())
                    row.update({"sites_method": f"structure ({model.stem[3:]}, p10 gate)", "sites_kept": c["kept"],
                                "sites": f"{c['kept']} kept / {c['lost']} lost / {c['unresolved']} unresolved; lost "
                                         + (" ".join(l + v["aa"] for l, v in o["sites"].items() if v["state"] == "lost") or "none")})
                else:
                    row.update({"sites_method": "sequence (no model in data/p10_models)", "sites_kept": kept_seq,
                                "sites": sites_seq})
                row["intact_odc"] = row["sites_kept"] >= INTACT_MIN
            standalone.append(row)
        prot_rows += standalone
        missing = set(accs) - found
        st = [fstat[a] for a in accs]
        fus = "intact" if "intact" in st else "degraded" if "degraded" in st else "undetermined"
        ps = proteomes(taxid)
        bp = best_proteome(ps)
        complete_genome = bool(bp and bp["genome"] and (bp["busco"] or 0) >= BUSCO_MIN)
        n_odc = sum(r["odc_family"] for r in standalone)
        n_int = sum(r["intact_odc"] for r in standalone)
        if n_int:
            sa = "yes"
        elif complete_genome:
            sa = "no"
        else:
            sa = "cannot tell"
        sp_rows.append({
            "species": sp, "taxid": taxid, "fusions": ", ".join(f"{a} ({fstat[a]})" for a in accs),
            "fused_odc_half": fus, "census_found_own_fusion": "yes" if not missing else f"NO - missing {', '.join(sorted(missing))}",
            "standalone_plp_proteins": len(standalone), "standalone_odc_family": n_odc, "intact_standalone_odc": n_int,
            "standalone_odc_detail": "; ".join(f"{r['uniprot'] or r['upi']} id {r['identity_to_P11926']:.2f}, "
                                               f"{r['sites_kept']}/13 ({r['sites_method'].split(' (')[0]})"
                                               for r in standalone if r["odc_family"]) or "none",
            "has_intact_standalone_odc": sa,
            "proteome": (f"{bp['id']} {bp['type']}, {bp['proteins']} proteins, assembly {bp['assembly'] or 'none'}, "
                         f"BUSCO {bp['busco_txt']}") if bp else "no UniProt proteome",
            "proteome_source": ("genome" if bp["genome"] else "not genome-assembly based") if bp else "none",
            "complete_genome_proteome": complete_genome,
            "all_proteomes": "; ".join(f"{p['id']} {p['type']} {p['proteins']} {p['assembly'] or '-'} BUSCO {p['busco_txt']}" for p in ps)})
        print(f"{sp:45s} fusion {fus:12s} standalone PLP {len(standalone):2d}, ODC-family {n_odc}, intact {n_int} -> "
              f"{sa:11s} proteome {sp_rows[-1]['proteome'][:70]}", flush=True)
    write("p19_proteins.tsv", prot_rows)
    write("p19_redundancy.tsv", sp_rows)
    report(sp_rows, prot_rows)


def write(name, rows):
    with open(RESULTS / name, "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def table(rows):
    """2x2 counts: rows fused intact / degraded, columns standalone intact yes / no."""
    t = {(f, s): 0 for f in ("intact", "degraded") for s in ("yes", "no")}
    for r in rows:
        if r["fused_odc_half"] in ("intact", "degraded") and r["has_intact_standalone_odc"] in ("yes", "no"):
            t[(r["fused_odc_half"], r["has_intact_standalone_odc"])] += 1
    a, b, c, d = t[("intact", "yes")], t[("intact", "no")], t[("degraded", "yes")], t[("degraded", "no")]
    return (a, b, c, d), (fisher(a, b, c, d) if a + b + c + d else None)


def report(sp, prots):
    (a, b, c, d), pA = table(sp)
    restricted = [r for r in sp if r["complete_genome_proteome"]]
    (a2, b2, c2, d2), pB = table(restricted)
    und = [r for r in sp if r["fused_odc_half"] == "undetermined"]
    cant = [r for r in sp if r["fused_odc_half"] != "undetermined" and r["has_intact_standalone_odc"] == "cannot tell"]
    used = a + b + c + d
    f = lambda p: "not computable" if p is None else f"{p:.3g}"
    L = ["# p19: is the fused ODC half intact only where it is the organism's only ODC?", "",
         "Generated by `scripts/p19_redundancy.py`; definitions in its docstring.", "",
         "## Result", "",
         "Test A - all species with a decided status:", "",
         "| fused ODC half | intact standalone ODC: yes | no |", "|---|---|---|",
         f"| intact | {a} | {b} |", f"| degraded | {c} | {d} |", "",
         f"Fisher exact, two-sided: **p = {f(pA)}** (n = {used} species).", "",
         "Test B - only species with a complete genome-derived proteome "
         f"(genome assembly, BUSCO >= {BUSCO_MIN:.0%}):", "",
         "| fused ODC half | intact standalone ODC: yes | no |", "|---|---|---|",
         f"| intact | {a2} | {b2} |", f"| degraded | {c2} | {d2} |", "",
         f"Fisher exact, two-sided: **p = {f(pB)}** (n = {a2 + b2 + c2 + d2} species).", "",
         f"**Excluded from Test A: {len(und) + len(cant)} of {len(sp)} species** - {len(und)} because the fused ODC "
         f"half is undetermined in p9 ({', '.join(r['species'] for r in und) or 'none'}), {len(cant)} as \"cannot "
         f"tell\": no intact standalone ODC found and no complete genome-derived proteome, so absence cannot be told "
         f"from missing data ({', '.join(r['species'] for r in cant) or 'none'}).", "",
         "## What the counts show, and premises that do not hold", "",
         f"- The prediction was that intact fused ODC halves occur only where there is no intact standalone ODC. "
         f"Of the {a + b} species with an intact fused half in Test A, {a} do have an intact standalone ODC.",
         "- *Gracilariopsis chorda* does not carry an intact standalone ODC by this task's rule (>= 11/13): its only "
         "ODC-family protein, A0A2V3IYX3, is 10 kept / 1 lost (C360R) / 2 unresolved. So the statement that all four "
         "degraded-fusion species carry intact standalone ODCs does not hold; three do (Galdieria partita at 11/13, "
         "Percolomonas at 12/13 by sequence, Rhodosorus at 13/13).",
         f"- The green-lineage fusions are not all intact: p9 leaves {len(und)} fusion species undetermined (listed "
         "above), and they are excluded, not counted as intact.",
         f"- Only {c + d} species have a degraded fused half, so neither test can detect anything but a very strong "
         "association; p = 1 here means no association was detected, not that none exists.",
         "- Proteome information comes from the UniProt proteomes API (cached in data/p19/); data/q_proteomes/ holds "
         "only six FASTA files and no completeness data.",
         "- Short alignments: identities computed over fewer than 100 aligned residues are unreliable (Adiantum "
         "A0A9D4ZKC7 0.381 over 42, Tetraselmis A0A6T8DES7 0.333 over 48, Pycnococcus UPI003D4D0C8E 0.244 over 41, "
         "Auxenochlorella A0A1D2A4Z4 0.224 over 58). None changes a species call: a fragment that short cannot "
         "cover 11 of the 13 sites, so it could not count as an intact ODC either way.",
         "- Sites were read from structure only for the five red-algal proteins with a model already in "
         "data/p10_models/; every green-lineage protein was read by sequence.", "",
         "## Per species", "",
         "| species | fusion(s) (p9 gated) | fused ODC half | standalone PLP | ODC-family | intact (>= 11/13) | intact standalone ODC? | proteome | complete genome proteome |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in sp:
        L.append(f"| *{r['species']}* | {r['fusions']} | {r['fused_odc_half']} | {r['standalone_plp_proteins']} | "
                 f"{r['standalone_odc_family']} | {r['intact_standalone_odc']} | **{r['has_intact_standalone_odc']}** | "
                 f"{r['proteome']} | {r['complete_genome_proteome']} |")
    bad = [r for r in sp if r["census_found_own_fusion"] != "yes"]
    L += ["", "Census check (each species' UniParc census contains its own fusion): "
          + ("passed for every species." if not bad else "FAILED for " + "; ".join(f"{r['species']} ({r['census_found_own_fusion']})" for r in bad) + "."), "",
          f"## Every standalone PLP-decarboxylase (ODC family = identity to P11926 >= {ODC_ID:.1%})", "",
          "Identity is from a local alignment to human ODC; sites by structure where a model was already in "
          "data/p10_models/, otherwise by sequence. Proteins near the threshold are the borderline calls.", "",
          "| species | UniProt / UPI | length | identity to P11926 | aligned | ODC family | sites | method | intact |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(prots, key=lambda r: (r["species"], -r["identity_to_P11926"])):
        L.append(f"| *{r['species']}* | {r['uniprot'] or r['upi']} | {r['length']} | {r['identity_to_P11926']:.3f} | "
                 f"{r['aligned_residues']} | {'yes' if r['odc_family'] else 'no'} | {r['sites'] or '-'} | "
                 f"{r['sites_method'] or '-'} | {'yes' if r['intact_odc'] else ''} |")
    ids = sorted(r["identity_to_P11926"] for r in prots)
    lo = max((x for x in ids if x < ODC_ID), default=None)
    hi = min((x for x in ids if x >= ODC_ID), default=None)
    L += ["", f"Closest identities on either side of the threshold: {lo} below, {hi} above.", ""]
    (RESULTS / "p19_report.md").write_text("\n".join(L), encoding="utf8")


if __name__ == "__main__":
    main()
