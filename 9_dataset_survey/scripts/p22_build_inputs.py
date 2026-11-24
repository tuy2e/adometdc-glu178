"""p22: build input files only (no analysis is run).

1. tblastn junction queries for A0A9C7PQ51 (Galdieria partita) and A0A1Y1HLG4 (Klebsormidium nitens), p14 recipe:
   last 60 residues of the AdoMetDC Pfam match + the whole linker + first 60 residues of the PLP match.
   Boundaries from results/p7_linkers.tsv (adometdc_end, plp_start). Sequence read from CA atoms of Marc's models
   (read only) and checked against the cached UniProt sequence.
2. ColabFold heterodimer inputs, SEQ1:SEQ2 on one line (as data/p16_colabfold):
   (a) A0A9C7PQ51 + A0A9C7PSP6, (b) A0A2V3IZY9 + A0A2V3IYX3. Fusions from Marc's models, standalones from the
   p17 UniProt cache.
3. Species census assembled from existing results: p19_redundancy.tsv / p19_proteins.tsv (census, standalone scores,
   proteome), p21b_rescore.tsv (fused ODC half, core rule), data/p21/taxonomy_ranks.json (lineage).

Writes data/p22_tblastn_queries/*.txt + README.md, data/p22_colabfold/*_COLAB.txt + README.md,
results/p22_census.tsv, results/p22_census.md
"""
import csv
import json
import sys
from collections import defaultdict

from p2_fetch_annotations import DATA, RESULTS
from p3_residue_grid import RUN

sys.stdout.reconfigure(encoding="utf8")
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}
TB = DATA / "p22_tblastn_queries"
CF = DATA / "p22_colabfold"
PHYLA = {"Rhodophyta", "Chlorophyta", "Streptophyta", "Heterolobosea", "Apicomplexa"}


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def model_seq(acc):
    res = {}
    with open(RUN / "models" / f"{acc}.pdb", encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                res[int(ln[22:26])] = THREE[ln[17:20]]
    return "".join(res[k] for k in sorted(res))


def uniprot_seq(acc, folder="uniprot"):
    d = json.loads((DATA / folder / f"{acc}.json").read_text(encoding="utf8"))
    e = d.get("entry", d)
    return (e.get("sequence") or {}).get("value")


def fusion_seq(acc):
    m = model_seq(acc)
    u = uniprot_seq(acc)
    return m, ("model sequence identical to cached UniProt" if u == m else
               "no cached UniProt sequence (entry inactive); model sequence used" if not u else
               f"MODEL DIFFERS from cached UniProt ({len(m)} vs {len(u)} aa); model used")


def main():
    TB.mkdir(parents=True, exist_ok=True)
    CF.mkdir(parents=True, exist_ok=True)
    out = []
    link = {r["accession"]: r for r in read("p7_linkers.tsv")}
    q = []
    for n, (acc, org, taxid) in enumerate([("A0A9C7PQ51", "Galdieria partita", 83374),
                                           ("A0A1Y1HLG4", "Klebsormidium nitens", 105231)], 1):
        s, note = fusion_seq(acc)
        a_end, p_start = int(link[acc]["adometdc_end"]), int(link[acc]["plp_start"])
        lin = p_start - a_end - 1
        assert s[a_end:p_start - 1] == link[acc]["linker_seq"], acc
        query = s[a_end - 60:p_start - 1 + 60]
        name = f"{n:02d}_{org.replace(' ', '_')}_{acc}.txt"
        (TB / name).write_text("\n".join(query[i:i + 60] for i in range(0, len(query), 60)) + "\n", encoding="utf8")
        out.append(TB / name)
        q.append((name, acc, org, taxid, len(query), lin, a_end, p_start, note))
        print(acc, len(query), note)
    genome = defaultdict(set)
    with open(RESULTS.parent / "results_81" / "q3_assemblies.tsv", encoding="utf8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["study_accession"]:
                genome[r["species"]].add(r["study_accession"])
    R = ["# p22 tblastn junction queries: the intron-dependent fusions", "",
         "Built by `scripts/p22_build_inputs.py` from Marc's models (read only; sequence from CA atoms, checked against "
         "UniProt) and the domain boundaries in `results/p7_linkers.tsv`. Same recipe as `data/p14_tblastn_queries/`: "
         "the **last 60 residues of the AdoMetDC Pfam match**, the **whole linker**, and the **first 60 residues of the "
         "PLP-decarboxylase Pfam match**. Positions below are in the full protein.", "",
         "## Why these two", "",
         "p21f lists 10 fusions whose gene model puts a predicted intron in the linker. For them, the gene structure "
         "says nothing for or against a real fusion, because a gene-prediction merge of two neighbouring genes looks "
         "exactly like that. A transcript read across the junction is the test.", "",
         "| file | protein | type into NCBI's organism box | query length | AdoMetDC / linker / PLP residues | sequence source | genome BioProject(s), from q3 |",
         "|---|---|---|---|---|---|---|"]
    for name, acc, org, taxid, L, lin, a_end, p_start, note in q:
        R.append(f"| `{name}` | {acc} | *{org}* (taxid {taxid}) | {L} aa | {a_end - 59}-{a_end} / {a_end + 1}-{p_start - 1} "
                 f"({lin} aa) / {p_start}-{p_start + 59} | {note} | {', '.join(sorted(genome.get(org, []))) or 'none in q3'} |")
    R += ["", "## How to run", "",
          "blast.ncbi.nlm.nih.gov -> tblastn -> paste one file -> Database: **Transcriptome Shotgun Assembly (TSA)** -> "
          "Organism box: the name above (let NCBI autocomplete it) -> BLAST. If nothing comes back, try the genus, then "
          "the database **Nucleotide collection (nr/nt)**.", "",
          "From p15: through the URL API, restricting `tsa_nt` with `ENTREZ_QUERY` left searches waiting for over 50 "
          "minutes; the web form restricts TSA through its organism box instead, which is why the web form is "
          "suggested here.", "",
          "## What a positive looks like", "",
          "- **One subject sequence (a TSA contig or an mRNA) aligning across the whole query in a single HSP**: starting "
          "inside the AdoMetDC part, running through the entire linker, and ending inside the PLP part, with no break "
          "at the linker. That is one RNA molecule carrying both halves, and here it also shows the predicted linker "
          "intron is really spliced out as predicted.",
          "- Hits covering only the AdoMetDC part or only the PLP part are **not** evidence: they can come from separate "
          "genes.",
          "- No hit means nothing on its own: the species may have little or no RNA data in TSA (not checked here).",
          "- If the hit's TSA BioProject is the genome's BioProject (table above), the same project built both, so it is "
          "**not independent** evidence; record that with the hit.", "",
          "## Note on the p14 choice", "",
          "The task calls the p14 choice of single-exon genes backwards. Partly: a single-exon fusion cannot be an "
          "intron-spanning merge, so these two intron-dependent genes do need transcripts more. But the single-exon "
          "cases are not settled either: `trackB_fusion_validation.md` lists read-through of adjacent genes and "
          "assembly errors (a lost stop codon, an indel) as ways an uninterrupted ORF can still be artefactual. "
          "Transcript evidence helps both groups.", ""]
    (TB / "README.md").write_text("\n".join(R), encoding="utf8")
    out.append(TB / "README.md")

    # 2. ColabFold heterodimers
    jobs = [("a", "Galdieria_partita", "A0A9C7PQ51", "A0A9C7PSP6"), ("b", "Gracilariopsis_chorda", "A0A2V3IZY9", "A0A2V3IYX3")]
    sizes = []
    for k, sp, fus, sa in jobs:
        f, fnote = fusion_seq(fus)
        s = uniprot_seq(sa, "p17_standalone")
        name = f"{k}_{sp}_{fus}_fusion+{sa}_standalone_COLAB.txt"
        (CF / name).write_text(f"{f}:{s}\n", encoding="utf8")
        out.append(CF / name)
        sizes.append((k, sp, fus, len(f), sa, len(s), name, fnote))
    C = ["# p22 ColabFold heterodimers: fused ODC half + same-species standalone ODC", "",
         "Each `*_COLAB.txt` is `FUSION:STANDALONE` on one line (chain A = fusion, chain B = standalone). Settings as in "
         "`data/p16_colabfold/README.md` (alphafold2_multimer_v3, unpaired_paired, 4 seeds, save_all). Fusion sequences "
         "are from Marc's models (read only), standalones from the UniProt records cached in `data/p17_standalone/`.", "",
         "| job | chain A (fusion) | aa | chain B (standalone) | aa | total residues | file |", "|---|---|---|---|---|---|---|"]
    for k, sp, fus, lf, sa, ls_, name, _ in sizes:
        C.append(f"| {k} | {fus} | {lf} | {sa} | {ls_} | **{lf + ls_}** | `{name}` |")
    C += ["", "## (a) Galdieria partita: A0A9C7PQ51 + A0A9C7PSP6", "",
          "- A0A9C7PQ51: ODC half degraded (p9 gated: 1 kept / 7 lost / 5 unresolved; lost K69R H197S E274N P275V G276D "
          "R277S C360G).",
          "- A0A9C7PSP6: ODC-family standalone, called intact in p19 by the >= 11/13 rule. It is **not 13/13**: by "
          "structure 11 kept / 1 lost (Y331S) / 1 unresolved (p10); Y331 is outside p21b's core list.",
          "- The question: does the fused ODC half pair with the standalone ODC at all (inter-chain PAE between the "
          "fusion's ODC half and chain B), and in the ODC dimer orientation?", "",
          "## (b) Gracilariopsis chorda: A0A2V3IZY9 + A0A2V3IYX3", "",
          "**The test: whether K69 of the standalone and C360 of the fusion come together the way they do in human ODC "
          "1D7K.** The two halves are broken at opposite positions, and in a real ODC dimer those two residues come "
          "from different chains (each active site takes Lys69 from one subunit and Cys360 from the other; in 1D7K, "
          "K69 of one chain to C360 of the other is 8.7 / 8.5 A CA-CA, p13). So one of the two active sites of this "
          "heterodimer could be complete and the other doubly broken.", "",
          "| residue (human ODC numbering) | fusion A0A2V3IZY9 | standalone A0A2V3IYX3 |", "|---|---|---|",
          "| K69 | **Arg429** (lost; p9 structure, pLDDT 86) | **Lys97** (kept; p17 sequence alignment) |",
          "| C360 | **Cys680** (kept, but 3.5 A off in p9, close to the 4 A limit) | **Arg380** (lost; p17 sequence, p10 structure) |", "",
          "What to measure: CA-CA distance Lys97(B) to Cys680(A), and Arg429(A) to Arg380(B), against 8.5-8.7 A in 1D7K, "
          "plus the inter-chain PAE between the two ODC parts. The standalone positions come from a sequence alignment "
          "(p17), not a structural one; check them on the predicted model before measuring.", "",
          "Caveat: the fusion's ODC half has other losses too (H197E, E274V, P275V, R277A, p9), so a complete K69 + C360 "
          "pair would not by itself make a working active site.", ""]
    (CF / "README.md").write_text("\n".join(C), encoding="utf8")
    out.append(CF / "README.md")

    # 3. census
    sp = read("p19_redundancy.tsv")
    prots = defaultdict(list)
    for r in read("p19_proteins.tsv"):
        if r["odc_family"] == "True":
            prots[r["species"]].append(r)
    core = {r["accession"]: r["call_core"] for r in read("p21b_rescore.tsv")}
    tax = json.loads((DATA / "p21" / "taxonomy_ranks.json").read_text(encoding="utf8"))
    p10lin = {(x["in_our_23"] or ""): x["lineage"] for x in read("p10_fusion_search.tsv") if x["in_our_23"]}
    rows = []
    for r in sp:
        accs = [x.split(" (")[0] for x in r["fusions"].split(", ")]
        calls = [core[a] for a in accs]
        fused = "intact" if "intact" in calls else "degraded" if "degraded" in calls else "undetermined"
        t = tax.get(r["taxid"], {})
        ranked = [n for n, rk in t.get("lineage", []) if rk in ("phylum", "class", "order", "family")]
        if len(ranked) >= 3:     # UniProt ranks: order them phylum -> family
            lin = "; ".join(ranked if ranked[0] in PHYLA else ranked[::-1])
        else:                    # taxid not in the p21 cache, or few ranked levels: p10's lineage (already root first)
            lin = "; ".join(p10lin.get(accs[0], "").split("; ")[-5:-1])
        best = max(prots[r["species"]], key=lambda p: (int(p["sites_kept"] or 0), p["sites_method"].startswith("structure")),
                   default=None)
        rows.append({"species": r["species"], "lineage": lin,
                     "fusions": ", ".join(f"{a} ({c})" for a, c in zip(accs, calls)),
                     "fused_odc_half_core_rule": fused,
                     "standalone_odc_family": r["standalone_odc_family"],
                     "best_standalone": (f"{best['uniprot'] or best['upi']} {best['sites_kept']}/13 "
                                         f"({best['sites_method'].split(' (')[0]})") if best else "none",
                     "intact_standalone_odc": r["has_intact_standalone_odc"],
                     "proteome": r["proteome"], "complete_genome_proteome": r["complete_genome_proteome"]})
    with open(RESULTS / "p22_census.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    out.append(RESULTS / "p22_census.tsv")
    M = ["# p22: species census, fused ODC half versus standalone ODCs", "",
         "Assembled by `scripts/p22_build_inputs.py` from existing results; nothing was recomputed.", "",
         "- **Fused ODC half**: p21b core rule (12 sites, Y331 dropped) on p9's gated structural reads; species with "
         "several fusions are intact if any fusion is intact, degraded if any is degraded and none intact.",
         "- **Standalones**: the p19 UniParc census; ODC family = identity to human ODC >= 36.5 %; site score = kept of "
         "the 13 sites (structure where a model was in data/p10_models, else sequence); **intact standalone = >= 11/13** "
         "(the p19 rule, scored on 13 sites, not the p21b core). \"cannot tell\" = none found and no complete "
         "genome-derived proteome.",
         "- **Proteome**: the UniProt proteome with the best genome/BUSCO evidence (p19); complete = genome assembly "
         "and BUSCO >= 80 %.", "",
         "| species | lineage | fusion(s) (core call) | fused ODC half | ODC-family standalones | best standalone | intact standalone ODC | proteome | complete genome proteome |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["fused_odc_half_core_rule"], r["species"])):
        M.append(f"| *{r['species']}* | {r['lineage']} | {r['fusions']} | {r['fused_odc_half_core_rule']} | "
                 f"{r['standalone_odc_family']} | {r['best_standalone']} | {r['intact_standalone_odc']} | {r['proteome']} | "
                 f"{r['complete_genome_proteome']} |")
    M += ["", "Two scales are in use, and they are labelled: the fused half is called on the 12-site core with p9's "
          "strict all-kept rule; standalones are scored on 13 sites with p19's >= 11 rule. They are not the same test.", ""]
    (RESULTS / "p22_census.md").write_text("\n".join(M), encoding="utf8")
    out.append(RESULTS / "p22_census.md")
    print("\n".join(str(p) for p in out))


if __name__ == "__main__":
    main()
