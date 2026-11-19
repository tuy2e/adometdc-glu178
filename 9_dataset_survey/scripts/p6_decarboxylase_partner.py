"""Track B, question 5: which pyridoxal-dependent decarboxylase is fused to AdoMetDC?

PF02784/PF00278 cover ornithine, arginine, diaminopimelate (lysine) and carboxynorspermidine
decarboxylases, and the inactive antizyme inhibitors. For every protein whose architecture contains
PF01536 followed by PF02784/PF00278, the region after the last AdoMetDC Pfam match is compared with
reference enzymes of known activity (pairwise local alignment, BLOSUM62, gap open 11 / extend 1), and
the two residues that make ODC catalytically active are read from the alignment to human ODC:
Lys69 (PLP Schiff base) and Cys360 (catalytic acid). Antizyme inhibitors lack these.

This is nearest-reference similarity, not a substrate assay: an ODC-fold enzyme can decarboxylate
arginine (e.g. chlorella virus PBCV-1 A207R), so "closest to ODC" does not by itself establish the
substrate. The table says which reference is closest, by how much, and whether the catalytic residues
are present.

Writes results/p6_decarboxylase_partner.tsv
Run (seconds; reference sequences cached in data/references/decarboxylases.fasta):
    py -3.12 p6_decarboxylase_partner.py
"""
import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

from p2_fetch_annotations import DATA, RESULTS, get, read_candidates

REFS = {  # accession: (label, activity basis)
    "P11926": ("ODC human", "reviewed, EC 4.1.1.17"),
    "P07805": ("ODC Trypanosoma brucei", "reviewed, EC 4.1.1.17"),
    "P08432": ("ODC yeast", "reviewed, EC 4.1.1.17"),
    "P93351": ("ODC tobacco (plant)", "reviewed, EC 4.1.1.17"),
    "Q8IJ77": ("AdoMetDC/ODC Plasmodium falciparum", "bifunctional; ODC activity shown (Muller 2000)"),
    "Q84527": ("ADC chlorella virus PBCV-1 (ODC-fold)", "arginine decarboxylase (Shah 2004)"),
    "Q9SI64": ("ADC1 Arabidopsis", "reviewed, EC 4.1.1.19"),
    "O23141": ("ADC2 Arabidopsis", "reviewed, EC 4.1.1.19"),
    "P21170": ("ADC SpeA E. coli", "reviewed, EC 4.1.1.19"),
    "P00861": ("DAPDC LysA E. coli", "reviewed, EC 4.1.1.20"),
    "Q56575": ("CANSDC Vibrio alginolyticus", "reviewed, EC 4.1.1.96"),
    "A8FNH9": ("CANSDC Campylobacter jejuni", "reviewed, EC 4.1.1.96"),
    "O14977": ("AZIN1 human (inactive ODC paralog)", "reviewed, no ODC activity"),
    "Q96A70": ("AZIN2 human (ODC-like)", "reviewed"),
}
ADOMETDC = {"PF01536", "PF02675"}
PLP_DC = {"PF02784", "PF00278"}


def references():
    p = DATA / "references" / "decarboxylases.fasta"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        st, _, body = get("https://rest.uniprot.org/uniprotkb/accessions?accessions="
                          + ",".join(REFS) + "&format=fasta")
        p.write_text(body, encoding="utf8")
    return {r.id.split("|")[1]: str(r.seq) for r in SeqIO.parse(str(p), "fasta")}


def main():
    refs = references()
    seqs = read_candidates()
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    with open(RESULTS / "architecture_per_protein.tsv", encoding="utf8") as fh:
        arch = [r for r in csv.DictReader(fh, delimiter="\t")
                if set(r["arch_pfam_ids"].split(" > ")) & ADOMETDC and set(r["arch_pfam_ids"].split(" > ")) & PLP_DC]
    rows = []
    for r in arch:
        doms = [(d.split()[0], *map(int, d.split()[1].split("-"))) for d in r["domains_with_positions"].split("; ")]
        start = max(e for pf, s, e in doms if pf in ADOMETDC) + 1
        q = seqs[r["accession"]][start - 1:]
        hits = []
        for acc, (label, basis) in REFS.items():
            a = al.align(q, refs[acc])[0]
            ident = cols = 0
            for (qs, qe), (ts, te) in zip(*a.aligned):
                for i, j in zip(range(qs, qe), range(ts, te)):
                    cols += 1
                    ident += q[i] == refs[acc][j]
            hits.append({"acc": acc, "label": label, "score": a.score, "ident": ident / cols if cols else 0,
                         "cov": cols / len(q), "aln": a})
        hits.sort(key=lambda h: -h["score"])
        # catalytic residues via human ODC
        h = next(x for x in hits if x["acc"] == "P11926")
        qmap = {}
        for (qs, qe), (ts, te) in zip(*h["aln"].aligned):
            for i, j in zip(range(qs, qe), range(ts, te)):
                qmap[j + 1] = (q[i], i + start)
        k69, c360 = qmap.get(69, ("-", "")), qmap.get(360, ("-", ""))
        best, second = hits[0], hits[1]
        best_other_family = next(x for x in hits[1:] if x["label"].split()[0] != best["label"].split()[0])
        rows.append({"accession": r["accession"], "organism": r["organism"], "fold_ok": r["fold_ok"],
                     "source": r["source"], "architecture": r["arch_pfam_ids"], "partner_region": f"{start}-{len(seqs[r['accession']])}",
                     "closest_reference": best["label"], "closest_acc": best["acc"],
                     "closest_identity": round(best["ident"], 3), "closest_coverage": round(best["cov"], 3),
                     "closest_score": round(best["score"]),
                     "second": f"{second['label']} ({second['ident']:.2f} id, score {second['score']:.0f})",
                     "best_other_family": f"{best_other_family['label']} ({best_other_family['ident']:.2f} id, score {best_other_family['score']:.0f})",
                     "score_margin_vs_other_family": round(best["score"] - best_other_family["score"]),
                     "ODC_Lys69": f"{k69[0]}{k69[1]}", "ODC_Cys360": f"{c360[0]}{c360[1]}",
                     "all_scores": "; ".join(f"{x['label']}={x['score']:.0f}/{x['ident']:.2f}" for x in hits)})
    with open(RESULTS / "p6_decarboxylase_partner.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for x in rows:
        print(f"{x['accession']:<11} {x['organism'][:30]:<30} {x['architecture']:<38} -> {x['closest_reference']:<38} "
              f"id {x['closest_identity']:.2f} cov {x['closest_coverage']:.2f} | margin vs {x['best_other_family'][:40]:<40} "
              f"{x['score_margin_vs_other_family']:>4} | K69 {x['ODC_Lys69']:<6} C360 {x['ODC_Cys360']}")


if __name__ == "__main__":
    main()
