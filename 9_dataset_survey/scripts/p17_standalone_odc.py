"""Track B p17: are the standalone PLP-decarboxylases of the degraded-fusion species intact?

SEQUENCE alignment only (not the p9 structural TM-align + CA-offset + pLDDT gate). Each protein is aligned to
human ODC P11926 (local, BLOSUM62, gap open 11 / extend 1 - the aligner used in p10) and the aligned residue
at each of the 13 human-ODC functional positions is reported. A position falling in a gap or outside the
local alignment is "-" (not aligned). Calls are not merged into any p9 table.

Sequences: UniProt REST JSON, cached in data/p17_standalone/ (no email or identifier sent).
A0AAV8UQI0 is reused from data/p16_colabfold/06_*_COLAB.txt (first chain of the ColabFold A:B input).
P11926 from data/references/decarboxylases.fasta (p6 cache).

Writes results/p17_standalone_odc.tsv (one row per protein per site) and results/p17_report.md
Run:  py -3.12 p17_standalone_odc.py
"""
import csv
import json
import sys

from Bio.Align import PairwiseAligner, substitution_matrices

from p2_fetch_annotations import DATA, RESULTS, get
from p6_decarboxylase_partner import references
from p9_self_processing import ODC_SITE_LABELS

sys.stdout.reconfigure(encoding="utf8")
OUT = DATA / "p17_standalone"
OUT.mkdir(parents=True, exist_ok=True)
PROTEINS = [("Galdieria partita", "A0A9C7PSP6"), ("Galdieria partita", "A0A9C7PSB6"),
            ("Galdieria partita", "A0A9C7Q3X8"), ("Galdieria partita", "A0A9C7PZK6"),
            ("Rhodosorus marinus", "A0AAV8UYK6"), ("Rhodosorus marinus", "A0AAV8UQI0"),
            ("Gracilariopsis chorda", "A0A2V3IYX3"), ("Gracilariopsis chorda", "A0A2V3IKW4"),
            ("Gracilariopsis chorda", "A0A2V3IW11")]
REUSED = {"A0AAV8UQI0": "06_*_COLAB.txt"}
LOCAL = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                        open_gap_score=-11, extend_gap_score=-1)


def fetch(acc):
    if acc in REUSED:
        path = next((DATA / "p16_colabfold").glob(REUSED[acc]))
        return path.read_text(encoding="utf8").split(":")[0].strip(), f"reused {path.relative_to(DATA.parent)}", "", ""
    p = OUT / f"{acc}.json"
    if not p.exists():
        st, _, body = get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
        if st != 200:
            raise RuntimeError(f"{acc}: HTTP {st}")
        p.write_text(body, encoding="utf8")
    d = json.loads(p.read_text(encoding="utf8"))
    org = d.get("organism", {}).get("scientificName", "")
    name = (d.get("proteinDescription", {}).get("recommendedName") or
            (d.get("proteinDescription", {}).get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    return d["sequence"]["value"], f"UniProt REST (cached {p.relative_to(DATA.parent)})", org, name


def align(q, ref):
    a = LOCAL.align(q, ref)[0]
    pairs = [(i, j) for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))]
    ident = sum(q[i] == ref[j] for i, j in pairs)
    return {j + 1: (i + 1, q[i]) for i, j in pairs}, ident, len(pairs), a.score


def main():
    odc = references()["P11926"]
    rows, summ = [], []
    for species, acc in PROTEINS:
        seq, src, org, name = fetch(acc)
        if org and org != species:
            print(f"WARNING {acc}: UniProt organism {org!r} differs from {species!r}")
        m, ident, n, score = align(seq, odc)
        q_aligned = sorted(i for i, _ in m.values())
        r_aligned = sorted(m)
        kept = lost = unal = 0
        calls = []
        for lab in ODC_SITE_LABELS:
            pos, hum = int(lab[1:]), lab[0]
            qpos, aa = m.get(pos, ("", "-"))
            state = "not aligned" if aa == "-" else "kept" if aa == hum else "substituted"
            kept += state == "kept"
            lost += state == "substituted"
            unal += state == "not aligned"
            calls.append(f"{lab}{'=' if state == 'kept' else ':' + aa}")
            rows.append({"species": species, "accession": acc, "length": len(seq), "site": lab, "human_aa": hum,
                         "aligned_aa": aa, "protein_position": qpos, "state": state,
                         "identity_to_P11926_pct": round(100 * ident / n, 1), "aligned_pairs": n,
                         "method": "sequence alignment (local BLOSUM62 11/1) to P11926; NOT the p9 structural gate"})
        summ.append({"species": species, "acc": acc, "len": len(seq), "name": name, "src": src,
                     "ident": round(100 * ident / n, 1), "n": n,
                     "qspan": f"{q_aligned[0]}-{q_aligned[-1]}", "rspan": f"{r_aligned[0]}-{r_aligned[-1]}",
                     "kept": kept, "lost": lost, "unal": unal, "calls": " ".join(calls),
                     "c360": next(r for r in rows if r["accession"] == acc and r["site"] == "C360")})
        print(f"{acc} {species}: {ident}/{n} identical ({100 * ident / n:.1f} %), {kept} kept / {lost} substituted / {unal} not aligned")
    with open(RESULTS / "p17_standalone_odc.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    report(summ)


def report(summ):
    L = ["# p17: standalone PLP-decarboxylases of the degraded-fusion species", "",
         "Generated by `scripts/p17_standalone_odc.py`.", "",
         "**Method: sequence alignment, not structure.** Each protein is aligned to human ODC P11926 (local "
         "alignment, BLOSUM62, gap open 11 / extend 1, as in p10) and the residue aligned to each of the 13 human-ODC "
         "functional positions is read off. This is **not** the p9 test (TM-align, CA within 4 A, pLDDT >= 70); "
         "the calls are in their own file and are not merged into any p9 table. A call is only as good as the "
         "alignment: the lower the identity to P11926, the less a single-residue call can be trusted. "
         "\"Substituted\" here means a different residue is aligned at that position, nothing more.", "",
         "| species | accession | length | identity to P11926 | aligned pairs | protein span / P11926 span | kept / substituted / not aligned | sites |",
         "|---|---|---|---|---|---|---|---|"]
    for s in summ:
        L.append(f"| *{s['species']}* | {s['acc']} | {s['len']} | {s['ident']} % | {s['n']} | {s['qspan']} / {s['rspan']} | "
                 f"{s['kept']} / {s['lost']} / {s['unal']} | {s['calls']} |")
    L += ["", "Sequence source: UniProt REST, cached in `data/p17_standalone/`; A0AAV8UQI0 reused from "
          "`data/p16_colabfold/06_*_COLAB.txt`.", "",
          "## Gracilariopsis chorda: one standalone ODC-family protein or three?", ""]
    g = [s for s in summ if s["species"] == "Gracilariopsis chorda"]
    for s in g:
        c = s["c360"]
        L.append(f"- **{s['acc']}** ({s['len']} aa; UniProt name \"{s['name'] or '-'}\"): {s['ident']} % identical to "
                 f"human ODC over {s['n']} aligned residues; position C360 = **{c['aligned_aa']}**"
                 f"{' (protein residue ' + str(c['protein_position']) + ')' if c['protein_position'] else ''}.")
    L += ["", "**Which is right.** Both statements count different things, and both are correct for what they count:", "",
          "- `p7_same_species.tsv` lists **three** standalone proteins because it counts every protein with a "
          "PLP-decarboxylase domain (PF02784) and no AdoMetDC domain (\"decarboxylase only\"). That column does not "
          "say ODC.",
          "- Only **one** of the three, A0A2V3IYX3, is ODC-family: it is 45.7 % identical to human ODC, in the same "
          "range as the Galdieria and Rhodosorus ODC-family proteins in the table above (45-48 %), and UniProt names "
          "it ornithine decarboxylase. A0A2V3IKW4 and A0A2V3IW11 are 27-28 % identical to human ODC, in the range of "
          "the proteins above that UniProt or p10 call diaminopimelate decarboxylase (25-28 %), and UniProt names "
          "both diaminopimelate decarboxylase.",
          "- So the earlier note is right if read as \"one standalone **ODC-family** protein\", and that protein has "
          "**R at C360** by this sequence alignment (as in p10). The other two have **C at C360**, but they are not "
          "ODC-family, and these 13 positions alone do not show what they decarboxylate.", "",
          "## Note on structures", "",
          "All nine proteins do have AlphaFold DB models (downloaded in p10, `data/p10_models/AF-<accession>.pdb`), "
          "and p10 already made structural reads of them (`results/p10_species_enzymes.tsv`). They are not used "
          "here, because this task asked for a sequence alignment.", ""]
    (RESULTS / "p17_report.md").write_text("\n".join(L), encoding="utf8")


if __name__ == "__main__":
    main()
