"""Phase 1, step 4: build master_candidates.tsv for all 1,000 FOLDTRACE candidates.

Reads (never writes) Marc's FOLDTRACE run folder:
    foldseek_hits.m8      Foldseek web-server format (21 columns) + appended taxonomy report
    foldtrace_calls.tsv   per-site calls; its foldseek_* columns are mislabelled, so all
                          Foldseek values here are taken from the raw m8 instead
    models/*.pdb          AlphaFold models (pLDDT in the B-factor column)

Writes:
    9_dataset_survey/results/master_candidates.tsv
    9_dataset_survey/results/master_candidates_summary.txt

Run from anywhere:
    py -3.12 p1_build_master_table.py
"""
import csv
import re
import statistics as st
from collections import Counter
from pathlib import Path

RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
OUT = Path(__file__).resolve().parent.parent / "results"
SITES = ["Ser68", "Ser229", "His243", "Cys82", "Glu247"]
AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
           LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
# Foldseek web-server m8 columns
M8 = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
      "prob", "evalue", "bits", "qlen", "tlen", "qaln", "taln", "tca", "tseq", "taxid", "taxname"]


def accession(s: str) -> str:
    return re.search(r"AF-([A-Z0-9]+)-F\d+", s).group(1)


def read_m8() -> dict:
    hits = {}
    with open(RUN / "foldseek_hits.m8", encoding="utf8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) != 21:          # skip the appended taxonomy report
                continue
            r = dict(zip(M8, row))
            r["rank"] = len(hits) + 1   # order in the Foldseek output
            r["description"] = row[1].split(" ", 1)[1] if " " in row[1] else ""
            hits[accession(row[1])] = r
    return hits


def read_model(acc: str) -> dict:
    """Per-residue pLDDT keyed by residue number, from CA atoms."""
    plddt = {}
    with open(RUN / "models" / f"{acc}.pdb", encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                plddt[int(ln[22:26])] = float(ln[60:66])
    return plddt


def main() -> None:
    hits = read_m8()
    with open(RUN / "foldtrace_calls.tsv", encoding="utf8") as fh:
        calls = {accession(c["candidate_id"]): c for c in csv.DictReader(fh, delimiter="\t")}
    models = {p.stem for p in (RUN / "models").glob("*.pdb")}
    assert set(hits) == set(calls) == models, "accession sets differ between m8, calls and models"

    rows = []
    for acc, h in sorted(hits.items(), key=lambda kv: kv[1]["rank"]):
        c = calls[acc]
        pl = read_model(acc)
        vals = list(pl.values())
        qlen, tlen = int(h["qlen"]), int(h["tlen"])
        row = {
            "accession": acc,
            "foldseek_rank": h["rank"],
            "foldtrace_rank": c["rank"],
            "description": h["description"],
            "fragment_in_name": "(Fragment)" in h["description"],
            "organism": h["taxname"],
            "taxid": h["taxid"],
            "length": tlen,
            "plddt_mean": round(st.mean(vals), 1),
            "plddt_frac_below70": round(sum(v < 70 for v in vals) / len(vals), 3),
            "fs_prob": h["prob"],
            "fs_evalue": h["evalue"],
            "fs_bits": h["bits"],
            "fs_fident_pct": h["fident"],
            "fs_alnlen": h["alnlen"],
            "fs_qstart": h["qstart"], "fs_qend": h["qend"],
            "fs_tstart": h["tstart"], "fs_tend": h["tend"],
            "fs_qcov": round((int(h["qend"]) - int(h["qstart"]) + 1) / qlen, 3),
            "fs_tcov": round((int(h["tend"]) - int(h["tstart"]) + 1) / tlen, 3),
            "residues_before_hit": int(h["tstart"]) - 1,
            "residues_after_hit": tlen - int(h["tend"]),
            "tm_norm_ref": c["tmalign_tm_norm_ref"],
            "rmsd": c["tmalign_rmsd"],
            "fold_ok": c["fold_ok"],
            "verdict": c["verdict"],
        }
        for s in SITES:
            obs, off, state = c[f"{s}_obs"], c[f"{s}_offsetA"], c[f"{s}_state"]
            if c["fold_ok"] != "true":
                reason = "below_tm_gate"
            elif state != "unresolved":
                reason = ""
            elif obs in ("", "-"):
                reason = "no_aligned_residue"
            else:
                reason = "aligned_residue_too_far"
            num = re.search(r"\d+", obs or "")
            row[f"{s}_obs"] = obs
            row[f"{s}_offsetA"] = off
            row[f"{s}_state"] = state
            row[f"{s}_unresolved_reason"] = reason
            row[f"{s}_plddt"] = pl.get(int(num.group())) if num else ""
        rows.append(row)

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "master_candidates.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    ok = [r for r in rows if r["fold_ok"] == "true"]
    lines = [f"candidates: {len(rows)}   fold_ok: {len(ok)}   fold-failing: {len(rows) - len(ok)}",
             f"verdict (fold_ok): {dict(Counter(r['verdict'] for r in ok))}",
             f"fragment in name: {sum(r['fragment_in_name'] for r in rows)}",
             f"mean pLDDT < 70: {sum(r['plddt_mean'] < 70 for r in rows)}",
             f">= 100 residues outside the Foldseek hit (possible fusions): "
             f"{sum(r['residues_before_hit'] + r['residues_after_hit'] >= 100 for r in rows)}"]
    for s in SITES:
        lines.append(f"{s} (fold_ok): states {dict(Counter(r[f'{s}_state'] for r in ok))}; "
                     f"unresolved reasons {dict(Counter(r[f'{s}_unresolved_reason'] for r in ok if r[f'{s}_unresolved_reason']))}")
    (OUT / "master_candidates_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf8")
    print("\n".join(lines))
    print(f"\nwrote {OUT / 'master_candidates.tsv'}")


if __name__ == "__main__":
    main()
