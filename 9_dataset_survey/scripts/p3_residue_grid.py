"""Phase 3: residue-by-residue map of every candidate onto the human reference frame (1MSV chain A).

Uses the same structural alignment as FOLDTRACE (tmtools TM-align on CA traces, loaded with
foldtrace's own structure loader), so positions agree with Marc's site calls.

For every candidate, every human position gets one call:
    identical     an aligned residue, within 4.0 A of the human position, same amino acid as human
    substituted   an aligned residue, within 4.0 A, different amino acid
    displaced     an aligned residue, but more than 4.0 A away (structure differs locally)
    deleted       no residue aligned, inside the part of the protein that aligns (internal gap)
    not_covered   no residue aligned, before the first / after the last aligned position
                  (a short or fragmentary model, NOT a deletion)
"Human amino acid" is wild-type UniProt P17707, so position 68 counts as Ser (1MSV holds S68A).

Between consecutive aligned positions, the script also records indel events:
    insertion     candidate has extra residues, human has none there
    deletion      human has residues there, candidate has none
    replacement   both have unaligned residues there (a loop that differs in both)
    N_extension / C_extension   candidate residues before / after the aligned part
Events spanning residues 22-27 (no coordinates in 1MSV) are flagged at_reference_gap, because
a "difference" there may only mean the crystal is missing that loop.

Every call and event carries the AlphaFold pLDDT of that region in that model: the residue's own
pLDDT, the mean over inserted residues, or for a deletion the mean of the two flanking residues.
confidence = high if that pLDDT >= 70, else low.

RESUMABLE: each alignment is cached in data/p3_alignments/<acc>.json; re-runs skip cached ones.
Reads (never writes) the FOLDTRACE run folder.

Writes
    results/p3_residue_grid.tsv       one row per (fold-passing candidate, human position)
    results/p3_indel_events.tsv       one row per insertion / deletion / replacement / extension
    results/p3_protein_summary.tsv    one row per candidate (all 1,000; fold_ok kept)
Run (a few minutes):
    py -3.12 p3_residue_grid.py
"""
import csv
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
from tmtools import tm_align
from foldtrace.io import load_structure

RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
BASE = Path(__file__).resolve().parent.parent
RESULTS, ALN, REFS = BASE / "results", BASE / "data" / "p3_alignments", BASE / "data" / "references"
REF_PDB, REF_CHAIN = RUN / "1MSV_A.pdb", "A"
OFFSET, TM_GATE, PLDDT_OK = 4.0, 0.5, 70.0


def human_sequence() -> str:
    p = REFS / "P17707.fasta"
    if not p.exists():
        REFS.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/P17707.fasta", timeout=60) as r:
            p.write_bytes(r.read())
    return "".join(l.strip() for l in p.read_text().splitlines() if not l.startswith(">"))


def plddt_by_resnum(path: Path) -> dict:
    out = {}
    with open(path, encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                out[int(ln[22:26])] = float(ln[60:66])
    return out


def accessions_in_order() -> list:
    accs = []
    with open(RUN / "foldseek_hits.m8", encoding="utf8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) == 21:
                accs.append(re.search(r"AF-([A-Z0-9]+)-F\d+", row[1]).group(1))
    return accs


def align(acc: str, ref, cand) -> dict:
    p = ALN / f"{acc}.json"
    if p.exists():
        return json.loads(p.read_text())
    r = tm_align(ref.ca, cand.ca, ref.seq, cand.seq)
    d = {"accession": acc, "tm_norm_ref": float(r.tm_norm_chain1), "tm_norm_cand": float(r.tm_norm_chain2),
         "rmsd": float(r.rmsd), "seqxA": r.seqxA, "seqyA": r.seqyA,
         "u": np.asarray(r.u).tolist(), "t": np.asarray(r.t).tolist()}
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(d))
    os.replace(tmp, p)
    return d


def conf(x):
    return "" if x is None else ("high" if x >= PLDDT_OK else "low")


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 1) if xs else None


def main() -> None:
    for d in (RESULTS, ALN):
        d.mkdir(parents=True, exist_ok=True)
    human = human_sequence()
    ref = load_structure(str(REF_PDB), chain=REF_CHAIN)
    mism = [f"{human[n - 1]}{n}{a}" for n, a in zip(ref.resnums, ref.seq) if human[n - 1] != a]
    print(f"reference {REF_PDB.name}: {len(ref.resnums)} residues {ref.resnums[0]}-{ref.resnums[-1]}; "
          f"differences from P17707: {mism}")
    ref_breaks = {(a, b) for a, b in zip(ref.resnums, ref.resnums[1:]) if b != a + 1}

    accs = accessions_in_order()
    if len(sys.argv) > 2 and sys.argv[1] == "--limit":
        accs = accs[: int(sys.argv[2])]
    gfh = open(RESULTS / "p3_residue_grid.tsv", "w", newline="", encoding="utf8")
    efh = open(RESULTS / "p3_indel_events.tsv", "w", newline="", encoding="utf8")
    sfh = open(RESULTS / "p3_protein_summary.tsv", "w", newline="", encoding="utf8")
    G = csv.writer(gfh, delimiter="\t")
    E = csv.writer(efh, delimiter="\t")
    S = csv.writer(sfh, delimiter="\t")
    G.writerow(["accession", "ref_pos", "human_aa", "call", "cand_resnum", "cand_aa", "ca_offset",
                "plddt", "plddt_source", "confidence"])
    E.writerow(["accession", "event", "after_ref_pos", "before_ref_pos", "ref_positions_missing",
                "n_ref_missing", "cand_start", "cand_end", "cand_len", "cand_seq", "plddt_mean", "plddt_min",
                "plddt_flank_mean", "confidence", "at_reference_gap"])
    S.writerow(["accession", "fold_ok", "tm_norm_ref", "tm_norm_cand", "rmsd", "cand_length",
                "n_aligned", "ref_coverage", "n_in_place", "seq_identity_in_place",
                "n_identical", "n_substituted", "n_displaced", "n_deleted", "n_not_covered",
                "n_diff_high", "n_diff_low", "n_insertions", "inserted_residues", "n_deletion_events",
                "n_replacements", "N_ext_len", "C_ext_len", "plddt_mean_aligned", "plddt_mean_whole"])

    for k, acc in enumerate(accs, 1):
        mpath = RUN / "models" / f"{acc}.pdb"
        cand = load_structure(str(mpath))
        pl = plddt_by_resnum(mpath)
        a = align(acc, ref, cand)
        u, t = np.asarray(a["u"]), np.asarray(a["t"])
        ref_on_cand = ref.ca @ u.T + t
        fold_ok = a["tm_norm_ref"] >= TM_GATE

        # walk the alignment
        pairs, i, j = [], 0, 0             # (ref_idx, cand_idx) aligned pairs, in order
        ref_unal, cand_unal = [], []        # per-gap bookkeeping
        segs = []                           # (ref_idxs_skipped, cand_idxs_skipped) between aligned pairs
        cur_r, cur_c = [], []
        for x, y in zip(a["seqxA"], a["seqyA"]):
            if x != "-" and y != "-":
                segs.append((cur_r, cur_c))
                cur_r, cur_c = [], []
                pairs.append((i, j))
            elif x != "-":
                cur_r.append(i)
            elif y != "-":
                cur_c.append(j)
            i += x != "-"
            j += y != "-"
        tail = (cur_r, cur_c)
        cr = cand.resnums

        def pld(idx):
            return pl.get(cr[idx])

        calls = {}
        evc = {"insertion": 0, "insertion_res": 0, "deletion": 0, "deletion_res": 0,
               "replacement": 0, "replacement_res": 0}
        if pairs:
            first_r, last_r = pairs[0][0], pairs[-1][0]
            for ri, ci in pairs:
                off = float(np.linalg.norm(ref_on_cand[ri] - cand.ca[ci]))
                h = human[ref.resnums[ri] - 1]
                c = "displaced" if off > OFFSET else ("identical" if cand.seq[ci] == h else "substituted")
                calls[ri] = (c, cr[ci], cand.seq[ci], round(off, 2), pld(ci), "residue")
            # gaps between aligned pairs: segs[n] precedes pairs[n]
            for n in range(1, len(pairs)):
                rgap, cgap = segs[n]
                left_r, left_c = pairs[n - 1]
                right_r, right_c = pairs[n]
                flank = mean([pld(left_c), pld(right_c)])
                for ri in rgap:
                    calls[ri] = ("deleted", "", "", "", flank, "flanks")
                lo, hi = ref.resnums[left_r], ref.resnums[right_r]
                n_missing = hi - lo - 1               # human positions between the anchors
                at_gap = any(lo <= p and q <= hi for p, q in ref_breaks)
                if n_missing == 0 and not cgap:
                    continue
                cp = [pld(c) for c in cgap]
                ev = "insertion" if n_missing == 0 else ("deletion" if not cgap else "replacement")
                q = mean(cp) if cgap else flank
                evc[ev] += 1
                evc[ev + "_res"] += len(cgap)
                E.writerow([acc, ev, lo, hi, f"{lo + 1}-{hi - 1}" if n_missing else "", n_missing,
                            cr[cgap[0]] if cgap else "", cr[cgap[-1]] if cgap else "", len(cgap),
                            "".join(cand.seq[c] for c in cgap), mean(cp) if cgap else "",
                            min(x for x in cp if x is not None) if cgap and any(x is not None for x in cp) else "",
                            flank, conf(q), at_gap])
            for ri in range(len(ref.resnums)):
                if ri < first_r or ri > last_r:
                    calls[ri] = ("not_covered", "", "", "", None, "")
            # terminal extensions
            for ev, idxs, anchor in (("N_extension", segs[0][1], pairs[0]), ("C_extension", tail[1], pairs[-1])):
                if idxs:
                    cp = [pld(c) for c in idxs]
                    pos = ref.resnums[anchor[0]]
                    E.writerow([acc, ev, "" if ev == "N_extension" else pos, pos if ev == "N_extension" else "",
                                "", 0, cr[idxs[0]], cr[idxs[-1]], len(idxs), "".join(cand.seq[c] for c in idxs),
                                mean(cp), min(cp), "", conf(mean(cp)), False])
        n_ext = len(segs[0][1]) if pairs else 0
        c_ext = len(tail[1]) if pairs else 0

        cnt = {c: 0 for c in ("identical", "substituted", "displaced", "deleted", "not_covered")}
        diff_hi = diff_lo = 0
        for ri in range(len(ref.resnums)):
            c, rn, aa, off, p, src = calls.get(ri, ("not_covered", "", "", "", None, ""))
            cnt[c] += 1
            if c in ("substituted", "displaced", "deleted"):
                if p is not None and p >= PLDDT_OK:
                    diff_hi += 1
                else:
                    diff_lo += 1
            if fold_ok:
                G.writerow([acc, ref.resnums[ri], human[ref.resnums[ri] - 1], c, rn, aa, off,
                            "" if p is None else round(p, 1), src, conf(p)])
        in_place = cnt["identical"] + cnt["substituted"]
        al_pl = [pld(ci) for _, ci in pairs]
        S.writerow([acc, fold_ok, round(a["tm_norm_ref"], 4), round(a["tm_norm_cand"], 4), round(a["rmsd"], 3),
                    len(cand.seq), len(pairs), round(len(pairs) / len(ref.resnums), 3), in_place,
                    round(cnt["identical"] / in_place, 3) if in_place else "",
                    cnt["identical"], cnt["substituted"], cnt["displaced"], cnt["deleted"], cnt["not_covered"],
                    diff_hi, diff_lo, evc["insertion"], evc["insertion_res"], evc["deletion"], evc["replacement"],
                    n_ext, c_ext, mean(al_pl), mean(list(pl.values()))])
        if k % 50 == 0 or k == len(accs):
            print(f"  {k}/{len(accs)}", flush=True)
    for fh in (gfh, efh, sfh):
        fh.close()
    print(f"wrote {RESULTS / 'p3_residue_grid.tsv'}\n      {RESULTS / 'p3_indel_events.tsv'}\n"
          f"      {RESULTS / 'p3_protein_summary.tsv'}")


if __name__ == "__main__":
    sys.exit(main())
