"""Job 1 (handoff 19 Nov): the fold-passers with two or more PF01536 (AdoMetDC) Pfam matches in one chain.

Question: is each one ONE AdoMetDC domain that Pfam split into pieces (annotation artefact), or TWO AdoMetDC
units fused in one chain (and if so, is one of them dead: a single-chain prozyme)?

Set: fold-passers whose Pfam layout has >= 2 AdoMetDC matches (PF01536/PF02675). That is the 81 "AdoMetDC
repeated" proteins (76 with two matches, 5 with three) plus 6 that also carry another Pfam domain; column
`arch_class` keeps them apart and the report's headline numbers use the 81. A0A8C5BDJ7 is excluded (its
UniProt sequence and model differ). Everything was fixed here before any result was read.

Method (same aligner and site rules as FOLDTRACE and p3: tmtools TM-align on CA traces, reference 1MSV
chain A, foldtrace's own structure loader when it is importable)

1. Units.
   pass 1   the whole chain is TM-aligned to 1MSV A. This is FOLDTRACE's own alignment, so its site calls
            are reproduced and checked against master_candidates.tsv (column whole_chain_matches_foldtrace).
            "Used" residues = aligned pairs within 5 A after superposition (TM-align's own ':' cutoff);
            unit 1 = first..last used residue.
   pass 2   every stretch of >= 80 residues before or after unit 1 is TM-aligned to 1MSV A on its own.
            A second unit exists if one scores TM >= 0.5 (normalised by 1MSV A, the FOLDTRACE fold gate).
   per Pfam each PF01536 match is also aligned on its own: TM-score (by 1MSV A and by its own length) and
            the span of human positions it covers. One domain split in two gives two partial matches that
            cover different halves of human numbering; two units give two matches that each cover most of it.
   class    two_units: pass-2 TM >= 0.5 | partial_second_copy: 0.3 <= pass-2 TM < 0.5 |
            one_domain: otherwise.
2. Sites per unit. Ser68, Ser229, His243, Cys82, Glu247 mapped on each unit's own alignment, exactly as
   FOLDTRACE does it (sites_used.tsv): no aligned residue or CA offset > 4.0 A -> unresolved; the expected
   residue -> retained; a residue listed as tolerated (Ser229 T, Glu247 Q) -> altered; anything else -> lost.
   Then the project's pLDDT rule: if that residue's pLDDT < 70 the call becomes unresolved (low pLDDT).
   (FOLDTRACE's own calls do not apply this rule; master_candidates.tsv has "lost" calls below 70.)
3. Dimer fit (two_units / partial_second_copy only). Template: the human dimer in 1MSV (chains A and B,
   unprocessed S68A proenzyme, PISA dimer; the same crystal FOLDTRACE uses). Processed-form check: 1JEN
   (chains B+A and D+C, processed PYR68, apo). Not 1I7B: per the repo's reference table its chains A and B are
   the alpha and beta chains of ONE protomer, so it holds no dimer. Unit residues are matched to template
   residues through their human numbering (pass-1/pass-2 alignments). Reported: each unit's own fit (RMSD);
   partner offset = superpose unit 1 on protomer 1, then RMSD of unit 2 to protomer 2 (best of the two ways
   round); pair TM-score over both units with one rigid superposition, normalised by the template dimer
   length; residues in contact between the two units (CA-CA < 8 A) vs the template.
   Thresholds, fixed now: partner offset <= 4 A packs like the dimer, 4-8 A similar but shifted, > 8 A different.
4. PAE (two_units / partial_second_copy only), AlphaFold DB. The model must be the same one: the AFDB
   sequence must equal the local model's, and the AFDB coordinates are compared with the local file. Mean PAE
   between the two units' used residues (both directions) and within each unit.
   Thresholds, fixed now: mean inter-unit PAE <= 10 A arrangement supported, 10-15 A uncertain, > 15 A
   unsupported (and nothing more is said about the arrangement).
5. Taxonomy (all) from data/taxonomy.json. Same species (two_units / partial only): other fold-passers of the
   same taxid in this run, and active UniProt entries for the species name with an AdoMetDC Pfam. Deleted
   entries are not searchable, so "none" means none among active entries.
6. The 16% site-loss confound. (a) Loss counted with the pLDDT rule, 81 vs 590 standalone (FOLDTRACE whole-
   chain calls, Fisher). (b) For every split point inside one domain: where it falls in human numbering, how
   many extra residues the protein has there (insertion) and their pLDDT, and whether lost sites sit inside
   the unmatched stretch. If splits are mostly insertions and losses are mostly low-pLDDT, the high rate is an
   annotation + model-confidence effect, not biology.

Reads (never writes) the FOLDTRACE run folder and lead 2's results/ and data/ files. Caches AlphaFold DB
files in data/q_afdb/ (with its own .gitignore, so the large PAE files are not committed) and UniProt
answers in data/q_uniprot/. Writes results_81/q2_*.tsv and results_81/q2_summary.txt.

Run (a few minutes; the network part is only the few two-unit proteins):
    py -3.12 q2_double_domain.py
"""
import csv
import json
import math
import os
import sys
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from tmtools import tm_align

BASE = Path(__file__).resolve().parent.parent
RESULTS, DATA, OUT = BASE / "results", BASE / "data", BASE / "results_81"
RUN = Path(os.environ.get("Q_RUN", r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search"
                                    r"\aim1_workflow\runs\stuy__adometdc"))
AFDB_DIR, UP_DIR = DATA / "q_afdb", DATA / "q_uniprot"
OFFLINE = os.environ.get("Q_OFFLINE") == "1"          # testing only: skip every web request

ADOMETDC = {"PF01536", "PF02675"}
EXCLUDE = {"A0A8C5BDJ7"}
SITES = [("Ser68", 68, "S", ""), ("Ser229", 229, "S", "T"), ("His243", 243, "H", ""),
         ("Cys82", 82, "C", ""), ("Glu247", 247, "E", "Q")]
OFFSET, TM_GATE, PLDDT_OK, CLOSE = 4.0, 0.5, 70.0, 5.0
MIN_SEG, PARTIAL = 80, 0.3
CONTACT = 8.0
INSERTION_MIN = 10

try:
    from foldtrace.io import load_structure as _ft_load
    LOADER = "foldtrace.io.load_structure"
except ImportError:                                   # only for testing outside Tui's environment
    _ft_load = None
    LOADER = "fallback PDB reader (foldtrace not importable)"


class Struct:
    def __init__(self, ca, seq, resnums):
        self.ca, self.seq, self.resnums = np.asarray(ca, float), seq, list(resnums)


THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
         "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
         "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}


def _fallback_load(path, chain=None):
    res, order = {}, []
    with open(path, encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ENDMDL"):
                break
            if not ln.startswith("ATOM") or ln[16] not in " A":
                continue
            ch = ln[21]
            if chain and ch != chain:
                continue
            key = (ch, int(ln[22:26]))
            if key not in res:
                res[key] = {"name": ln[17:20], "atoms": {}}
                order.append(key)
            res[key]["atoms"][ln[12:16].strip()] = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
    keep = [k for k in order if {"N", "CA", "C"} <= set(res[k]["atoms"]) and res[k]["name"] in THREE]
    return Struct([res[k]["atoms"]["CA"] for k in keep], "".join(THREE[res[k]["name"]] for k in keep),
                  [k[1] for k in keep])


def load(path, chain=None):
    if _ft_load:
        s = _ft_load(str(path), chain=chain) if chain else _ft_load(str(path))
        return Struct(s.ca, s.seq, s.resnums)
    return _fallback_load(path, chain)


def sub(s, idx):
    return Struct(s.ca[idx], "".join(s.seq[i] for i in idx), [s.resnums[i] for i in idx])


def plddt(path):
    out = {}
    with open(path, encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                out[int(ln[22:26])] = float(ln[60:66])
    return out


def read(name):
    with open(RESULTS / name, encoding="utf8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    if not rows:
        (OUT / name).write_text("(no rows)\n", encoding="utf8")
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT / name, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def fmt(x, n=3):
    return "" if x is None else (round(x, n) if isinstance(x, float) else x)


# ------------------------------------------------------------------ alignment
def align(ref, cand):
    """TM-align ref (1MSV A) onto cand, exactly as p3. Returns dict with correspondences in cand frame."""
    r = tm_align(ref.ca, cand.ca, ref.seq, cand.seq)
    u, t = np.asarray(r.u), np.asarray(r.t)
    ref_on = ref.ca @ u.T + t
    pairs, i, j = [], 0, 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-":
            pairs.append((i, j, float(np.linalg.norm(ref_on[i] - cand.ca[j]))))
        i += x != "-"
        j += y != "-"
    close = [(i, j, d) for i, j, d in pairs if d <= CLOSE]
    hum = [ref.resnums[i] for i, _, _ in close]
    return {"tm_ref": float(r.tm_norm_chain1), "tm_self": float(r.tm_norm_chain2), "rmsd": float(r.rmsd),
            "pairs": pairs, "close": close, "n_close": len(close),
            "human_from": min(hum) if hum else None, "human_to": max(hum) if hum else None,
            "ref_on": ref_on}


def call_sites(ref, cand, a, pl):
    by_ref = {i: (j, d) for i, j, d in a["pairs"]}
    out = {}
    for label, pos, expected, tolerated in SITES:
        if pos not in ref.resnums:
            out[label] = dict(state="unresolved", reason="not in reference", obs="", offset=None, plddt=None)
            continue
        ri = ref.resnums.index(pos)
        if ri not in by_ref:
            out[label] = dict(state="unresolved", reason="no aligned residue", obs="", offset=None, plddt=None)
            continue
        j, d = by_ref[ri]
        aa, num = cand.seq[j], cand.resnums[j]
        p = pl.get(num)
        if d > OFFSET:
            st, why = "unresolved", "offset > 4 A"
        else:
            st = "retained" if aa == expected else ("altered" if aa in tolerated else "lost")
            why = ""
        out[label] = dict(state=st, reason=why, obs=f"{aa}{num}", offset=d, plddt=p)
    for v in out.values():
        v["state_plddt_rule"] = (v["state"] if v["state"] == "unresolved" or (v["plddt"] or 0) >= PLDDT_OK
                                 else "unresolved")
        if v["state"] != "unresolved" and (v["plddt"] or 0) < PLDDT_OK:
            v["reason"] = f"pLDDT {v['plddt']:.1f} < 70 (FOLDTRACE call: {v['state']})"
    return out


def site_summary(calls, key="state_plddt_rule"):
    c = Counter(v[key] for v in calls.values())
    return f"{c['retained']} retained, {c['altered']} altered, {c['lost']} lost, {c['unresolved']} unresolved"


# ------------------------------------------------------------------ superposition helpers
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def rmsd(P, Q, R, t):
    return float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))


def tm_fixed(P, Q, L):
    """TM-score of fixed correspondences P->Q (template length L), superposition searched from several seeds."""
    d0 = max(0.5, 1.24 * (L - 15) ** (1 / 3) - 1.8)
    best = 0.0
    n = len(P)
    seeds = [np.arange(n)] + [np.arange(s, min(n, s + 20)) for s in range(0, max(1, n - 19), 10)]
    for idx in seeds:
        if len(idx) < 3:
            continue
        sel = idx
        for _ in range(20):
            R, t = kabsch(P[sel], Q[sel])
            d = np.linalg.norm(P @ R.T + t - Q, axis=1)
            tm = float((1 / (1 + (d / d0) ** 2)).sum() / L)
            best = max(best, tm)
            new = np.where(d < max(d0, 4.0))[0]
            if len(new) < 3 or (len(new) == len(sel) and np.all(new == sel)):
                break
            sel = new
    return best


def template(path, chain_sets):
    """Each protomer = the listed chains concatenated, keyed by human residue number."""
    prot = []
    for chains in chain_sets:
        d = {}
        for ch in chains:
            s = load(path, ch)
            for n, x in zip(s.resnums, s.ca):
                d[n] = x
        prot.append(d)
    return prot


def dimer_fit(units, ref, tmpl):
    """units: two (cand Struct, alignment) in sequence order. Returns metrics vs one template dimer."""
    def corr(cand, a, P):
        h = [(ref.resnums[i], cand.ca[j]) for i, j, _ in a["close"]]
        h = [(n, x) for n, x in h if n in P]
        return np.array([x for _, x in h]), np.array([P[n] for n, _ in h])
    L = len(tmpl[0]) + len(tmpl[1])
    best = None
    for p1, p2 in ((0, 1), (1, 0)):
        X1, Y1 = corr(*units[0], tmpl[p1])
        X2, Y2 = corr(*units[1], tmpl[p2])
        if len(X1) < 10 or len(X2) < 10:
            continue
        R, t = kabsch(X1, Y1)
        off = rmsd(X2, Y2, R, t)
        R2, t2 = kabsch(X2, Y2)
        res = {"n_pairs_unit1": len(X1), "n_pairs_unit2": len(X2),
               "unit1_own_rmsd": rmsd(X1, Y1, R, t), "unit2_own_rmsd": rmsd(X2, Y2, R2, t2),
               "partner_offset": off,
               "pair_tm": tm_fixed(np.vstack([X1, X2]), np.vstack([Y1, Y2]), L)}
        if best is None or res["partner_offset"] < best["partner_offset"]:
            best = res
    return best


def contacts(A, B):
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    return int((D < CONTACT).any(axis=1).sum())


# ------------------------------------------------------------------ web (AlphaFold DB, UniProt)
def http_get(url, binary=False):
    from p2_fetch_annotations import get           # same retry logic as lead 2's downloads
    st, _, body = get(url)
    return st, body


def afdb(acc):
    """(api_entry, pae_matrix, afdb_ca) from AlphaFold DB, cached in data/q_afdb/."""
    AFDB_DIR.mkdir(parents=True, exist_ok=True)
    gi = AFDB_DIR / ".gitignore"
    if not gi.exists():
        gi.write_text("# PAE matrices and models downloaded by q2_double_domain.py (large; re-downloadable)\n*\n"
                      "!.gitignore\n", encoding="utf8")
    api_p, pae_p, pdb_p = AFDB_DIR / f"{acc}_api.json", AFDB_DIR / f"{acc}_pae.npz", AFDB_DIR / f"{acc}.pdb"
    if not api_p.exists():
        st, body = http_get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
        if st != 200:
            return None, None, None
        api_p.write_text(body, encoding="utf8")
    entries = json.loads(api_p.read_text(encoding="utf8"))
    e = next((x for x in entries if x.get("uniprotAccession") == acc), entries[0] if entries else None)
    if e is None:
        return None, None, None
    if not pae_p.exists():
        st, body = http_get(e["paeDocUrl"])
        doc = json.loads(body)
        doc = doc[0] if isinstance(doc, list) else doc
        if "predicted_aligned_error" in doc:
            M = np.asarray(doc["predicted_aligned_error"], dtype=np.float32)
        else:                                          # old format: residue1, residue2, distance
            n = max(doc["residue1"])
            M = np.zeros((n, n), dtype=np.float32)
            M[np.asarray(doc["residue1"]) - 1, np.asarray(doc["residue2"]) - 1] = doc["distance"]
        np.savez_compressed(pae_p, pae=M.astype(np.float16))
    if not pdb_p.exists():
        st, body = http_get(e["pdbUrl"])
        pdb_p.write_text(body, encoding="utf8")
    return e, np.load(pae_p)["pae"].astype(np.float32), load(pdb_p)


def uniprot_species(species, taxid):
    UP_DIR.mkdir(parents=True, exist_ok=True)
    key = (species or taxid).replace(" ", "_").replace("/", "_")
    p = UP_DIR / f"{key}.tsv"
    if not p.exists():
        q = (f'taxonomy_name:"{species}"' if species else f"taxonomy_id:{taxid}") + \
            " AND (xref:pfam-PF01536 OR xref:pfam-PF02675)"
        st, body = http_get("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) +
                            "&format=tsv&fields=accession,length,xref_pfam,organism_name,reviewed,sequence&size=200")
        p.write_text(body if st == 200 else "", encoding="utf8")
    lines = p.read_text(encoding="utf8").strip().splitlines()[1:]
    return [l.split("\t") for l in lines if l.strip()]


def identity(a, b):
    """Percent identity over aligned columns, global alignment with free end gaps (as p7). None if no Biopython."""
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError:
        return None
    al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1, end_gap_score=0)
    x = al.align(a, b)[0]
    cols = [(a[i], b[j]) for (s1, e1), (s2, e2) in zip(*x.aligned) for i, j in zip(range(s1, e1), range(s2, e2))]
    return round(100 * sum(p == q for p, q in cols) / min(len(a), len(b)), 1) if cols else 0.0


def gene_model(acc, source, unit1_end, unit2_start):
    """ENA gene-model evidence and introns between the two units, with lead 2's p7 helpers (cached in data/ena/)."""
    try:
        import p7_fusion_validation as p7
        from p7_fusion_validation import sources, ena_record, parse_cds, exon_boundaries_protein, kind_of
        p7.ENA = DATA / "q_ena"                        # keep this job's ENA cache out of lead 2's data/ena/
    except Exception as ex:                            # noqa: BLE001
        return {"gene_model": f"p7 helpers not importable ({type(ex).__name__})"}
    srcs, _ = sources(acc, source)
    out = []
    for pid, mol in srcs:
        if not pid or pid == ".":
            continue
        rec = ena_record(pid)
        idl, loc, lens = parse_cds(rec)
        kind = kind_of(mol, idl)
        notes = sorted({l.split('"')[1] for l in rec.splitlines()
                        if "/note=\"source:" in l or "/note=\"Evidence" in l})
        b = exon_boundaries_protein(lens) if lens else []
        between = [x for x in b if unit1_end <= x <= unit2_start]
        out.append(f"{pid}: {kind}, {len(lens)} exons" + (f" ({', '.join(notes)})" if notes else "") +
                   (f"; introns between the units at residue {','.join(map(str, between))}" if between else
                    "; no intron between the units" if lens else ""))
    return {"gene_model": " | ".join(out) or "(no EMBL record)"}


# ------------------------------------------------------------------ main
def main():
    from p5_loss_vs_fusion import fisher            # lead 2's exact Fisher test
    master = {r["accession"]: r for r in read("master_candidates.tsv")}
    arch = {r["accession"]: r for r in read("architecture_per_protein.tsv")}
    tax = json.loads((DATA / "taxonomy.json").read_text(encoding="utf8"))
    ref = load(RUN / "1MSV_A.pdb", "A")
    tmpl_msv = template(RUN / "1MSV.pdb", [["A"], ["B"]])
    tmpl_jen = template(DATA / "references" / "1JEN.pdb", [["B", "A"], ["D", "C"]])
    msv_contacts = contacts(np.array(list(tmpl_msv[0].values())), np.array(list(tmpl_msv[1].values())))
    print(f"loader: {LOADER}; reference 1MSV A {len(ref.resnums)} residues; "
          f"1MSV dimer A/B contact residues (CA-CA < {CONTACT} A): {msv_contacts}")

    def ad_spans(a):
        return [(d.split()[0], *map(int, d.split()[1].split("-")))
                for d in arch[a]["domains_with_positions"].split("; ") if d and d.split()[0] in ADOMETDC]

    accs = [a for a, r in arch.items() if r["fold_ok"] == "true" and a not in EXCLUDE and len(ad_spans(a)) >= 2]
    only = [x for x in sys.argv[1:] if not x.startswith("-")]
    if only:
        accs = [a for a in accs if a in only]
    if "--available" in sys.argv:                      # testing: only models present locally
        accs = [a for a in accs if (RUN / "models" / f"{a}.pdb").exists()]
    print(f"{len(accs)} proteins with >= 2 AdoMetDC Pfam matches "
          f"({Counter(arch[a]['adometdc_order'] for a in accs)})")

    unit_rows, region_rows, site_rows, split_rows = [], [], [], []
    concord = Counter()
    for k, acc in enumerate(accs, 1):
        mp = RUN / "models" / f"{acc}.pdb"
        cand = load(mp)
        pl = plddt(mp)
        m, ar = master[acc], arch[acc]
        t = tax.get(m["taxid"], {"ranks": {}})
        rk = t.get("ranks", {})
        spans = ad_spans(acc)

        # pass 1: whole chain
        a1 = align(ref, cand)
        used = sorted(j for _, j, _ in a1["close"])
        u1_lo, u1_hi = cand.resnums[used[0]], cand.resnums[used[-1]]
        s1 = call_sites(ref, cand, a1, pl)
        cmap = {cand.resnums[j]: ref.resnums[i] for i, j, _ in a1["close"]}   # model residue -> human position
        for lab, _, _, _ in SITES:                     # reproduce FOLDTRACE
            ft_state, ft_obs = m[f"{lab}_state"], m[f"{lab}_obs"]
            ok = (s1[lab]["state"] == ft_state) and (ft_state == "unresolved" or s1[lab]["obs"] == ft_obs)
            concord["match" if ok else "differ"] += 1
            if not ok:
                print(f"  [check] {acc} {lab}: here {s1[lab]['state']} {s1[lab]['obs']}, "
                      f"FOLDTRACE {ft_state} {ft_obs}")

        # pass 2: stretches outside unit 1
        segs = []
        idx_before = [i for i, n in enumerate(cand.resnums) if n < u1_lo]
        idx_after = [i for i, n in enumerate(cand.resnums) if n > u1_hi]
        for side, idx in (("N-side", idx_before), ("C-side", idx_after)):
            if len(idx) >= MIN_SEG:
                sc = sub(cand, idx)
                a2 = align(ref, sc)
                segs.append((side, sc, a2))
        best2 = max(segs, key=lambda x: x[2]["tm_ref"]) if segs else None
        tm2 = best2[2]["tm_ref"] if best2 else None
        cls = ("two_units" if tm2 is not None and tm2 >= TM_GATE else
               "partial_second_copy" if tm2 is not None and tm2 >= PARTIAL else "one_domain")

        units = [("pass1", cand, a1, s1)]
        if cls != "one_domain":
            s2 = call_sites(ref, best2[1], best2[2], pl)
            units.append(("pass2 " + best2[0], best2[1], best2[2], s2))
        # order units by position in the chain
        units.sort(key=lambda u: min(u[1].resnums[j] for _, j, _ in u[2]["close"]) if u[2]["close"] else 0)
        for ui, (src, sc, a, s) in enumerate(units, 1):
            span = [sc.resnums[j] for _, j, _ in a["close"]]
            for lab, pos, _, _ in SITES:
                v = s[lab]
                site_rows.append({"accession": acc, "organism": m["organism"], "unit_class": cls, "unit": ui,
                                  "unit_source": src, "unit_span": f"{min(span)}-{max(span)}" if span else "",
                                  "site": lab, "obs": v["obs"], "ca_offset": fmt(v["offset"], 2),
                                  "plddt": fmt(v["plddt"], 1), "call_foldtrace_rule": v["state"],
                                  "call_with_plddt_rule": v["state_plddt_rule"], "note": v["reason"]})

        # Pfam regions on their own
        for n_, (pf, s_, e_) in enumerate(spans, 1):
            idx = [i for i, n in enumerate(cand.resnums) if s_ <= n <= e_]
            if len(idx) < 10:
                continue
            ar_ = align(ref, sub(cand, idx))
            in_unit = [ui for ui, (_, sc, a, _) in enumerate(units, 1)
                       if a["close"] and min(sc.resnums[j] for _, j, _ in a["close"]) <= e_
                       and max(sc.resnums[j] for _, j, _ in a["close"]) >= s_]
            region_rows.append({"accession": acc, "match": n_, "pfam": pf, "start": s_, "end": e_, "length": e_ - s_ + 1,
                                "tm_norm_1MSV_A": fmt(ar_["tm_ref"]), "tm_norm_region": fmt(ar_["tm_self"]),
                                "n_close_pairs": ar_["n_close"],
                                "human_positions_covered_alone": f"{ar_['human_from']}-{ar_['human_to']}" if ar_["close"] else "",
                                "human_positions_in_whole_chain_alignment": (lambda h: f"{min(h)}-{max(h)} ({len(h)} pos.)" if h else "")(
                                    [cmap[n] for n in range(s_, e_ + 1) if n in cmap]),
                                "plddt_mean": fmt(float(np.mean([pl[n] for n in range(s_, e_ + 1) if n in pl])), 1),
                                "in_unit": ",".join(map(str, in_unit))})

        # split points between consecutive Pfam matches inside the same unit (pass-1 frame)
        cl = sorted(cmap)
        for (pfa, sa, ea), (pfb, sb, eb) in zip(spans, spans[1:]):
            if not (u1_lo <= ea and sb <= u1_hi):
                continue                               # the gap is between units, not a split
            left = [n for n in cl if n <= ea]
            right = [n for n in cl if n >= sb]
            if not left or not right:
                continue
            a_, b_ = left[-1], right[0]
            ins = (b_ - a_ - 1) - (cmap[b_] - cmap[a_] - 1)
            gap = [n for n in range(ea + 1, sb)]
            lost_in_gap = [lab for lab, _, _, _ in SITES
                           if s1[lab]["state"] == "lost" and s1[lab]["obs"] and ea < int(s1[lab]["obs"][1:]) < sb]
            split_rows.append({"accession": acc, "organism": m["organism"], "arch_class": ar["adometdc_order"],
                               "match_left_end": ea, "match_right_start": sb, "gap_residues": max(0, sb - ea - 1),
                               "anchor_left": f"{a_}~human {cmap[a_]}", "anchor_right": f"{b_}~human {cmap[b_]}",
                               "extra_residues_vs_human": ins,
                               "split_cause": "insertion" if ins >= INSERTION_MIN else "divergent stretch (no insertion)",
                               "gap_frac_aligned": fmt(sum(n in cmap for n in gap) / len(gap), 2) if gap else "",
                               "gap_plddt_mean": fmt(float(np.mean([pl[n] for n in gap if n in pl])), 1) if gap else "",
                               "foldtrace_lost_sites_in_gap": ",".join(lost_in_gap)})

        row = {"accession": acc, "organism": m["organism"], "taxid": m["taxid"],
               "kingdom": rk.get("kingdom", ""), "phylum": rk.get("phylum", ""), "class": rk.get("class", ""),
               "order": rk.get("order", ""), "species": rk.get("species", ""),
               "arch_class": ar["adometdc_order"], "architecture": ar["arch_pfam_ids"], "source": ar["source"],
               "uniprot_status": ar["uniprot_status"], "length": len(cand.seq),
               "adometdc_matches": "; ".join(f"{s_}-{e_}" for _, s_, e_ in spans),
               "whole_chain_tm": fmt(a1["tm_ref"]), "unit1_span": f"{u1_lo}-{u1_hi}",
               "unit1_human_covered": f"{a1['human_from']}-{a1['human_to']}",
               "pass2_segments": "; ".join(f"{sd} {sc.resnums[0]}-{sc.resnums[-1]} TM {a['tm_ref']:.3f}"
                                            for sd, sc, a in segs) or "(none >= 80 aa)",
               "pass2_best_tm": fmt(tm2), "unit_class": cls, "n_units": len(units),
               "foldtrace_verdict": m["verdict"],
               "foldtrace_lost_sites": ",".join(l for l, _, _, _ in SITES if m[f"{l}_state"] == "lost"),
               "lost_with_plddt_rule": ",".join(l for l, _, _, _ in SITES if s1[l]["state_plddt_rule"] == "lost")}
        for ui, (_, sc, a, s) in enumerate(units, 1):
            row[f"unit{ui}_sites"] = site_summary(s)

        # dimer fit + PAE + same species: only where a second unit exists
        if cls != "one_domain":
            pair = [(u[1], u[2]) for u in units[:2]]
            uspan = [(min(sc.resnums[j] for _, j, _ in a["close"]), max(sc.resnums[j] for _, j, _ in a["close"]))
                     for sc, a in pair]
            useq = [cand.seq[cand.resnums.index(lo):cand.resnums.index(hi) + 1] for lo, hi in uspan]
            row["unit_spans"] = "; ".join(f"{lo}-{hi}" for lo, hi in uspan)
            row["unit1_vs_unit2_identity"] = identity(*useq)
            for name, tm_ in (("1MSV", tmpl_msv), ("1JEN", tmpl_jen)):
                f = dimer_fit(pair, ref, tm_)
                if f:
                    row.update({f"{name}_{kk}": fmt(v, 2) for kk, v in f.items()})
            UA = np.array([pair[0][0].ca[j] for _, j, _ in pair[0][1]["close"]])
            UB = np.array([pair[1][0].ca[j] for _, j, _ in pair[1][1]["close"]])
            row["model_contact_residues_unit1"] = contacts(UA, UB)
            row["template_contact_residues"] = msv_contacts
            if not OFFLINE:
                try:
                    e, M, afs = afdb(acc)
                    if e is None:
                        row["pae_note"] = "not in AlphaFold DB"
                    else:
                        same_seq = e.get("sequence", e.get("uniprotSequence")) == cand.seq
                        common = [n for n in cand.resnums if n in afs.resnums]
                        dev = (float(np.sqrt(((cand.ca[[cand.resnums.index(n) for n in common]] -
                                               afs.ca[[afs.resnums.index(n) for n in common]]) ** 2).sum(1).mean()))
                               if common else None)
                        ia = [pair[0][0].resnums[j] - 1 for _, j, _ in pair[0][1]["close"]]
                        ib = [pair[1][0].resnums[j] - 1 for _, j, _ in pair[1][1]["close"]]
                        inter = np.concatenate([M[np.ix_(ia, ib)].ravel(), M[np.ix_(ib, ia)].ravel()])
                        mi = float(inter.mean())
                        row.update({"afdb_version": e.get("latestVersion", ""), "afdb_same_sequence": same_seq,
                                    "afdb_vs_local_ca_rmsd": fmt(dev, 2),
                                    "pae_inter_mean": fmt(mi, 1), "pae_inter_median": fmt(float(np.median(inter)), 1),
                                    "pae_inter_frac_below_10": fmt(float((inter < 10).mean()), 2),
                                    "pae_intra_unit1_mean": fmt(float(M[np.ix_(ia, ia)].mean()), 1),
                                    "pae_intra_unit2_mean": fmt(float(M[np.ix_(ib, ib)].mean()), 1),
                                    "arrangement": ("model mismatch: PAE not used" if not same_seq or (dev or 0) > 1.0
                                                    else "supported" if mi <= 10 else
                                                    "uncertain" if mi <= 15 else "unsupported")})
                except Exception as ex:                # noqa: BLE001  keep going; report the failure
                    row["pae_note"] = f"AFDB failed: {type(ex).__name__}: {ex}"[:200]
                try:
                    hits = uniprot_species(rk.get("species", ""), m["taxid"])
                    row["uniprot_species_adometdc"] = "; ".join(
                        f"{h[0]} ({h[1]} aa, {arch[h[0]]['adometdc_order'] if h[0] in arch else 'not in run'}"
                        + (f"; identity to unit 1 {identity(useq[0], h[5])}%, unit 2 {identity(useq[1], h[5])}%"
                           if len(h) > 5 and h[5] else "") + ")"
                        for h in hits if h[0] != acc) or "(no other active entry)"
                except Exception as ex:                # noqa: BLE001
                    row["uniprot_species_adometdc"] = f"UniProt failed: {type(ex).__name__}"
                try:
                    row.update(gene_model(acc, ar["source"], uspan[0][1], uspan[1][0]))
                except Exception as ex:                # noqa: BLE001
                    row["gene_model"] = f"ENA failed: {type(ex).__name__}: {ex}"[:200]
            same = [b for b, r in master.items() if b != acc and r["taxid"] == m["taxid"] and r["fold_ok"] == "true"]
            row["same_taxid_in_run"] = "; ".join(f"{b} ({arch[b]['adometdc_order']}, {arch[b]['length']} aa)"
                                                 for b in same) or "(none)"
        unit_rows.append(row)
        if k % 10 == 0 or k == len(accs):
            print(f"  {k}/{len(accs)}", flush=True)

    write("q2_units.tsv", unit_rows)
    write("q2_pfam_regions.tsv", region_rows)
    write("q2_sites_per_unit.tsv", site_rows)
    write("q2_splits.tsv", split_rows)

    # ------------------------------------------------------------------ summary
    L = []
    rep = [r for r in unit_rows if r["arch_class"] == "AdoMetDC repeated"]
    L.append(f"q2_double_domain: {len(unit_rows)} proteins ({len(rep)} 'AdoMetDC repeated'); loader {LOADER}")
    L.append(f"whole-chain site calls reproduce FOLDTRACE: {concord['match']} match, {concord['differ']} differ")
    L.append(f"unit class (all): {dict(Counter(r['unit_class'] for r in unit_rows))}")
    L.append(f"unit class (the 81): {dict(Counter(r['unit_class'] for r in rep))}")
    L.append("pass-2 best TM (all proteins with a >= 80 aa stretch outside unit 1): " +
             ", ".join(f"{r['accession']} {r['pass2_best_tm']}" for r in sorted(
                 [r for r in unit_rows if r["pass2_best_tm"] != ""], key=lambda r: -r["pass2_best_tm"])[:15]))
    L.append("\nproteins with a second unit:")
    for r in unit_rows:
        if r["unit_class"] != "one_domain":
            keys = ["accession", "organism", "phylum", "arch_class", "length", "adometdc_matches", "whole_chain_tm",
                    "pass2_segments", "unit_class", "unit_spans", "unit1_vs_unit2_identity", "unit1_sites", "unit2_sites", "1MSV_partner_offset",
                    "1MSV_pair_tm", "1JEN_partner_offset", "1JEN_pair_tm", "model_contact_residues_unit1",
                    "afdb_same_sequence", "afdb_vs_local_ca_rmsd", "pae_inter_mean", "pae_intra_unit1_mean",
                    "pae_intra_unit2_mean", "arrangement", "pae_note", "gene_model", "same_taxid_in_run",
                    "uniprot_species_adometdc"]
            L.append("  " + " | ".join(f"{x}={r.get(x, '')}" for x in keys))

    # Pfam match coverage, one_domain proteins of the 81
    regs = defaultdict(list)
    for x in region_rows:
        regs[x["accession"]].append(x)
    od = [r for r in rep if r["unit_class"] == "one_domain"]
    n_any = sum(any((x["tm_norm_1MSV_A"] or 0) >= TM_GATE for x in regs[r["accession"]]) for r in od)
    n_two = sum(sum((x["tm_norm_1MSV_A"] or 0) >= TM_GATE for x in regs[r["accession"]]) >= 2 for r in od)
    L.append(f"\none_domain among the 81: {len(od)}; with >= 1 Pfam match scoring TM >= 0.5 on its own: {n_any}; "
             f"with >= 2: {n_two}")

    # splits
    sp81 = [s for s in split_rows if s["arch_class"] == "AdoMetDC repeated"]
    L.append(f"\nsplit points inside one domain (the 81): {len(sp81)}; cause: "
             f"{dict(Counter(s['split_cause'] for s in sp81))}")
    pos = Counter()
    for s in sp81:
        h = int(s["anchor_left"].split("human ")[1])
        pos[(h // 25) * 25] += 1
    L.append("  left anchor in human numbering, 25-residue bins: " +
             ", ".join(f"{b}-{b + 24}: {n}" for b, n in sorted(pos.items())))
    ins = [s["extra_residues_vs_human"] for s in sp81]
    if ins:
        L.append(f"  extra residues vs human at the split: median {sorted(ins)[len(ins) // 2]}, "
                 f"range {min(ins)} to {max(ins)}")
    king = {r["accession"]: r["kingdom"] or r["phylum"] or "(unknown)" for r in unit_rows}
    L.append("  the 81 by kingdom: " + ", ".join(f"{k}: {n}" for k, n in Counter(king[r['accession']] for r in rep).most_common()))
    bycause = defaultdict(Counter)
    for s in sp81:
        bycause[king[s["accession"]]][s["split_cause"]] += 1
    L.append("  split cause by kingdom: " + "; ".join(f"{k}: {dict(v)}" for k, v in bycause.items()))
    L.append(f"  lost FOLDTRACE sites lying inside an unmatched stretch: "
             f"{sum(bool(s['foldtrace_lost_sites_in_gap']) for s in sp81)} split points")

    # confound (a): loss with the pLDDT rule, 81 vs standalone, FOLDTRACE whole-chain calls
    def loss(r, rule):
        for lab, _, _, _ in SITES:
            if r[f"{lab}_state"] == "lost":
                p = float(r[f"{lab}_plddt"]) if r[f"{lab}_plddt"] else 0
                if not rule or p >= PLDDT_OK:
                    return True
        return False
    stand = [a for a, r in arch.items() if r["fold_ok"] == "true" and r["adometdc_order"] == "AdoMetDC only"
             and a not in EXCLUDE]
    rep_acc = [a for a, r in arch.items() if r["fold_ok"] == "true" and r["adometdc_order"] == "AdoMetDC repeated"
               and a not in EXCLUDE]
    for rule in (False, True):
        a_ = sum(loss(master[a], rule) for a in rep_acc)
        c_ = sum(loss(master[a], rule) for a in stand)
        p = fisher(a_, len(rep_acc) - a_, c_, len(stand) - c_)
        L.append(f"\nloss of >= 1 site, {'only calls at pLDDT >= 70' if rule else 'all FOLDTRACE calls'}: "
                 f"repeated {a_}/{len(rep_acc)} ({100 * a_ / len(rep_acc):.1f}%) vs standalone {c_}/{len(stand)} "
                 f"({100 * c_ / len(stand):.1f}%), Fisher p = {p:.3g}")
    L.append("  repeated proteins with a loss at pLDDT >= 70: " + "; ".join(
        f"{a} {master[a]['organism']} ({','.join(l for l, _, _, _ in SITES if master[a][f'{l}_state'] == 'lost' and float(master[a][f'{l}_plddt'] or 0) >= PLDDT_OK)})"
        for a in rep_acc if loss(master[a], True)))
    txt = "\n".join(L) + "\n"
    (OUT / "q2_summary.txt").write_text(txt, encoding="utf8")
    print("\n" + txt)
    print(f"wrote {OUT}\\q2_units.tsv, q2_pfam_regions.tsv, q2_sites_per_unit.tsv, q2_splits.tsv, q2_summary.txt")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
