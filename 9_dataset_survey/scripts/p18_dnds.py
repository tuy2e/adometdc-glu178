"""Track B p18: is the degraded ODC half under purifying selection, or decaying? dN/dS per region, per pair.

Candidates (all from results/p10_fusion_search.tsv / p10_all_fusions.tsv - nothing new is searched):
  Galdieria      G. partita A0A9C7PQ51 + A0A9C7PSW2, G. sulphuraria M2YAF1 + UPI00207E5336, G. yellowstonensis A0AAV9I966
  Gracilariaceae Gracilariopsis chorda A0A2V3IZY9, Gracilaria gracilis UPI00403B59D6, Agarophyton chilense UPI00403BC8E0
  Ceratopteris richardii, 8 fusion copies (green lineage, ODC half intact = positive control)
  Rhodosorus marinus A0AAV8UZF7 and Percolomonas cosmopolitus A0A7S1KLJ4: no congener fusion in p10 (no pair)
Pairs = every pair within a genus / family group above.

Step 1  protein identity per pair and per region (pairwise FAMSA alignment, identity = identical / gap-free columns).
Step 2  CDS: UniParc cross-references (rest.uniprot.org/uniparc/<UPI>/databases) -> ENA protein_id -> ENA coding
        record (ebi.ac.uk/ena/browser/api/embl/<protein_id>; its sequence block is the spliced CDS). Cached in
        data/p18_cds/. A CDS is kept only if its translation (standard code, /codon_start honoured, final stop
        removed) equals the UniParc protein exactly; otherwise it is dropped with the reason.
Step 3  codon alignment: protein pair aligned with FAMSA (pyfamsa; MAFFT is not installed), codons threaded onto
        it. Columns with a gap, stop or ambiguous base in either sequence are removed. dN, dS by Nei-Gojobori 1986
        with Jukes-Cantor correction (Biopython cal_dn_ds, method NG86), separately for
          AdoMetDC  first PF01536 start .. last PF01536 end
          linker    after the AdoMetDC Pfam end .. before the first PLP-decarboxylase Pfam start (PF02784/PF00278)
          ODC       first PLP Pfam start .. C-terminus
        Regions are read on the first sequence of each pair (Pfam spans per protein from p10_fusion_search.tsv;
        results/fusion_architectures.tsv has no per-protein positions). The ODC region runs to the C-terminus
        because the Pfam PF02784 match is truncated in several red-algal fusions (G. partita 434-512 only).
        dS > 1, or dS not computable (proportion of synonymous differences too high for the JC correction),
        = saturated: dN/dS is written but marked not interpretable.

Writes data/p18_cds/*, results/p18_dnds.tsv, results/p18_report.md
Run:  py -3.12 p18_dnds.py
"""
import csv
import itertools
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request
import warnings

from Bio import BiopythonExperimentalWarning
from Bio.Seq import Seq

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from pyfamsa import Aligner, Sequence

from p2_fetch_annotations import DATA, RESULTS, USER_AGENT

sys.stdout.reconfigure(encoding="utf8")
CDS = DATA / "p18_cds"
CDS.mkdir(parents=True, exist_ok=True)
GROUPS = {
    "Galdieria": ["A0A9C7PQ51", "A0A9C7PSW2", "M2YAF1", "UPI00207E5336", "A0AAV9I966"],
    "Gracilariaceae": ["A0A2V3IZY9", "UPI00403B59D6", "UPI00403BC8E0"],
    "Ceratopteris (control)": ["A0A8T2T6R3", "A0A8T2T887", "A0A8T2TAW1", "A0A8T2T5U0", "A0A8T2V5S6",
                               "A0A8T2V9N7", "A0A8T2V9B8", "A0A8T2VH64"],
}
SINGLETONS = {"A0AAV8UZF7": "Rhodosorus marinus", "A0A7S1KLJ4": "Percolomonas cosmopolitus"}
ADO, PLP = {"PF01536"}, {"PF02784", "PF00278"}
STOPS = {"TAA", "TAG", "TGA"}
LOW_DS = 0.05   # below this, too few synonymous changes for a ratio to mean anything (threshold chosen for p18)


def http(url, accept="application/json"):
    for k in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode("utf8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return ""
            time.sleep(2 ** k)
        except OSError:
            time.sleep(2 ** k)
    raise RuntimeError(url)


def cached(name, url, accept="application/json"):
    p = CDS / name
    if not p.exists():
        p.write_text(http(url, accept), encoding="utf8")
        time.sleep(0.4)
    return p.read_text(encoding="utf8")


# ---------------------------------------------------------------- candidates
def candidates():
    out = {}
    with open(RESULTS / "p10_fusion_search.tsv", encoding="utf8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            for key in [r["upi"]] + [a.split(".")[0] for a in r["uniprot_accessions"].split(",") if a]:
                out[key] = r
    return out


def regions(arch, length):
    doms = [(d.split()[0], *map(int, d.split()[1].split("-"))) for d in arch.split(" > ")]
    a0 = min(s for p, s, e in doms if p in ADO)
    a1 = max(e for p, s, e in doms if p in ADO)
    p0 = min(s for p, s, e in doms if p in PLP and s > a1)
    return {"AdoMetDC": (a0, a1), "linker": (a1 + 1, p0 - 1), "ODC": (p0, length)}


# ---------------------------------------------------------------- CDS
def parse_coding(rec):
    m = re.search(r"^SQ .*?\n(.*?)^//", rec, re.M | re.S)
    if not m:
        return None, None
    nt = re.sub(r"[^acgtnACGTN]", "", re.sub(r"\d", "", m.group(1))).upper()
    cs = re.search(r'/codon_start=(\d)', rec)
    return nt, int(cs.group(1)) if cs else 1


def cds_for(upi, prot):
    """(cds or None, note). Tries each ENA protein_id cross-reference of the UniParc entry until one translates exactly."""
    xr = json.loads(cached(f"{upi}.uniparc_xrefs.json", f"https://rest.uniprot.org/uniparc/{upi}/databases?format=json&size=100")
                    or '{"results": []}')
    pids = []
    for x in xr["results"]:
        if x.get("database", "").startswith("EMBL") and x.get("id"):
            pid = f"{x['id']}.{x['version']}" if x.get("version") and "." not in x["id"] else x["id"]
            if pid not in pids:
                pids.append(pid)
    if not pids:
        return None, "no ENA protein_id among the UniParc cross-references"
    notes = []
    for pid in pids:
        rec = cached(f"{pid}.embl", f"https://www.ebi.ac.uk/ena/browser/api/embl/{pid}", accept="text/plain")
        nt, start = parse_coding(rec)
        if not nt:
            notes.append(f"{pid}: no sequence in ENA record")
            continue
        nt = nt[start - 1:]
        nt = nt[:len(nt) // 3 * 3]
        tr = str(Seq(nt).translate()).rstrip("*")
        if nt[-3:] in STOPS:
            nt = nt[:-3]
        if tr == prot:
            return nt, f"{pid} (translation identical to UniParc)"
        notes.append(f"{pid}: translation differs ({len(tr)} aa vs {len(prot)} aa; "
                     f"{sum(a != b for a, b in zip(tr, prot))} mismatches over the shared length)")
    return None, "; ".join(notes)


# ---------------------------------------------------------------- alignment and dN/dS
def align(s1, s2):
    msa = Aligner(guide_tree="upgma").align([Sequence(b"a", s1.encode()), Sequence(b"b", s2.encode())])
    d = {s.id.decode(): s.sequence.decode() for s in msa}
    return d["a"], d["b"]


def columns(a1, a2):
    """[(pos1, pos2, aa1, aa2)] with 1-based residue numbers, None in a gap."""
    out, i, j = [], 0, 0
    for x, y in zip(a1, a2):
        i += x != "-"
        j += y != "-"
        out.append((i if x != "-" else None, j if y != "-" else None, x, y))
    return out


def ng86(c1, c2):
    if not c1:
        return None, None
    try:
        dn, ds = cal_dn_ds(CodonSeq("".join(c1)), CodonSeq("".join(c2)), method="NG86")
    except (ValueError, ZeroDivisionError):
        return None, None
    ok = lambda v: v is not None and not math.isnan(v) and v >= 0     # Biopython returns -1 when not computable
    return (dn if ok(dn) else None), (ds if ok(ds) else None)


def pair_type(x, y, meta):
    if meta[x]["organism"] != meta[y]["organism"]:
        return "different species"
    if meta[x].get("locus") and meta[x]["locus"] == meta[y]["locus"]:
        return "isoforms of one locus (not an independent pair)"
    if meta[x].get("strain") != meta[y].get("strain"):
        return "same species, different strains (allelic)"
    return "same genome, different loci (paralogs)"


def r3(x):
    return "" if x is None else round(x, 3)


def main():
    cand = candidates()
    seqs, cds, cds_note, meta = {}, {}, {}, {}
    for acc in sorted({a for g in GROUPS.values() for a in g} | set(SINGLETONS)):
        r = cand[acc]
        upi = r["upi"]
        ent = json.loads(cached(f"{upi}.uniparc.json", f"https://rest.uniprot.org/uniparc/{upi}?format=json"))
        seqs[acc] = ent["sequence"]["value"]
        meta[acc] = {"organism": r["organism"], "upi": upi, "arch": r["pfam_architecture"],
                     "regions": regions(r["pfam_architecture"], len(seqs[acc]))}
        if acc in SINGLETONS:
            continue
        cds[acc], cds_note[acc] = cds_for(upi, seqs[acc])
        if cds[acc]:
            rec = (CDS / f"{cds_note[acc].split()[0]}.embl").read_text(encoding="utf8")
            for k in ("locus_tag", "strain"):
                m = re.search(rf'/{k}="([^"]+)"', rec)
                meta[acc][k.replace("_tag", "")] = m.group(1) if m else ""
        print(f"{acc:14s} {r['organism']:28s} CDS {'kept' if cds[acc] else 'DROPPED'}: {cds_note[acc]}", flush=True)

    rows = []
    for group, accs in GROUPS.items():
        for x, y in itertools.combinations(accs, 2):
            a1, a2 = align(seqs[x], seqs[y])
            cols = columns(a1, a2)
            reg = meta[x]["regions"]
            row = {"group": group, "seq1": x, "organism1": meta[x]["organism"], "seq2": y,
                   "organism2": meta[y]["organism"],
                   "pair_type": pair_type(x, y, meta),
                   "cds1": "kept" if cds.get(x) else "dropped", "cds2": "kept" if cds.get(y) else "dropped"}
            for name, (lo, hi) in [("all", (1, len(seqs[x])))] + list(reg.items()):
                cc = [c for c in cols if c[0] is not None and lo <= c[0] <= hi and c[1] is not None]
                row[f"{name}_aligned_residues"] = len(cc)
                row[f"{name}_pct_identity"] = round(100 * sum(c[2] == c[3] for c in cc) / len(cc), 1) if cc else ""
                if name == "all" or not (cds.get(x) and cds.get(y)):
                    continue
                c1, c2 = [], []
                for p1, p2, _, _ in cc:
                    k1, k2 = cds[x][3 * (p1 - 1):3 * p1], cds[y][3 * (p2 - 1):3 * p2]
                    if len(k1) == len(k2) == 3 and not ({k1, k2} & STOPS) and set(k1 + k2) <= set("ACGT"):
                        c1.append(k1)
                        c2.append(k2)
                dn, ds = ng86(c1, c2)
                sat = ds is None or ds > 1
                low = ds is not None and ds < LOW_DS
                row[f"{name}_codons"] = len(c1)
                row[f"{name}_dN"] = r3(dn)
                row[f"{name}_dS"] = r3(ds) if ds is not None else "not computable"
                row[f"{name}_dNdS"] = r3(dn / ds) if dn is not None and ds and ds > 0 else ""
                row[f"{name}_differing_codons"] = sum(a != b for a, b in zip(c1, c2))
                row[f"{name}_saturated"] = sat
                row[f"{name}_too_little_divergence"] = low
                row[f"{name}_interpretable"] = not sat and not low
            rows.append(row)
            print(f"{group:24s} {x} vs {y}: identity {row['all_pct_identity']} %", flush=True)
    cols_ = []
    for r in rows:
        cols_ += [c for c in r if c not in cols_]
    with open(RESULTS / "p18_dnds.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, cols_, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    grac_ok = any(r.get("ODC_interpretable") is True for r in rows if "A0A2V3IZY9" in (r["seq1"], r["seq2"]))
    fb = [] if grac_ok else fallback(meta)
    report(rows, cds_note, meta, fb, grac_ok)
    json.dump({"rows": rows, "cds_note": cds_note, "meta": {k: {**v, "regions": {a: list(b) for a, b in v["regions"].items()}}
                                                               for k, v in meta.items()}},
              open(CDS / "_p18_summary.json", "w", encoding="utf8"), indent=1)


# ---------------------------------------------------------------- step 4 fallback (protein level, weaker)
def fallback(meta):
    """Per-region identity of each fusion to the consensus of all p10 AdoMetDC-then-PLP fusions (that region aligned
    with FAMSA; consensus = majority residue in columns with <= 50 % gaps; the protein itself is included)."""
    ents = json.loads(next((DATA / "p10").glob("uniparc__feature_id_PF01536*.json")).read_text(encoding="utf8"))
    seq_of = {e["uniParcId"]: e["sequence"]["value"] for e in ents}
    fam = {}
    with open(RESULTS / "p10_fusion_search.tsv", encoding="utf8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["adometdc_then_plp"] == "True" and r["upi"] in seq_of:
                try:
                    fam[r["upi"]] = (seq_of[r["upi"]], regions(r["pfam_architecture"], len(seq_of[r["upi"]])))
                except ValueError:
                    pass
    targets = [a for g in GROUPS.values() for a in g] + list(SINGLETONS)
    keep_cer = {"A0A8T2T6R3", "A0A8T2V5S6"}                       # one isoform per Ceratopteris locus
    targets = [a for a in targets if not a.startswith("A0A8T2") or a in keep_cer]
    out = {a: {"accession": a, "organism": meta[a]["organism"], "n_family": len(fam)} for a in targets}
    for name in ("AdoMetDC", "linker", "ODC"):
        ids = [u for u in fam if fam[u][1][name][1] >= fam[u][1][name][0]]
        msa = Aligner(guide_tree="upgma").align(
            [Sequence(u.encode(), fam[u][0][fam[u][1][name][0] - 1:fam[u][1][name][1]].encode()) for u in ids])
        aln = {x.id.decode(): x.sequence.decode() for x in msa}
        L = len(next(iter(aln.values())))
        cons = []
        for k in range(L):
            res = [aln[u][k] for u in ids if aln[u][k] != "-"]
            cons.append(max(set(res), key=res.count) if len(res) >= len(ids) / 2 else None)
        for a in targets:
            t = aln.get(meta[a]["upi"])
            if t is None:
                out[a][f"{name}_pct_identity_to_consensus"], out[a][f"{name}_columns"] = "", 0
                continue
            pairs = [(c, t[k]) for k, c in enumerate(cons) if c and t[k] != "-"]
            out[a][f"{name}_pct_identity_to_consensus"] = round(100 * sum(c == x for c, x in pairs) / len(pairs), 1) if pairs else ""
            out[a][f"{name}_columns"] = len(pairs)
    return list(out.values())


# ---------------------------------------------------------------- report
def fmt(r, k):
    return "" if r.get(k) in (None, "") else str(r[k])


REP_CER = {"A0A8T2T6R3", "A0A8T2V5S6"}          # one isoform per Ceratopteris locus in the report tables


def shown(r):
    return all(not a.startswith("A0A8T2") or a in REP_CER for a in (r["seq1"], r["seq2"]))


def report(rows, cds_note, meta, fb, grac_ok):
    L = ["# p18: is the degraded ODC half under purifying selection, or decaying?", "",
         "Generated by `scripts/p18_dnds.py`.", "",
         "## Could the question be answered?", ""]
    if grac_ok:
        L.append("Yes, for Gracilariopsis: at least one Gracilariopsis pair has an interpretable ODC-half dN/dS.")
    else:
        L += ["**No.** The question needs a dN/dS for the *Gracilariopsis chorda* ODC half, and no usable nucleotide "
              "pair exists for it. Its only relatives with a fusion in p10 are *Gracilaria gracilis* and *Agarophyton "
              "chilense* (same family, Gracilariaceae); against both, synonymous sites are saturated in every region "
              "(dS > 1 or not computable). So Gracilariopsis cannot be compared with the species that have spare ODCs "
              "by dN/dS from public data. A protein-level fallback is at the end, clearly labelled; it does not answer "
              "the question.", "",
              "The positive control is not usable either: the two Ceratopteris loci are saturated against each other "
              "(dS 1.1-1.8), and copies within a locus are splice isoforms of one gene. The only interpretable ratios "
              "are from *Galdieria partita* vs *G. sulphuraria*, whose ODC halves are both degraded (p10), so there is "
              "no intact-ODC-half ratio to compare them with.", ""]
    L += ["**Premises in the task that do not match the project's own results** (p10_species_enzymes.tsv, p17):", "",
          "- *Galdieria partita*'s standalone ODC-family proteins (A0A9C7PSP6, A0A9C7Q3X8) are 11/13, not fully intact "
          "(Y331 and D332 changed, by sequence; p10 structure 11 kept / 1 lost / 1 unresolved).",
          "- *Percolomonas cosmopolitus*'s two ODC-family TSA proteins are 12/13 by sequence (Y331 to H or F), with no "
          "structure.",
          "- Only *Rhodosorus marinus* has a 13/13 standalone ODC (A0AAV8UQI0, and the same protein in TSA).", "",
          "## Step 1 - which pairs exist", "",
          "From p10_fusion_search.tsv / p10_all_fusions.tsv. No fusion congener exists in p10 for *Rhodosorus marinus* "
          "or *Percolomonas cosmopolitus*, so neither has any pair. Pairs within each group, protein identity from a "
          "pairwise FAMSA alignment. Ceratopteris is shown with one isoform per locus (A0A8T2T6R3, A0A8T2V5S6); "
          "all isoform combinations are in the TSV.", "",
          "| group | pair | type | identity, whole / AdoMetDC / linker / ODC (%) |", "|---|---|---|---|"]
    for r in filter(shown, rows):
        L.append(f"| {r['group']} | {r['seq1']} ({r['organism1']}) - {r['seq2']} ({r['organism2']}) | {r['pair_type']} | "
                 f"{r['all_pct_identity']} / {r['AdoMetDC_pct_identity']} / {r['linker_pct_identity']} / {r['ODC_pct_identity']} |")
    L += ["", "Ceratopteris: the eight fusion entries are four splice isoforms each of two loci "
          f"({meta['A0A8T2T6R3'].get('locus', '?')} and {meta['A0A8T2V5S6'].get('locus', '?')}). Isoform pairs are not "
          "independent comparisons; there is one real paralog pair.", "",
          "## Step 2 - CDS", ""]
    for a, n in sorted(cds_note.items()):
        L.append(f"- {a} ({meta[a]['organism']}; locus {meta[a].get('locus') or '?'}; strain {meta[a].get('strain') or '-'}): {n}")
    dropped = [a for a, n in cds_note.items() if "identical" not in n]
    L += ["", f"Dropped: {', '.join(dropped) if dropped else 'none - every CDS translates to the exact UniParc protein'}.", "",
          "## Step 3 - dN/dS per region", "",
          "Nei-Gojobori 1986 with Jukes-Cantor correction (Biopython NG86), codons threaded on a FAMSA protein alignment "
          "(MAFFT is not installed; for two sequences FAMSA does the same job), gap/stop/ambiguous columns removed. "
          "AdoMetDC = PF01536 span; linker = between AdoMetDC and the first PLP-decarboxylase Pfam match; ODC = first "
          "PLP Pfam start to the C-terminus (per-protein Pfam spans from p10_fusion_search.tsv, because "
          "fusion_architectures.tsv has no per-protein positions; the ODC region runs to the C-terminus because the "
          "PF02784 match is truncated in several red-algal fusions, e.g. G. partita 434-512). "
          f"**saturated** = dS > 1 or not computable; **too little divergence** = dS < {LOW_DS} (my threshold: too few "
          "synonymous changes for a ratio). Neither is interpretable. No significance test was run (NG86 gives point "
          "estimates only). One isoform per Ceratopteris locus is shown; the rest are in the TSV.", "",
          "| pair | type | region | codons | differing codons | dN | dS | dN/dS | status |", "|---|---|---|---|---|---|---|---|---|"]
    for r in filter(shown, rows):
        if r["pair_type"].startswith("isoforms"):
            continue
        for g in ("AdoMetDC", "linker", "ODC"):
            if f"{g}_dN" not in r:
                continue
            st = ("saturated" if r[f"{g}_saturated"] else "too little divergence" if r[f"{g}_too_little_divergence"]
                  else "**interpretable**")
            w = fmt(r, f"{g}_dNdS")
            L.append(f"| {r['seq1']} - {r['seq2']} | {r['pair_type']} | {g} | {r[f'{g}_codons']} | {r[f'{g}_differing_codons']} | "
                     f"{fmt(r, f'{g}_dN')} | {fmt(r, f'{g}_dS')} | {'**' + w + '**' if r[f'{g}_interpretable'] else w} | {st} |")
    good = [r for r in rows if r.get("ODC_interpretable") is True]
    L += ["", "**Interpretable ODC-half ratios:** " + (
        "; ".join(f"{r['seq1']} - {r['seq2']} ({r['organism1']} - {r['organism2']}): ODC {r['ODC_dNdS']}, "
                  f"AdoMetDC {r['AdoMetDC_dNdS']}, linker {r['linker_dNdS']}" for r in good) or "none") + ".", "",
          "These four are effectively one comparison: the two *G. partita* copies are 99.6 % identical and the two "
          "*G. sulphuraria* proteins (strains 074W and NIES-3638) 99.9 %. Every interpretable ratio is below 1 in all "
          "three regions; no test of whether any differs from 1, or from another region, was run.", ""]
    if fb:
        L += ["## Step 4 - FALLBACK, protein level (does NOT answer the question)", "",
              "Per-region identity of each fusion to the consensus of all p10 AdoMetDC-then-PLP fusions "
              f"(n = {fb[0]['n_family']}; each region aligned with FAMSA; consensus = majority residue in columns with "
              "<= 50 % gaps; the protein itself is in the family). One isoform per Ceratopteris locus.", "",
              "**Why this is weaker:** identity to a consensus mixes how long ago a lineage diverged with how "
              "constrained it is; without synonymous sites as a neutral baseline it cannot separate selection from time.", "",
              "| accession | organism | AdoMetDC | linker | ODC |", "|---|---|---|---|---|"]
        for f in fb:
            L.append(f"| {f['accession']} | *{f['organism']}* | {f['AdoMetDC_pct_identity_to_consensus']} % ({f['AdoMetDC_columns']}) | "
                     f"{f['linker_pct_identity_to_consensus']} % ({f['linker_columns']}) | {f['ODC_pct_identity_to_consensus']} % ({f['ODC_columns']}) |")
        L += ["", "Numbers in brackets = consensus columns compared. Linkers vary in length and sequence, so the "
              "linker consensus is especially weak. *Agarophyton chilense* has only a PF00278 match (603-695) after its "
              "AdoMetDC domain, so its ODC region starts late and covers few columns. These values are also written to results/p18_fallback_identity.tsv.", ""]
        with open(RESULTS / "p18_fallback_identity.tsv", "w", encoding="utf8", newline="") as fh:
            w = csv.DictWriter(fh, list(fb[0]), delimiter="\t")
            w.writeheader()
            w.writerows(fb)
    (RESULTS / "p18_report.md").write_text("\n".join(L), encoding="utf8")


if __name__ == "__main__":
    main()
