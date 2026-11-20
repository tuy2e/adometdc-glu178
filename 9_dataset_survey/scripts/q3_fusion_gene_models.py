"""Two gene-model checks on the 23 AdoMetDC-ODC fusions (Track B), prompted by the sandwich protein (q1) turning out
to be a probable gene merge. Fixed before any result was read.

A. Junction intron size. For every genome gene model with a predicted intron in the linker (between the last
   AdoMetDC Pfam match and the first PLP-decarboxylase Pfam match, p7 positions), that intron's length is compared
   with (i) the other introns of the same gene and (ii) a random sample of introns from the same genome
   annotation: 400 CDS records sharing the gene's locus_tag prefix, fetched from ENA (seed 1).
   A merge of two neighbouring genes makes the software predict an "intron" across the intergenic DNA, so it tends
   to be much longer than normal introns (q1: 1,125 and 1,456 nt against a typical 39-62 nt).
   Flags: outlier = larger than every other intron of the same gene AND above the genome sample's 95th
   percentile; strong outlier = the same with the 99th percentile. Reported separately so either half can be read.
   A normal-sized junction intron is NOT evidence that the fusion is real: a merge can also make exons out of
   intergenic or UTR DNA, so the intron itself need not span the gap. Only a large one is a warning sign.

B. Other assemblies and annotations of the same species (or strain, where the taxon is a strain or "sp.").
   1. ENA assemblies under the species taxon (tax_tree): accession, name, version, level, strain, study.
   2. UniProt proteomes under the species taxon (each is one annotation of one assembly).
   3. If there are 2 or more proteomes, every one is downloaded from UniParc (which keeps deleted entries) and
      scanned with the fusion's own two halves (AdoMetDC part = residues 1..AdoMetDC end; PLP part = PLP start..end):
      3-mer prefilter, then local alignment (BLOSUM62, gap 11/1). A hit needs >= 50% coverage of that half and
      >= 30% identity (AdoMetDC) or >= 25% (PLP part; the ODC family is more divergent). Per proteome:
      "fusion" if one protein hits both halves in AdoMetDC-then-PLP order; "separate" if both halves are hit only
      by different proteins; "AdoMetDC only"/"PLP only"/"neither" otherwise. The proteome holding the identical
      fusion sequence is the positive control.
   Assemblies without a UniProt proteome have no protein set to scan; they are listed but not judged.

Reads results/p7_gene_models.tsv, data/ena/*.embl (p7's cache, read only), data/taxonomy.json and the FOLDTRACE
model sequences. Caches in data/q_ena_genome/, data/q_proteomes/ (gitignored: large) and data/q_tax/.
Writes results_81/q3_junction_introns.tsv, q3_genome_intron_samples.tsv, q3_assemblies.tsv,
q3_proteome_scan.tsv, q3_summary.txt.

Run (10-30 minutes, mostly downloads; resumable):  py -3.12 q3_fusion_gene_models.py
"""
import csv
import gzip
import json
import os
import random
import re
import statistics as st
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
RESULTS, DATA, OUT = BASE / "results", BASE / "data", BASE / "results_81"
ENA_P7 = DATA / "ena"
GEN_DIR, PROT_DIR, TAX_DIR = DATA / "q_ena_genome", DATA / "q_proteomes", DATA / "q_tax"
N_SAMPLE, BATCH, SEED = 400, 50, 1
COV, ID_AD, ID_PLP, KMER_MIN = 0.5, 30.0, 25.0, 8
UA = "adometdc-dataset-survey/1.0 (Python urllib; academic research)"


# ------------------------------------------------------------------ web
def fetch(url, tries=6):
    """GET as text. 404/400 return ''. ENA refuses Accept: application/json, so Accept is left open."""
    for k in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read().decode("utf8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404, 410):
                return ""
            wait = int(e.headers.get("Retry-After") or 0) or 2 ** k
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError):
            wait = 2 ** k
        if k < tries:
            time.sleep(wait)
    raise RuntimeError(f"gave up: {url}")


def cached(path, url, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        body = fetch(url)
        tmp = path.with_suffix(path.suffix + ".tmp")
        if gz:
            with gzip.open(tmp, "wt", encoding="utf8") as fh:
                fh.write(body)
        else:
            tmp.write_text(body, encoding="utf8")
        os.replace(tmp, path)
    if gz:
        with gzip.open(path, "rt", encoding="utf8") as fh:
            return fh.read()
    return path.read_text(encoding="utf8")


def gitignore(d, note):
    d.mkdir(parents=True, exist_ok=True)
    if not (d / ".gitignore").exists():
        (d / ".gitignore").write_text(f"# {note}\n*\n!.gitignore\n", encoding="utf8")


# ------------------------------------------------------------------ EMBL parsing
def cds_location(rec):
    ft, grab = [], False
    for l in rec.splitlines():
        if l.startswith("FT   CDS"):
            grab = True
            ft.append(l[21:].strip())
        elif grab and l.startswith("FT                   ") and not l[21:].startswith("/"):
            ft.append(l[21:].strip())
        elif grab:
            break
    return "".join(ft)


def introns(loc):
    """Introns in transcript order: [(length_nt, protein residue after which it falls)]. None if not usable."""
    if not loc or loc.count(":") and len(set(re.findall(r"([A-Za-z0-9_.]+):[<>]?\d", loc))) > 1:
        return None
    spans = [(int(a), int(b)) for a, b in re.findall(r"[<>]?(\d+)\.\.[<>]?(\d+)", loc)]
    if not spans:
        return None
    gaps = [spans[i + 1][0] - spans[i][1] - 1 for i in range(len(spans) - 1)]
    lens = [b - a + 1 for a, b in spans]
    if loc.startswith("complement"):
        lens, gaps = lens[::-1], gaps[::-1]
    out, c = [], 0
    for L, g in zip(lens[:-1], gaps):
        c += L
        out.append((g, round(c / 3, 1)))
    return out


def qual(rec, name):
    m = re.search(rf'/{name}="([^"]+)"', rec)
    return m.group(1) if m else ""


def pct_rank(x, xs):
    return round(100 * sum(v < x for v in xs) / len(xs), 1) if xs else None


def quantile(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else None


# ------------------------------------------------------------------ A. genome intron sample
def genome_sample(prefix, exclude):
    gitignore(GEN_DIR, "ENA CDS records sampled by q3_fusion_gene_models.py (re-downloadable)")
    p = GEN_DIR / f"{prefix}_sample.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf8"))
    q = urllib.parse.quote(f'locus_tag="{prefix}_*"')
    lst = fetch(f"https://www.ebi.ac.uk/ena/portal/api/search?result=coding&query={q}&fields=accession&limit=0&format=tsv")
    accs = sorted({l.split("\t")[0] for l in lst.strip().splitlines()[1:] if l.strip()} - {exclude})
    random.Random(SEED).shuffle(accs)
    pick = accs[:N_SAMPLE]
    lengths, n_genes, n_multi = [], 0, 0
    for i in range(0, len(pick), BATCH):
        body = fetch("https://www.ebi.ac.uk/ena/browser/api/embl/" + ",".join(pick[i:i + BATCH]))
        for rec in body.split("\n//"):
            if "FT   CDS" not in rec:
                continue
            ints = introns(cds_location(rec))
            if ints is None:
                continue
            n_genes += 1
            n_multi += bool(ints)
            lengths += [g for g, _ in ints]
        time.sleep(0.3)
    d = {"prefix": prefix, "n_cds_in_annotation": len(accs) + 1, "n_sampled": len(pick), "n_genes_parsed": n_genes,
         "n_multi_exon": n_multi, "intron_lengths": lengths}
    p.write_text(json.dumps(d), encoding="utf8")
    return d


# ------------------------------------------------------------------ B. assemblies, proteomes, scan
def species_taxon(name, taxid):
    """(species-level name, taxid). Strains and 'sp.' taxa are used as they are."""
    if not name or " sp." in name:
        return name, taxid
    t = cached(TAX_DIR / f"{name.replace(' ', '_')}.json",
               "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/" + urllib.parse.quote(name))
    try:
        j = json.loads(t)
        return name, str(j[0]["taxId"])
    except (ValueError, IndexError, KeyError):
        return name, taxid


def assemblies(tx):
    q = urllib.parse.quote(f"tax_tree({tx})")
    t = cached(TAX_DIR / f"assemblies_{tx}.tsv",
               f"https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query={q}&fields=accession,assembly_name,"
               f"version,assembly_level,strain,study_accession,wgs_set,last_updated,scientific_name&limit=0&format=tsv")
    lines = t.strip().splitlines()
    if len(lines) < 2:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def proteomes(tx):
    t = cached(TAX_DIR / f"proteomes_{tx}.tsv",
               f"https://rest.uniprot.org/proteomes/search?query=taxonomy_id:{tx}"
               f"&fields=upid,organism,organism_id,protein_count,genome_assembly,cpd,busco&format=tsv&size=100")
    lines = t.strip().splitlines()
    if len(lines) < 2:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def proteome_fasta(upid):
    gitignore(PROT_DIR, "UniParc proteome FASTA files downloaded by q3_fusion_gene_models.py (large)")
    t = cached(PROT_DIR / f"{upid}.fasta.gz",
               f"https://rest.uniprot.org/uniparc/stream?query=proteome:{upid}&format=fasta", gz=True)
    seqs, name, buf = {}, None, []
    for l in t.splitlines():
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def kmers(s, k=3):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def scan(seqs, ad, plp):
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    out = {}
    sk = {name: kmers(s) for name, s in seqs.items()}
    for label, q, idmin in (("AdoMetDC", ad, ID_AD), ("PLP", plp, ID_PLP)):
        qk = kmers(q)
        for name, s in seqs.items():
            if len(qk & sk[name]) < KMER_MIN:
                continue
            a = al.align(q, s)[0]
            cols = [(q[i], s[j]) for (q1, q2), (s1, s2) in zip(*a.aligned) for i, j in zip(range(q1, q2), range(s1, s2))]
            cov = (a.aligned[0][-1][1] - a.aligned[0][0][0]) / len(q)
            ident = 100 * sum(x == y for x, y in cols) / len(cols) if cols else 0
            if cov >= COV and ident >= idmin:
                out.setdefault(name, {"len": len(s), "seq": s})[label] = (round(ident, 1), round(cov, 2),
                                                                            a.aligned[1][0][0] + 1, a.aligned[1][-1][1])
    return out


# ------------------------------------------------------------------ main
def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from p2_fetch_annotations import read_candidates          # model sequences from the FOLDTRACE run
    OUT.mkdir(parents=True, exist_ok=True)
    seqs = read_candidates()
    tax = json.loads((DATA / "taxonomy.json").read_text(encoding="utf8"))
    with open(RESULTS / "p7_gene_models.tsv", encoding="utf8") as fh:
        gm = list(csv.DictReader(fh, delimiter="\t"))
    L = [f"q3_fusion_gene_models: {len(gm)} fusions (p7_gene_models.tsv)"]

    # ---------------- A
    jrows, grows, done = [], [], {}
    for r in gm:
        if r["evidence_type"] != "genome":
            continue
        pid = next(x.split(":")[0] for x in r["records"].split("; ") if x.endswith(":genome"))
        rec = (ENA_P7 / f"{pid}.embl").read_text(encoding="utf8")
        ints = introns(cds_location(rec)) or []
        ae, ps = int(r["adometdc_end"]), int(r["plp_start"])
        link = [(g, p) for g, p in ints if ae < p < ps]
        if not link:
            continue
        other = [g for g, p in ints if not ae < p < ps]
        prefix = qual(rec, "locus_tag").split("_")[0]
        if prefix not in done:
            try:
                done[prefix] = genome_sample(prefix, pid.split(".")[0])
            except Exception as ex:                    # noqa: BLE001
                done[prefix] = {"prefix": prefix, "error": f"{type(ex).__name__}: {ex}"[:150], "intron_lengths": []}
        gs = done[prefix]["intron_lengths"]
        for g, p in link:
            largest = all(g > o for o in other) and all(g >= x for x, _ in link)
            p95, p99 = quantile(gs, 0.95), quantile(gs, 0.99)
            flag = ("strong outlier" if largest and p99 is not None and g > p99 else
                    "outlier" if largest and p95 is not None and g > p95 else "not an outlier")
            jrows.append({"accession": r["accession"], "organism": r["organism"], "protein_id": pid,
                          "locus_tag": qual(rec, "locus_tag"), "gene_model_source": qual(rec, "note"),
                          "n_exons": len(ints) + 1, "adometdc_end": ae, "plp_start": ps,
                          "junction_intron_nt": g, "after_residue": p,
                          "other_introns_same_gene_nt": ",".join(map(str, other)),
                          "gene_median_other": st.median(other) if other else "",
                          "gene_max_other": max(other) if other else "",
                          "largest_in_gene": largest,
                          "genome_sample_introns": len(gs), "genome_median": st.median(gs) if gs else "",
                          "genome_p95": p95 if gs else "", "genome_p99": p99 if gs else "",
                          "genome_percentile_of_junction": pct_rank(g, gs) if gs else "",
                          "flag": flag if gs else ("largest in gene; no genome sample" if largest else "not an outlier (gene only)")})
    for prefix, d in done.items():
        gs = d["intron_lengths"]
        grows.append({"locus_tag_prefix": prefix, "cds_in_annotation": d.get("n_cds_in_annotation", ""),
                      "sampled": d.get("n_sampled", ""), "genes_parsed": d.get("n_genes_parsed", ""),
                      "multi_exon_genes": d.get("n_multi_exon", ""), "introns": len(gs),
                      "median": st.median(gs) if gs else "", "p05": quantile(gs, 0.05), "p95": quantile(gs, 0.95),
                      "p99": quantile(gs, 0.99), "max": max(gs) if gs else "", "error": d.get("error", "")})
    write("q3_junction_introns.tsv", jrows)
    write("q3_genome_intron_samples.tsv", grows)
    L.append(f"\nA. junction introns: {len(jrows)} linker introns in {len({x['accession'] for x in jrows})} genes; "
             f"flags {dict(Counter(x['flag'] for x in jrows))}")
    for x in jrows:
        L.append(f"  {x['accession']} {x['organism']}: {x['junction_intron_nt']} nt after residue {x['after_residue']}; "
                 f"same gene median {x['gene_median_other']}, max {x['gene_max_other']}; genome median "
                 f"{x['genome_median']}, p95 {x['genome_p95']}, percentile {x['genome_percentile_of_junction']} -> {x['flag']}")

    # ---------------- B
    arows, srows = [], []
    seen = set()
    for r in gm:
        t = tax.get(r["taxid"], {}).get("ranks", {})
        name, tx = t.get("species") or r["organism"], r["taxid"]
        fus = seqs.get(r["accession"], "")
        ae, ps = int(r["adometdc_end"]), int(r["plp_start"])
        try:
            name, tx = species_taxon(name, tx)
            key = tx
            asm = assemblies(tx) if key not in seen else []
            prs = proteomes(tx)
        except Exception as ex:                        # noqa: BLE001
            arows.append({"accession": r["accession"], "organism": r["organism"], "species_taxid": tx,
                          "error": f"{type(ex).__name__}: {ex}"[:150]})
            continue
        if key not in seen:
            for a in asm:
                arows.append({"fusion_accession": r["accession"], "species": name, "species_taxid": tx, **a})
            if not asm:
                arows.append({"fusion_accession": r["accession"], "species": name, "species_taxid": tx,
                              "accession": "(no ENA assembly)"})
        seen.add(key)
        if len(prs) >= 2 and fus:
            for pr in prs:
                upid = pr.get("Proteome Id", "")
                try:
                    ps_ = proteome_fasta(upid)
                    hits = scan(ps_, fus[:ae], fus[ps - 1:])
                except Exception as ex:                # noqa: BLE001
                    srows.append({"fusion_accession": r["accession"], "proteome": upid, "error": f"{type(ex).__name__}"})
                    continue
                both = [n for n, h in hits.items() if "AdoMetDC" in h and "PLP" in h and h["AdoMetDC"][3] < h["PLP"][2] + 50]
                ad_only = [n for n, h in hits.items() if "AdoMetDC" in h and n not in both]
                plp_only = [n for n, h in hits.items() if "PLP" in h and n not in both]
                own = any(h["seq"] == fus for h in hits.values())
                verdict = ("fusion" if both else "separate" if ad_only and plp_only else
                           "AdoMetDC only" if ad_only else "PLP only" if plp_only else "neither")
                srows.append({"fusion_accession": r["accession"], "organism": r["organism"], "proteome": upid,
                              "proteome_organism": pr.get("Organism", ""), "assembly": pr.get("Genome assembly ID", ""),
                              "protein_count": pr.get("Protein count", ""), "proteins_scanned": len(ps_),
                              "contains_identical_fusion": own, "verdict": verdict,
                              "fusion_hits": "; ".join(f"{n} ({hits[n]['len']} aa)" for n in both),
                              "adometdc_only_hits": "; ".join(f"{n} ({hits[n]['len']} aa, {hits[n]['AdoMetDC'][0]}% id)" for n in ad_only),
                              "plp_only_hits": "; ".join(f"{n} ({hits[n]['len']} aa, {hits[n]['PLP'][0]}% id)" for n in plp_only)})
    write("q3_assemblies.tsv", arows)
    write("q3_proteome_scan.tsv", srows)
    per = Counter(x["fusion_accession"] for x in arows if x.get("accession") and not x["accession"].startswith("("))
    L.append(f"\nB. ENA assemblies per fusion's species: " + ", ".join(f"{a} {n}" for a, n in sorted(per.items())))
    multi = [a for a, n in per.items() if n >= 2]
    L.append(f"  species with >= 2 ENA assemblies: {len(multi)}")
    for x in arows:
        if x.get("fusion_accession") in multi:
            L.append(f"    {x['fusion_accession']} {x.get('species', '')}: {x.get('accession', '')} {x.get('assembly_name', '')} "
                     f"v{x.get('version', '')} {x.get('assembly_level', '')} strain={x.get('strain', '')} "
                     f"study={x.get('study_accession', '')} updated={x.get('last_updated', '')}")
    L.append(f"  proteome scans: {len(srows)}")
    for x in srows:
        L.append(f"    {x['fusion_accession']} {x.get('proteome', '')} {x.get('assembly', '')} own={x.get('contains_identical_fusion', '')} "
                 f"-> {x.get('verdict', x.get('error', ''))} | fusion: {x.get('fusion_hits', '')} | AdoMetDC only: "
                 f"{x.get('adometdc_only_hits', '')} | PLP only: {x.get('plp_only_hits', '')}")
    txt = "\n".join(L) + "\n"
    (OUT / "q3_summary.txt").write_text(txt, encoding="utf8")
    print(txt)


def write(name, rows):
    if not rows:
        (OUT / name).write_text("(no rows)\n", encoding="utf8")
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT / name, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
