"""Track B step 3: the red-algal pattern in the ODC half of the AdoMetDC-ODC fusions.

Pattern under test: the fusions with confident substitutions at the 13 human-ODC functional residues
(p9 gated: residue in place, pLDDT >= 70, not the human residue) are Galdieria, Rhodosorus, Gracilariopsis
(Rhodophyta) and Percolomonas (Heterolobosea); no green-lineage fusion has one.

1. GROW THE SAMPLE. UniParc searched directly (no structure prefilter) for sequences carrying a Pfam
   AdoMetDC (PF01536/PF02675) AND a PLP-decarboxylase (PF02784/PF00278) feature; UniProtKB queried the same
   way as a cross-check. Everything not in the 23 is listed with taxonomy. New ones with an AlphaFold DB
   model get the p9 gated AdoMetDC sites (Marc's five) and the p7 partner integrity + p9 gated ODC sites.
   Models go to data/p10_models/ (never Marc's folder).
2. ONE EVENT OR MANY. ODC-half region (first PLP Pfam start to last PLP Pfam end) of every fusion with a
   region >= 250 aa, plus standalone references, aligned with FAMSA (pyfamsa), columns > 50% gaps removed,
   Kimura protein distance, neighbour joining, 1,000 column bootstraps (seed 1). Support is reported for the
   exact splits of interest (red-algal fusions; red + Percolomonas; the degraded set; green fusions).
   This is a distance tree, not maximum likelihood. Separately, the 13 sites of every degraded protein,
   residue by residue, and for each pair whether their substitutions are the same residue.
3. RED-ALGAL DECARBOXYLASES. For Galdieria partita, Rhodosorus marinus, Gracilariopsis chorda and
   Percolomonas cosmopolitus: every UniParc sequence of the taxon with a PLP-decarboxylase, pyruvoyl ADC
   (PF01862), agmatine iminohydrolase (PF04371), arginase/agmatinase (PF00491), nitrilase/CN-hydrolase
   (PF00795, the NCPAH family) or AdoMetDC feature. PLP decarboxylases are classified by the p6 method
   (local alignment, BLOSUM62, gap 11/1, nearest reference of known activity) against p6's references plus
   LDC (O50657 Selenomonas LDC/ODC, A0A8F1NN86 Phlegmariurus LDC2) and Arabidopsis DAPDC1 (Q949X7); the 13 ODC
   residues are read from the alignment to human ODC, and, where an AlphaFold model exists, structurally with
   the pLDDT gate. UniProt proteome size and BUSCO are recorded, because absence from an incomplete proteome
   is not absence from the genome.
4. RETAINED OR ROTTING. Partner length, TM to human ODC, ODC positions aligned, longest unaligned stretch,
   mean pLDDT and % of partner residues below 70, by ODC-half status.
5. CORRELATES across the 23 (descriptive only, no tests): AdoMetDC half (p9), linker length, protein length,
   putrescine-site residues (p8), gene-model source (p7).

ODC-half status (used throughout):
    degraded     >= 1 gated loss (confident in-place substitution) among the 13
    incomplete   no gated loss; partner < 300 aa or >= 100 human-ODC positions unaligned in one stretch
    intact       13/13 kept
    undetermined otherwise (some sites unresolved)

Needs pyfamsa for part 2:  py -3.12 -m pip install pyfamsa
Run:  py -3.12 p10_red_algal_pattern.py   (network: UniParc, UniProt, AlphaFold DB; ~60 models; then the
      bootstrap; caches in data/p10/ and data/p10_models/, so re-runs are fast)
Writes results/p10_*.tsv, p10_odc_tree.nwk, p10_odc_alignment.fasta, p10_report.md
"""
import csv
import json
import math
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from foldtrace.io import load_structure

import p6_decarboxylase_partner as p6
import p7_fusion_validation as p7
import p8_putrescine_site as p8
from p2_fetch_annotations import DATA, RESULTS, get, read_candidates
from p3_residue_grid import REF_PDB, REF_CHAIN, PLDDT_OK, RUN, human_sequence, plddt_by_resnum
from p9_self_processing import ODC_SITE_LABELS, marc_sites

sys.stdout.reconfigure(encoding="utf8")
CACHE, MODELS = DATA / "p10", DATA / "p10_models"
CACHE.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)
ADOMETDC, PLP = {"PF01536", "PF02675"}, {"PF02784", "PF00278"}
ROUTE = {"PF01862": "pyruvoyl ADC (PdaD family)", "PF04371": "agmatine iminohydrolase (AIH) family",
         "PF00491": "arginase/agmatinase family", "PF00795": "nitrilase/CN-hydrolase family (NCPAH candidate)"}
SPECIES = {83374: "Galdieria partita", 101924: "Rhodosorus marinus", 448386: "Gracilariopsis chorda",
           63605: "Percolomonas cosmopolitus"}
EXTRA_REFS = {"O50657": ("LDC/ODC Selenomonas ruminantium", "reviewed, EC 4.1.1.18 + 4.1.1.17"),
              "A0A8F1NN86": ("LDC2 Phlegmariurus tetrastichus (plant)", "reviewed, EC 4.1.1.18"),
              "Q949X7": ("DAPDC1 Arabidopsis", "reviewed, EC 4.1.1.20")}
NCPAH_REF = "Q8VYF5"   # Arabidopsis NLP1, N-carbamoylputrescine amidase, reviewed EC 3.5.1.53
FOCUS = ["Rhodophyta", "Glaucocystophyceae", "Cryptophyceae", "Heterolobosea", "Haptophyta", "Stramenopiles",
         "Alveolata", "Rhizaria", "Viridiplantae", "Metazoa", "Fungi", "Bacteria", "Archaea"]
ODC_FAMILY = ("ODC human", "ODC Trypanosoma", "ODC yeast", "ODC tobacco", "AdoMetDC/ODC Plasmodium")
N_BOOT, MIN_REGION = 1000, 250
SMOKE = "--smoke" in sys.argv   # check run: new red-algal fusions only, 20 bootstraps, tree skipped without pyfamsa
if SMOKE:
    N_BOOT = 20
PUTRESCINE = [13, 15, 111, 174, 176, 178, 256, 285]   # p8 positions


# ---------------------------------------------------------------- fetching (cached)
def cached(name, fn):
    p = CACHE / name
    if p.exists():
        return json.loads(p.read_text(encoding="utf8"))
    obj = fn()
    p.write_text(json.dumps(obj), encoding="utf8")
    return obj


def uniparc_search(query):
    def run():
        out, url = [], ("https://rest.uniprot.org/uniparc/search?" + urllib.parse.urlencode(
            {"query": query, "format": "json", "size": 500}))
        while url:
            st, hd, body = get(url)
            out += json.loads(body)["results"]
            m = re.search(r'<([^>]+)>;\s*rel="next"', hd.get("Link", "") or hd.get("link", ""))
            url = m.group(1) if m else None
        return out
    return cached("uniparc_" + re.sub(r"\W+", "_", query)[:120] + ".json", run)


def uniprot_tsv(query, fields):
    def run():
        st, hd, body = get("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": query, "format": "tsv", "fields": fields, "size": 500}))
        return body
    body = cached("uniprot_" + re.sub(r"\W+", "_", query)[:120] + ".json", run)
    return list(csv.DictReader(body.splitlines(), delimiter="\t"))


def lineage(taxid):
    def run():
        st, _, body = get(f"https://rest.uniprot.org/taxonomy/{taxid}?format=json")
        d = json.loads(body) if st == 200 else {}
        return {"name": d.get("scientificName", ""), "lineage": [x.get("scientificName", "") for x in d.get("lineage", [])]}
    return cached(f"tax_{taxid}.json", run)


def group_of(lin):
    names = set(lin["lineage"]) | {lin["name"]}
    return next((g for g in FOCUS if g in names), "other")


def proteome(taxid):
    def run():
        st, _, body = get("https://rest.uniprot.org/proteomes/search?" + urllib.parse.urlencode(
            {"query": f"organism_id:{taxid}", "format": "json"}))
        out = []
        for r in json.loads(body).get("results", []):
            b = ((r.get("proteomeCompletenessReport") or {}).get("buscoReport") or {})
            out.append({"id": r.get("id"), "type": r.get("proteomeType"), "proteins": r.get("proteinCount"),
                        "busco": f"C {b.get('complete')} / F {b.get('fragmented')} / M {b.get('missing')} of {b.get('total')} ({b.get('lineageDb')})" if b else ""})
        return out
    return cached(f"proteome_{taxid}.json", run)


def af_model(accs, seq):
    """First accession with an AlphaFold DB model whose sequence matches; None if none."""
    for a in dict.fromkeys(x.split(".")[0] for x in accs):   # UniParc lists inactive accessions as ACC.version
        p = MODELS / f"AF-{a}.pdb"
        miss = MODELS / f"AF-{a}.none"
        if miss.exists():
            continue
        if not p.exists():
            try:
                with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{a}", timeout=60) as r:
                    st, body = r.status, r.read().decode("utf8")
            except urllib.error.HTTPError as e:   # 400/404: no model for this accession
                st = e.code
            if st != 200:
                miss.write_text(f"HTTP {st}")
                continue
            ent = json.loads(body)[0]
            mseq = ent.get("uniprotSequence") or ent.get("sequence")
            if mseq and mseq != seq:
                miss.write_text("model sequence differs from this UniParc sequence")
                continue
            st, _, pdb = get(ent["pdbUrl"])
            p.write_text(pdb, encoding="utf8")
        return a, p
    return None, None


def pfam_spans(entry):
    out = []
    for f in entry.get("sequenceFeatures", []):
        if f.get("database") == "Pfam":
            for loc in f.get("locations", []):
                out.append((f["databaseId"], loc["start"], loc["end"]))
    return sorted(out, key=lambda x: x[1])


# ---------------------------------------------------------------- p6 nearest reference, extended
LOCAL = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                        open_gap_score=-11, extend_gap_score=-1)


def references():
    refs = p6.references()
    missing = [a for a in list(EXTRA_REFS) + [NCPAH_REF] if a not in refs]
    if missing:
        def run():
            st, _, body = get("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(missing)
                              + "&format=fasta")
            return body
        body = cached("extra_refs.json", run)
        acc = None
        for ln in body.splitlines():
            if ln.startswith(">"):
                acc = ln.split("|")[1]
                refs[acc] = ""
            elif acc:
                refs[acc] += ln.strip()
    return refs


def nearest(q, refs):
    labels = {**p6.REFS, **EXTRA_REFS}
    hits = []
    for acc, (label, _) in labels.items():
        a = LOCAL.align(q, refs[acc])[0]
        pairs = [(i, j) for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))]
        ident = sum(q[i] == refs[acc][j] for i, j in pairs)
        hits.append({"acc": acc, "label": label, "score": a.score, "ident": ident / len(pairs) if pairs else 0,
                     "cov": len(pairs) / len(q), "pairs": pairs})
    hits.sort(key=lambda h: -h["score"])
    best = hits[0]
    fam = lambda h: h["label"].split()[0].split("/")[-1]
    other = next(h for h in hits[1:] if fam(h) != fam(best))
    h = next(x for x in hits if x["acc"] == "P11926")
    odc_res = {j + 1: q[i] for i, j in h["pairs"]}
    return best, other, odc_res


def site_string(odc_res):
    return " ".join(f"{lab}{'=' if odc_res.get(int(lab[1:])) == lab[0] else ':' + odc_res.get(int(lab[1:]), '-')}"
                    for lab in ODC_SITE_LABELS)


# ---------------------------------------------------------------- structural reads with the pLDDT gate
def gate(call_aa, human_aa, off, p):
    if call_aa in (None, "-"):
        return "unresolved", "no_residue_aligned"
    if off > p8.OFFSET:
        return "unresolved", "displaced" + (";plddt_below_70" if p is not None and p < PLDDT_OK else "")
    if p is None or p < PLDDT_OK:
        return "unresolved", "plddt_below_70"
    return ("kept" if call_aa == human_aa else "lost"), ""


def odc_structural(odc, path, lo, hi=10 ** 6):
    tm, n_al, sites, maxgap = p7.structural(odc, path, lo, hi)
    pl = plddt_by_resnum(path)
    out = {}
    for lab in ODC_SITE_LABELS:
        aa, rn, off = sites[int(lab[1:])]
        p = pl.get(rn) if rn != "" else None
        st, why = gate(aa if aa != "-" else None, lab[0], off, p)
        out[lab] = {"aa": aa, "resnum": rn, "offset": off, "plddt": p, "state": st, "why": why}
    part = [v for n, v in pl.items() if lo <= n <= hi]
    return {"tm": tm, "n_aligned": n_al, "maxgap": maxgap, "sites": out,
            "partner_len": len(part), "partner_plddt": float(np.mean(part)),
            "partner_low_pct": 100 * sum(v < PLDDT_OK for v in part) / len(part)}


def odc_status(states, partner_len, maxgap):
    c = Counter(states)
    if c["lost"]:
        return "degraded"
    if partner_len < 300 or maxgap >= 100:
        return "incomplete"
    return "intact" if c["kept"] == len(states) else "undetermined"


# ---------------------------------------------------------------- distance tree
def nj(D, names):
    """Neighbour joining. Returns (newick without support, list of leaf-index frozensets, one per internal node)."""
    D = D.copy()
    nodes = [(frozenset([i]), names[i]) for i in range(len(names))]
    clusters = []
    while len(nodes) > 3:
        n = len(nodes)
        r = D.sum(1)
        Q = (n - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        li = max(0.0, 0.5 * D[i, j] + (r[i] - r[j]) / (2 * (n - 2)))
        lj = max(0.0, D[i, j] - li)
        s = nodes[i][0] | nodes[j][0]
        clusters.append(s)
        new = (s, f"({nodes[i][1]}:{li:.4f},{nodes[j][1]}:{lj:.4f})#{len(clusters) - 1}")
        dk = 0.5 * (D[i] + D[j] - D[i, j])
        keep = [k for k in range(n) if k not in (i, j)]
        D = np.vstack([np.hstack([D[np.ix_(keep, keep)], dk[keep][:, None]]), np.append(dk[keep], 0.0)])
        nodes = [nodes[k] for k in keep] + [new]
    a, b, c = nodes
    l = [max(0.0, x) for x in (0.5 * (D[0, 1] + D[0, 2] - D[1, 2]), 0.5 * (D[0, 1] + D[1, 2] - D[0, 2]),
                                0.5 * (D[0, 2] + D[1, 2] - D[0, 1]))]
    return f"({a[1]}:{l[0]:.4f},{b[1]}:{l[1]:.4f},{c[1]}:{l[2]:.4f});", clusters


def canon(s, n):
    return frozenset(s) if 0 not in s else frozenset(range(n)) - frozenset(s)


def distances(X, G, w):
    """Kimura protein distance from one-hot X (aa x n x L) and non-gap mask G (n x L), column weights w."""
    valid = (G * w) @ G.T
    same = sum((x * w) @ x.T for x in X)
    p = 1 - same / np.maximum(valid, 1)
    arg = 1 - p - 0.2 * p * p
    d = np.where(arg > 0.05, -np.log(np.maximum(arg, 0.05)), -math.log(0.05))
    np.fill_diagonal(d, 0)
    return d


def tree(seqs, labels):
    try:
        from pyfamsa import Aligner, Sequence
    except ImportError:
        if SMOKE:
            return None
        sys.exit("part 2 needs pyfamsa:  py -3.12 -m pip install pyfamsa")
    ids = list(seqs)
    msa = Aligner().align([Sequence(i.encode(), seqs[i].encode()) for i in ids])
    aln = {s.id.decode(): s.sequence.decode() for s in msa}
    A = np.array([list(aln[i]) for i in ids])
    keep = (A == "-").mean(0) <= 0.5
    A = A[:, keep]
    with open(RESULTS / "p10_odc_alignment.fasta", "w", encoding="utf8") as fh:
        for i in ids:
            fh.write(f">{i} {labels[i]}\n{aln[i]}\n")
    aas = sorted(set("".join(A.ravel())) - {"-"})
    X = [(A == a).astype(float) for a in aas]
    G = (A != "-").astype(float)
    n, L = A.shape
    names = [re.sub(r"[^\w.|-]", "_", f"{i}|{labels[i]}") for i in ids]
    nwk, clusters = nj(distances(X, G, np.ones(L)), names)
    splits = [canon(c, n) for c in clusters]
    counts = Counter()
    rng = np.random.default_rng(1)
    for b in range(N_BOOT):
        w = np.bincount(rng.integers(0, L, L), minlength=L).astype(float)
        _, cl = nj(distances(X, G, w), names)
        counts.update({canon(c, n) for c in cl})
        if (b + 1) % 200 == 0:
            print(f"  bootstrap {b + 1}/{N_BOOT}", flush=True)
    for k, s in enumerate(splits):
        nwk = nwk.replace(f")#{k}:", f"){100 * counts[s] / N_BOOT:.0f}:", 1).replace(f")#{k},", f"){100 * counts[s] / N_BOOT:.0f},", 1)
    nwk = re.sub(r"\)#\d+", ")", nwk)
    (RESULTS / "p10_odc_tree.nwk").write_text(nwk + "\n", encoding="utf8")
    return ids, splits, counts, n, L


# ---------------------------------------------------------------- main
def main():
    human = human_sequence()
    ref = load_structure(str(REF_PDB), chain=REF_CHAIN)
    ref_idx = {num: k for k, num in enumerate(ref.resnums)}
    p8.POSITIONS = [s[1] for s in marc_sites()]
    altered = {s[0]: s[4] for s in marc_sites()}
    refs = references()
    odc_s = load_structure(str(p7.af_model("P11926")))
    odc = p7.sub_structure(odc_s, 1, p7.ODC_ORDERED)
    our = {r["accession"]: r for r in p7.read("p6_decarboxylase_partner.tsv")}
    cand = read_candidates()
    our_seq = {cand[a]: a for a in our}
    p5 = {r["accession"]: r for r in p7.read("p5_loss_vs_fusion_proteins.tsv")}
    arch = {r["accession"]: r for r in p7.read("architecture_per_protein.tsv")}
    rep = ["# Track B step 3: the red-algal pattern (p10, generated tables)", "",
           "Generated by `scripts/p10_red_algal_pattern.py`; definitions in its docstring. No p-values anywhere.", ""]

    # ---- 1. search
    q_uniparc = "(feature_id:PF01536 OR feature_id:PF02675) AND (feature_id:PF02784 OR feature_id:PF00278)"
    hits = uniparc_search(q_uniparc)
    q_uniprot = "(xref:pfam-PF01536 OR xref:pfam-PF02675) AND (xref:pfam-PF02784 OR xref:pfam-PF00278)"
    up = uniprot_tsv(q_uniprot, "accession,organism_name,length")
    found = []
    for e in hits:
        accs = e.get("uniProtKBAccessions") or []
        tax = (e.get("commonTaxons") or [{}])[0]
        taxid = tax.get("commonTaxonId")
        lin = lineage(taxid) if taxid else {"name": "", "lineage": []}
        sp = pfam_spans(e)
        ae = max((s[2] for s in sp if s[0] in ADOMETDC), default=None)
        plp = [s for s in sp if s[0] in PLP]
        in23 = next((a for a in accs if a in our), None) or our_seq.get(e["sequence"]["value"])
        found.append({"upi": e["uniParcId"], "uniprot_accessions": ",".join(accs), "in_our_23": in23 or "",
                      "taxid": taxid, "organism": lin["name"] or tax.get("commonTaxon", ""), "group": group_of(lin),
                      "lineage": "; ".join(lin["lineage"][-8:]), "length": e["sequence"]["length"],
                      "pfam_architecture": " > ".join(f"{a} {b}-{c}" for a, b, c in sp),
                      "adometdc_then_plp": bool(ae and any(s[1] > ae for s in plp)),
                      "plp_start": min((s[1] for s in plp if ae and s[1] > ae), default=min((s[1] for s in plp), default="")),
                      "plp_end": max((s[2] for s in plp), default=""), "adometdc_end": ae or "",
                      "_seq": e["sequence"]["value"], "_accs": accs})
    cols1 = [k for k in found[0] if not k.startswith("_")]
    with open(RESULTS / "p10_fusion_search.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, cols1, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(found, key=lambda r: (r["group"], r["organism"])))
    up_accs = {r["Entry"] for r in up}
    covered = {a for r in found for a in r["_accs"]}
    new = [r for r in found if not r["in_our_23"]]
    rep += ["## 1. Grow the sample", "",
            f"UniParc query `{q_uniparc}`: {len(hits)} sequences, {len(found) - len(new)} of them in our 23, "
            f"{len(new)} new. UniProtKB query `{q_uniprot}`: {len(up)} active entries; "
            f"{len(up_accs - covered)} not among the UniParc hits ({', '.join(sorted(up_accs - covered)) or 'none'}).",
            "", "New sequences by group (AdoMetDC before PLP / other order):", ""]
    for g, c in Counter(r["group"] for r in new).most_common():
        o = Counter(r["adometdc_then_plp"] for r in new if r["group"] == g)
        rep.append(f"- {g}: {c} ({o[True]} / {o[False]})")
    rep += ["", "| UPI | UniProt | organism | group | length | AdoMetDC then PLP |", "|---|---|---|---|---|---|"]
    rep += [f"| {r['upi']} | {r['uniprot_accessions'] or '(none)'} | *{r['organism']}* | {r['group']} | {r['length']} | "
            f"{'yes' if r['adometdc_then_plp'] else 'no'} |" for r in sorted(new, key=lambda r: (FOCUS.index(r['group']) if r['group'] in FOCUS else 99, r['organism']))]

    # ---- analyse every fusion (our 23 from p9/p7 tables, new ones here)
    gated = defaultdict(dict)
    for r in p7.read("p9_odc_sites_gated.tsv"):
        gated[r["accession"]][r["odc_site"]] = {"aa": r["cand_aa"], "state": r["state"], "plddt": r["plddt"],
                                                "why": r["unresolved_reason"]}
    integ = {r["accession"]: r for r in p7.read("p7_partner_integrity.tsv")}
    both = {r["accession"]: r for r in p7.read("p9_both_halves.tsv")}
    fus = []
    for a in our:
        m = RUN / "models" / f"{a}.pdb"
        pl = plddt_by_resnum(m)
        lo = int(integ[a]["partner_from"])
        part = [v for n, v in pl.items() if n >= lo]
        doms = [(d.split()[0], *map(int, d.split()[1].split("-"))) for d in arch[a]["domains_with_positions"].split("; ")]
        ae = max(e for pf, s_, e in doms if pf in ADOMETDC)
        plp = [(s_, e) for pf, s_, e in doms if pf in PLP and s_ > ae]
        st = odc_status([gated[a][l]["state"] for l in ODC_SITE_LABELS], int(integ[a]["partner_len"]),
                        int(integ[a]["longest_unaligned_ODC_stretch"]))
        fus.append({"accession": a, "set": "our 23", "organism": our[a]["organism"], "group": group_of(lineage(int(p5[a]["taxid"]))),
                    "model": "Marc's run", "adometdc_half": both[a]["adometdc_half"], "odc_status": st,
                    "odc_k_l_u": f"{both[a]['odc_kept']}/{both[a]['odc_lost']}/{both[a]['odc_unresolved']}",
                    "partner_len": int(integ[a]["partner_len"]), "TM_to_human_ODC": float(integ[a]["TM_to_human_ODC"]),
                    "odc_positions_aligned": integ[a]["n_ODC_positions_aligned"],
                    "longest_unaligned": int(integ[a]["longest_unaligned_ODC_stretch"]),
                    "partner_plddt_mean": round(float(np.mean(part)), 1),
                    "partner_pct_below_70": round(100 * sum(v < PLDDT_OK for v in part) / len(part), 1),
                    "closest_reference": our[a]["closest_reference"], "sites": gated[a], "_seq": cand[a],
                    "_plp": (min(x[0] for x in plp), max(x[1] for x in plp)) if plp else ("", "")})
    for r in new:
        if not r["plp_start"] or (SMOKE and r["group"] != "Rhodophyta"):
            continue
        seq, lo = r["_seq"], int(r["plp_start"])
        q = seq[(r["adometdc_end"] or lo - 1):] if r["adometdc_then_plp"] else seq[lo - 1:int(r["plp_end"])]
        best, other, odc_res = nearest(q, refs)
        acc, path = af_model(r["_accs"], seq)
        row = {"accession": acc or (r["_accs"][0] if r["_accs"] else r["upi"]), "set": "new", "organism": r["organism"],
               "group": r["group"], "model": "AlphaFold DB" if path else "none (sequence only)",
               "closest_reference": f"{best['label']} ({best['ident']:.2f} id; next family {other['label']}, "
                                    f"margin {best['score'] - other['score']:.0f})",
               "sequence_sites": site_string(odc_res), "_seq": seq, "_plp": (r["plp_start"], r["plp_end"]),
               "adometdc_then_plp": r["adometdc_then_plp"]}
        if path:
            sites, _ = p8.map_sites(f"p10_{acc}", ref, ref_idx, human, mpath=path)
            ado_states = []
            for lab, pos, *_ in marc_sites():
                s = sites[pos]
                stt = "altered" if s["state"] == "lost" and s["cand_aa"] in altered[lab] else s["state"]
                ado_states.append(stt)
            c = Counter(ado_states)
            o = odc_structural(odc, path, max(1, lo - 25), 10 ** 6 if r["adometdc_then_plp"] else int(r["plp_end"]) + 25)
            row.update({"adometdc_half": "not intact" if c["lost"] else ("undetermined" if c["unresolved"] else "intact"),
                        "adometdc_sites": " ".join(f"{lab}:{s}" for (lab, *_), s in zip(marc_sites(), ado_states)),
                        "odc_status": odc_status([v["state"] for v in o["sites"].values()], o["partner_len"], o["maxgap"]),
                        "odc_k_l_u": "/".join(str(Counter(v["state"] for v in o["sites"].values())[k])
                                              for k in ("kept", "lost", "unresolved")),
                        "partner_len": o["partner_len"], "TM_to_human_ODC": round(o["tm"], 3),
                        "odc_positions_aligned": f"{o['n_aligned']}/{p7.ODC_ORDERED}", "longest_unaligned": o["maxgap"],
                        "partner_plddt_mean": round(o["partner_plddt"], 1),
                        "partner_pct_below_70": round(o["partner_low_pct"], 1), "sites": o["sites"]})
        fus.append(row)
    fcols = ["accession", "set", "organism", "group", "model", "adometdc_half", "adometdc_sites", "odc_status",
             "odc_k_l_u", "partner_len", "TM_to_human_ODC", "odc_positions_aligned", "longest_unaligned",
             "partner_plddt_mean", "partner_pct_below_70", "closest_reference", "sequence_sites", "adometdc_then_plp"]
    with open(RESULTS / "p10_all_fusions.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fcols + [f"ODC_{l}" for l in ODC_SITE_LABELS], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for f in fus:
            w.writerow({**{k: f.get(k, "") for k in fcols},
                        **{f"ODC_{l}": (f"{f['sites'][l]['aa']} {f['sites'][l]['state']}" if "sites" in f else "")
                           for l in ODC_SITE_LABELS}})
    nf = [f for f in fus if f["set"] == "new"]
    rep += ["", "New fusions analysed (AlphaFold DB model where one exists; same gates as p9):", "",
            "| accession | organism | group | AdoMetDC half | ODC status | ODC kept/lost/unresolved | closest reference |",
            "|---|---|---|---|---|---|---|"]
    rep += [f"| {f['accession']} | *{f['organism']}* | {f['group']} | {f.get('adometdc_half', 'no model')} | "
            f"{f.get('odc_status', 'no model')} | {f.get('odc_k_l_u', '')} | {f['closest_reference']} |" for f in nf]
    rep += ["", "ODC-half status by group, all fusions with a model:", ""]
    for g in sorted({f["group"] for f in fus if "odc_status" in f}):
        rep.append(f"- {g}: " + ", ".join(f"{k} {v}" for k, v in
                                          Counter(f["odc_status"] for f in fus if f["group"] == g and "odc_status" in f).most_common()))

    # ---- 2a. substitutions in the degraded proteins
    deg = [f for f in fus if f.get("odc_status") == "degraded"]
    rep += ["", "## 2. One event or many", "", "### Substitutions at the 13 ODC residues, degraded proteins", "",
            "Cell: residue in the fusion; upper case = confident (pLDDT >= 70, in place), lower case = unresolved, "
            "'.' = human residue kept, '-' = nothing aligned.", "",
            "| site | " + " | ".join(f"{f['organism']} {f['accession']}" for f in deg) + " | same substitute in >= 2 |",
            "|---|" + "---|" * len(deg) + "---|"]
    sub_rows = []
    for lab in ODC_SITE_LABELS:
        cells, subs = [], Counter()
        for f in deg:
            s = f["sites"][lab]
            aa = s["aa"] if s["aa"] not in ("", None) else "-"
            if s["state"] == "kept":
                cells.append(".")
            elif s["state"] == "lost":
                cells.append(aa)
                subs[aa] += 1
            else:
                cells.append(aa.lower() if aa != "-" else "-")
        shared = ", ".join(f"{a} x{c}" for a, c in subs.items() if c >= 2) or "none"
        rep.append(f"| {lab} | " + " | ".join(cells) + f" | {shared} |")
        sub_rows.append({"site": lab, **{f["accession"]: c for f, c in zip(deg, cells)}, "shared_confident": shared})
    pairs = []
    for i in range(len(deg)):
        for j in range(i + 1, len(deg)):
            a, b = deg[i], deg[j]
            both_lost = [l for l in ODC_SITE_LABELS if a["sites"][l]["state"] == b["sites"][l]["state"] == "lost"]
            same = [l for l in both_lost if a["sites"][l]["aa"] == b["sites"][l]["aa"]]
            pairs.append({"a": f"{a['organism']} {a['accession']}", "b": f"{b['organism']} {b['accession']}",
                          "sites_both_lost": len(both_lost), "same_substitute": len(same),
                          "same_at": " ".join(f"{l}{a['sites'][l]['aa']}" for l in same),
                          "different_at": " ".join(f"{l}:{a['sites'][l]['aa']}/{b['sites'][l]['aa']}"
                                                   for l in both_lost if l not in same)})
    with open(RESULTS / "p10_degraded_substitutions.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(sub_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sub_rows)
    with open(RESULTS / "p10_degraded_pairs.tsv", "w", newline="", encoding="utf8") as fh:
        if pairs:
            w = csv.DictWriter(fh, list(pairs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(pairs)
    rep += ["", "Pairs (sites where both are confidently substituted; how many with the same residue):", "",
            "| A | B | both lost | same residue | same at | different at |", "|---|---|---|---|---|---|"]
    rep += [f"| {p['a']} | {p['b']} | {p['sites_both_lost']} | {p['same_substitute']} | {p['same_at']} | {p['different_at']} |"
            for p in pairs]

    # ---- 3. the four species
    q_fam = " OR ".join(f"feature_id:{x}" for x in sorted(PLP | ADOMETDC | set(ROUTE)))
    sp_rows, species_plp = [], {}
    rep += ["", "## 3. Basic amino acid decarboxylases and arginine-route enzymes in the four species", ""]
    for taxid, name in SPECIES.items():
        ents = uniparc_search(f"organism_id:{taxid} AND ({q_fam})")
        prots = proteome(taxid)
        rep.append(f"### *{name}* (taxid {taxid})")
        rep.append("UniProt proteome: " + ("; ".join(f"{p['id']} {p['type']}, {p['proteins']} proteins, BUSCO {p['busco'] or 'n/a'}"
                                                      for p in prots) or "none") + f". UniParc sequences found: {len(ents)}.")
        rep += ["", "| UPI | UniProt | length | Pfam | class | nearest reference (identity; margin to next family) | "
                "13 ODC residues by sequence ('=' human) | structure, gated k/l/u |", "|---|---|---|---|---|---|---|---|"]
        for e in ents:
            sp = pfam_spans(e)
            ids = {s[0] for s in sp}
            seq = e["sequence"]["value"]
            accs = e.get("uniProtKBAccessions") or []
            row = {"taxid": taxid, "species": name, "upi": e["uniParcId"], "uniprot": ",".join(accs),
                   "length": len(seq), "pfam": " > ".join(f"{a} {b}-{c}" for a, b, c in sp)}
            if ids & ADOMETDC and ids & PLP:
                row["class"] = "AdoMetDC-PLP fusion"
            elif ids & ADOMETDC:
                row["class"] = "standalone AdoMetDC"
            elif ids & PLP:
                best, other, odc_res = nearest(seq, refs)
                k = sum(odc_res.get(int(l[1:])) == l[0] for l in ODC_SITE_LABELS)
                fam = best["label"]
                cls = ("ODC family" if fam.startswith(ODC_FAMILY) else
                       "ODC-fold ADC (PBCV-1-like)" if "PBCV" in fam else
                       "ADC (SpeA/plant ADC fold)" if fam.startswith(("ADC1", "ADC2", "ADC SpeA")) else
                       "DAPDC" if fam.startswith(("DAPDC", "LysA")) else
                       "LDC" if fam.startswith(("LDC", "LDC2")) else
                       "CANSDC" if fam.startswith("CANSDC") else "antizyme-inhibitor-like" if fam.startswith("AZIN") else fam)
                row.update({"class": cls, "nearest": f"{fam} ({best['ident']:.2f}; {best['score'] - other['score']:.0f} vs {other['label']})",
                            "odc_sites_seq": f"{k}/13 " + site_string(odc_res),
                            "K69_C360_seq": f"{odc_res.get(69, '-')}{odc_res.get(360, '-')}"})
                acc, path = af_model(accs, seq)
                if path:
                    o = odc_structural(odc, path, 1)
                    c = Counter(v["state"] for v in o["sites"].values())
                    row["structure"] = (f"{acc}: TM {o['tm']:.2f}, {c['kept']}/{c['lost']}/{c['unresolved']} "
                                        f"(lost {' '.join(l + v['aa'] for l, v in o['sites'].items() if v['state'] == 'lost') or 'none'})")
                    row["odc_13_structural_kept"] = c["kept"]
                species_plp.setdefault(taxid, []).append((e["uniParcId"], seq, sp, cls))
            else:
                fam = next(ROUTE[x] for x in ROUTE if x in ids)
                row["class"] = fam
                if "PF00795" in ids:
                    a = LOCAL.align(seq, refs[NCPAH_REF])[0]
                    pr = [(i, j) for (qs, qe), (ts, te) in zip(*a.aligned) for i, j in zip(range(qs, qe), range(ts, te))]
                    row["nearest"] = f"identity to Arabidopsis NLP1 (NCPAH) {sum(seq[i] == refs[NCPAH_REF][j] for i, j in pr) / max(1, len(pr)):.2f} over {len(pr)} aa"
            sp_rows.append(row)
            rep.append(f"| {row['upi']} | {row['uniprot'] or '(none)'} | {row['length']} | {row['pfam']} | {row['class']} | "
                       f"{row.get('nearest', '')} | {row.get('odc_sites_seq', '')} | {row.get('structure', '')} |")
        rep.append("")
    scols = ["taxid", "species", "upi", "uniprot", "length", "pfam", "class", "nearest", "K69_C360_seq",
             "odc_sites_seq", "structure", "odc_13_structural_kept"]
    with open(RESULTS / "p10_species_enzymes.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, scols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sp_rows)

    # ---- 4. retained or rotting
    rep += ["", "## 4. Retained or rotting", "",
            "| ODC status | n | partner length | TM to human ODC | longest unaligned | "
            "mean pLDDT | % partner < 70 |", "|---|---|---|---|---|---|---|"]
    for st in ("degraded", "intact", "undetermined", "incomplete"):
        fs = [f for f in fus if f.get("odc_status") == st]
        if not fs:
            continue
        rng = lambda k: f"{min(f[k] for f in fs)}-{max(f[k] for f in fs)} (median {np.median([f[k] for f in fs]):.3g})"
        rep.append(f"| {st} | {len(fs)} | {rng('partner_len')} | {rng('TM_to_human_ODC')} | {rng('longest_unaligned')} | "
                   f"{rng('partner_plddt_mean')} | {rng('partner_pct_below_70')} |")
    rep += ["", "Per protein (degraded first):", "",
            "| accession | organism | set | status | partner aa | TM | aligned | longest gap | pLDDT | % < 70 |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for f in sorted([f for f in fus if "odc_status" in f], key=lambda f: (f["odc_status"] != "degraded", f["odc_status"], f["organism"])):
        rep.append(f"| {f['accession']} | *{f['organism']}* | {f['set']} | {f['odc_status']} | {f['partner_len']} | "
                   f"{f['TM_to_human_ODC']} | {f['odc_positions_aligned']} | {f['longest_unaligned']} | "
                   f"{f['partner_plddt_mean']} | {f['partner_pct_below_70']} |")

    # ---- 5. correlates across the 23
    link = {r["accession"]: r for r in p7.read("p7_linkers.tsv")}
    gm = {r["accession"]: r for r in p7.read("p7_gene_models.tsv")}
    put = defaultdict(dict)
    for r in p7.read("p8_putrescine_site.tsv"):
        if r["accession"] in our:
            put[r["accession"]][int(r["human_pos"])] = r["cand_aa"] + ("" if r["state"] != "unresolved" else "?")
    corr = []
    for f in [f for f in fus if f["set"] == "our 23"]:
        a = f["accession"]
        corr.append({"accession": a, "organism": f["organism"], "group": f["group"], "odc_status": f["odc_status"],
                     "adometdc_half": f["adometdc_half"], "linker_len": link[a]["linker_len"],
                     "linker_plddt_mean": link[a]["linker_plddt_mean"], "protein_len": arch[a]["length"],
                     "putrescine_sites": " ".join(f"{p}{put[a].get(p) or '-'}" for p in PUTRESCINE),
                     "gene_model_evidence": gm[a]["evidence_type"], "n_exons": gm[a]["n_exons"],
                     "junction_vs_introns": gm[a]["junction_vs_introns"]})
    with open(RESULTS / "p10_correlates_23.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(corr[0]), delimiter="\t")
        w.writeheader()
        w.writerows(corr)
    rep += ["", "## 5. Correlates across the 23 (descriptive; no tests)", "",
            "| accession | organism | group | ODC status | AdoMetDC half | linker aa | protein aa | putrescine sites | gene model | exons | junction |",
            "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in sorted(corr, key=lambda c: (c["odc_status"] != "degraded", c["odc_status"], c["group"])):
        rep.append(f"| {c['accession']} | *{c['organism']}* | {c['group']} | {c['odc_status']} | {c['adometdc_half']} | "
                   f"{c['linker_len']} | {c['protein_len']} | {c['putrescine_sites']} | {c['gene_model_evidence']} | "
                   f"{c['n_exons']} | {c['junction_vs_introns']} |")
    # ---- 2b. tree
    seqs, labels, tags = {}, {}, {}
    for f in fus:
        s0, s1 = f["_plp"]
        if not s0 or not s1 or int(s1) - int(s0) + 1 < MIN_REGION:
            continue
        seqs[f["accession"]] = f["_seq"][int(s0) - 1:int(s1)]
        labels[f["accession"]] = f"{f['organism']}|{f['group']}|fusion|{f.get('odc_status', 'nomodel')}"
        tags[f["accession"]] = f
    excluded = [f["accession"] for f in fus if f["accession"] not in seqs]
    for acc, (label, _) in {**p6.REFS, **EXTRA_REFS}.items():
        if acc == "Q8IJ77" or acc in seqs:
            continue
        e = uniparc_search(f"uniprotkb:{acc}")
        sp = pfam_spans(e[0]) if e else []
        pl = [s for s in sp if s[0] in PLP]
        r = refs[acc]
        seqs[acc] = r[min(s[1] for s in pl) - 1:max(s[2] for s in pl)] if pl else r
        labels[acc] = f"{label}|reference"
    for taxid, lst in species_plp.items():
        for upi, s, sp, cls in lst:
            pl = [x for x in sp if x[0] in PLP]
            seqs[upi] = s[min(x[1] for x in pl) - 1:max(x[2] for x in pl)]
            labels[upi] = f"{SPECIES[taxid]}|standalone|{cls}"
    t = tree(seqs, labels)
    if t is None:
        (RESULTS / "p10_report_smoke.md").write_text("\n".join(rep) + "\n", encoding="utf8")
        print("\n".join(rep))
        print("smoke run: pyfamsa not installed, stopped before the tree")
        return
    ids, splits, counts, n, L = t
    idx = {a: k for k, a in enumerate(ids)}

    def support(members):
        s = canon({idx[m] for m in members if m in idx}, n)
        in_tree = s in splits
        return in_tree, 100 * counts[s] / N_BOOT

    def smallest_clade(members):
        m = {idx[x] for x in members if x in idx}
        best = None
        for s in splits:
            for side in (s, frozenset(range(n)) - s):
                if m <= side and (best is None or len(side) < len(best)):
                    best = side
        return [] if best is None else sorted(ids[k] for k in best - m)

    fus_in = [a for a in ids if a in tags]
    sets = {"red-algal fusions": [a for a in fus_in if tags[a]["group"] == "Rhodophyta"],
            "red-algal fusions + Percolomonas": [a for a in fus_in if tags[a]["group"] in ("Rhodophyta", "Heterolobosea")],
            "degraded ODC halves": [a for a in fus_in if tags[a].get("odc_status") == "degraded"],
            "green-lineage fusions": [a for a in fus_in if tags[a]["group"] == "Viridiplantae"],
            "Plasmodium/apicomplexan fusions": [a for a in fus_in if tags[a]["group"] == "Alveolata"]}
    rep += ["### Distance tree of the ODC half", "",
            f"{n} sequences, {L} alignment columns kept (<= 50% gaps), NJ on Kimura distances, {N_BOOT} bootstraps. "
            f"Fusions left out (PLP region < {MIN_REGION} aa or no PLP span): {', '.join(excluded) or 'none'}. "
            "Tree with bootstrap percentages: `p10_odc_tree.nwk`.", "",
            "| set | n | a clade in the NJ tree? | bootstrap % | others inside the smallest clade holding all of them |",
            "|---|---|---|---|---|"]
    trows = []
    for k, mem in sets.items():
        if len(mem) < 2:
            rep.append(f"| {k} | {len(mem)} | (fewer than 2) | | |")
            continue
        ok, sup = support(mem)
        intr = smallest_clade(mem)
        rep.append(f"| {k} | {len(mem)} | {'yes' if ok else 'no'} | {sup:.0f} | "
                   + (", ".join(f"{x} ({labels[x].split('|')[0]})" for x in intr[:12]) + (" ..." if len(intr) > 12 else "") or "none") + " |")
        trows.append({"set": k, "members": " ".join(mem), "clade_in_nj_tree": ok, "bootstrap_pct": round(sup, 1),
                      "intruders_in_smallest_clade": " ".join(intr)})
    with open(RESULTS / "p10_tree_support.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, list(trows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(trows)

    (RESULTS / "p10_report.md").write_text("\n".join(rep) + "\n", encoding="utf8")
    print("\n".join(rep))
    print("\nwrote results/p10_report.md and p10_*.tsv, p10_odc_tree.nwk, p10_odc_alignment.fasta")


if __name__ == "__main__":
    sys.exit(main())
