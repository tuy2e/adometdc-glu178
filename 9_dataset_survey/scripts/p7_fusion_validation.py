"""Track B checks: are the AdoMetDC - PLP-decarboxylase fusions real genes, and are the partner domains intact?

For the 23 proteins in results/p6_decarboxylase_partner.tsv:

A. Gene-model source. EMBL/ENA protein records behind each UniProt or UniParc entry: genome (WGS / CON /
   genomic DNA) or transcriptome (TSA / transcribed RNA). For genome records, the CDS exon structure
   (ENA join(...) location) is converted to protein coordinates, so the exon boundaries around the
   AdoMetDC / decarboxylase junction can be read.
B. Same species: active UniProt entries for the same taxon with an AdoMetDC Pfam (PF01536/PF02675) and/or
   a PLP-decarboxylase Pfam (PF02784/PF00278), classified as standalone AdoMetDC, standalone decarboxylase,
   or fusion. Deleted (inactive) entries are not searchable, so absence here is not absence from the genome.
C. Linker: residues between the end of the last AdoMetDC Pfam match and the start of the first PLP Pfam
   match, with its length, sequence and AlphaFold pLDDT. P. falciparum PfAdoMetDC/ODC (Q8IJ77) is measured
   the same way for comparison (its InterPro matches are fetched).
D. Partner-domain integrity, by sequence and by structure, against human ODC (P11926), whose UniProt-
   annotated functional residues are read from the alignment: K69 (PLP Schiff base), H197, S200, G237,
   E274-R277 (PLP), D331, D332 (substrate), C360 (catalytic), D361, Y389. Structure: the partner region of
   each AlphaFold model is TM-aligned to the AlphaFold model of human ODC (residues 1-422; the C-terminal
   tail is disordered); TM-score normalised by ODC, CA offset of each functional residue, pLDDT. The same
   is done for controls: T. brucei ODC (active), human AZIN1 (inactive ODC paralog), human AZIN2.
   Sequence-only extra reference: PBCV-1 arginine decarboxylase (no AlphaFold model).

Reads Marc's models folder read-only. Caches ENA records in data/ena/, AlphaFold reference models in
data/references/. Writes results/p7_gene_models.tsv, p7_same_species.tsv, p7_linkers.tsv,
p7_partner_integrity.tsv.
Run (a few minutes the first time, mostly ENA and UniProt requests):  py -3.12 p7_fusion_validation.py
"""
import csv
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from tmtools import tm_align
from foldtrace.io import load_structure

from p2_fetch_annotations import DATA, RESULTS, get, read_candidates, save_json

RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
ENA, REFS = DATA / "ena", DATA / "references"
ADOMETDC, PLP = {"PF01536", "PF02675"}, {"PF02784", "PF00278"}
ODC_SITES = {69: "K69 PLP Schiff base", 197: "H197", 200: "S200 PLP", 237: "G237 PLP", 274: "E274 PLP",
             275: "P275 PLP", 276: "G276 PLP", 277: "R277 PLP", 331: "Y331 substrate", 332: "D332 substrate",
             360: "C360 catalytic", 361: "D361 substrate", 389: "Y389"}
CONTROLS = {"P07805": "ODC T. brucei (active)", "O14977": "AZIN1 human (inactive)", "Q96A70": "AZIN2 human"}
ODC_ORDERED = 422


def read(name):
    with open(RESULTS / name, encoding="utf8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, rows):
    with open(RESULTS / name, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def plddt(path):
    out = {}
    with open(path, encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                out[int(ln[22:26])] = float(ln[60:66])
    return out


# ------------------------------------------------------------------ A. gene models
def ena_record(pid):
    ENA.mkdir(parents=True, exist_ok=True)
    p = ENA / f"{pid}.embl"
    if not p.exists():
        body = ""
        for attempt in range(4):      # ENA answers 406 to Accept: application/json, so plain request here
            try:
                req = urllib.request.Request(f"https://www.ebi.ac.uk/ena/browser/api/embl/{pid}",
                                             headers={"Accept": "text/plain", "User-Agent": "adometdc-dataset-survey/1.0"})
                body = urllib.request.urlopen(req, timeout=90).read().decode("utf8", "replace")
                break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    break
                time.sleep(2 ** attempt)
            except OSError:
                time.sleep(2 ** attempt)
        p.write_text(body, encoding="utf8")
    return p.read_text(encoding="utf8")


def parse_cds(rec):
    idl = next((l for l in rec.splitlines() if l.startswith("ID ")), "")
    ft, grab = [], False
    for l in rec.splitlines():
        if l.startswith("FT   CDS"):
            grab = True
            ft.append(l[21:].strip())
        elif grab and l.startswith("FT                   ") and not l[21:].startswith("/"):
            ft.append(l[21:].strip())
        elif grab:
            break
    loc = "".join(ft)
    spans = [(int(a), int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc)]
    lens = [b - a + 1 for a, b in spans]
    if loc.startswith("complement"):
        lens = lens[::-1]
    return idl, loc, lens


def exon_boundaries_protein(lens):
    """Protein residue after which each intron falls (cumulative CDS nt / 3)."""
    out, c = [], 0
    for L in lens[:-1]:
        c += L
        out.append(round(c / 3, 1))
    return out


def sources(acc, src):
    """[(protein_id, kind)] from the UniProt entry (active) or UniParc databases (inactive)."""
    res = []
    up = json.loads((DATA / "uniprot" / f"{acc}.json").read_text(encoding="utf8"))
    if up["status"] == "active":
        for x in up["entry"].get("uniProtKBCrossReferences", []):
            if x["database"] == "EMBL":
                pr = {p["key"]: p["value"] for p in x.get("properties", [])}
                res.append((pr.get("ProteinId", "").split(".")[0] + "." + pr.get("ProteinId", ".").split(".")[-1]
                            if pr.get("ProteinId") else "", pr.get("MoleculeType", "")))
        return res, up["entry"].get("extraAttributes", {}).get("uniParcId", "")
    upi = ((up.get("entry") or {}).get("extraAttributes") or {}).get("uniParcId")
    p = DATA / "uniparc" / f"{upi}_databases.json"
    if not p.exists():
        st, _, body = get(f"https://rest.uniprot.org/uniparc/{upi}/databases?format=json&size=100")
        p.write_text(body, encoding="utf8")
    for x in json.loads(p.read_text(encoding="utf8")).get("results", []):
        if x["database"].startswith("EMBL"):
            res.append((f"{x['id']}.{x.get('version', 1)}", x["database"]))
    return res, upi


def kind_of(mol, idline):
    s = f"{mol} {idline}"
    if "TSA" in s or "Transcribed_RNA" in s or "transcribed RNA" in s:
        return "transcriptome (TSA)"
    if "mRNA" in s:
        return "mRNA/cDNA"
    if "Genomic" in s or "genomic" in s or "WGS" in s or "CON" in s:
        return "genome"
    return "unknown"


# ------------------------------------------------------------------ D. integrity helpers
def af_model(acc):
    p = REFS / f"AF-{acc}.pdb"
    if not p.exists():
        st, _, body = get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
        url = json.loads(body)[0]["pdbUrl"]
        st, _, pdb = get(url)
        p.write_text(pdb, encoding="utf8")
    return p


def sub_structure(s, lo, hi):
    idx = [i for i, n in enumerate(s.resnums) if lo <= n <= hi]
    return s.ca[idx], "".join(s.seq[i] for i in idx), [s.resnums[i] for i in idx]


def structural(odc, cand_path, lo, hi):
    ca, seq, nums = sub_structure(load_structure(str(cand_path)), lo, hi)
    oca, oseq, onums = odc
    r = tm_align(oca, ca, oseq, seq)
    corr, i, j = {}, 0, 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-":
            corr[onums[i]] = j
        i += x != "-"
        j += y != "-"
    on = oca @ np.asarray(r.u).T + np.asarray(r.t)
    out = {}
    for pos in ODC_SITES:
        if pos in corr:
            jj = corr[pos]
            off = float(np.linalg.norm(on[onums.index(pos)] - ca[jj]))
            out[pos] = (seq[jj], nums[jj], off)
        else:
            out[pos] = ("-", "", None)
    # longest run of ODC positions with no aligned partner residue (internal deletion)
    runs, cur = [], 0
    for n in onums:
        if n in corr:
            if cur:
                runs.append(cur)
            cur = 0
        else:
            cur += 1
    return float(r.tm_norm_chain1), len(corr), out, max(runs) if runs else 0


def main():
    seqs = read_candidates()
    p6 = {r["accession"]: r for r in read("p6_decarboxylase_partner.tsv")}
    arch = {r["accession"]: r for r in read("architecture_per_protein.tsv")}
    master = {r["accession"]: r for r in read("master_candidates.tsv")}
    refseq = {r.id.split("|")[1]: str(r.seq) for r in SeqIO.parse(str(REFS / "decarboxylases.fasta"), "fasta")}

    def domains(a):
        return [(d.split()[0], *map(int, d.split()[1].split("-"))) for d in arch[a]["domains_with_positions"].split("; ")]

    gm_rows, sp_rows, ln_rows, it_rows = [], [], [], []
    for a in p6:
        d = domains(a)
        ad_end = max(e for pf, s, e in d if pf in ADOMETDC)
        plp_start = min(s for pf, s, e in d if pf in PLP and s > ad_end)
        # A. gene models
        srcs, upi = sources(a, arch[a]["source"])
        kinds, exon_txt, n_ex, junction = [], [], "", ""
        for pid, mol in srcs:
            if not pid or pid == ".":
                continue
            rec = ena_record(pid)
            idl, loc, lens = parse_cds(rec)
            k = kind_of(mol, idl)
            kinds.append(f"{pid}:{k}")
            if k == "genome" and lens and not exon_txt:
                b = exon_boundaries_protein(lens)
                n_ex = len(lens)
                exon_txt = ",".join(str(x) for x in b)
                near = [x for x in b if ad_end - 30 <= x <= plp_start + 30]
                junction = ("intron(s) at residue " + ",".join(map(str, near)) + " within 30 aa of the junction"
                            if near else ("domains on different exons, no intron within 30 aa of the junction"
                                          if any(ad_end < x < plp_start for x in b) else
                                          "junction inside one exon" if b else "single-exon CDS"))
        up_tax = ((json.loads((DATA / "uniprot" / f"{a}.json").read_text(encoding="utf8")).get("entry") or {})
                  .get("organism") or {}).get("taxonId")
        gm_rows.append({"accession": a, "organism": master[a]["organism"],
                        "taxid": str(up_tax or master[a]["taxid"]),   # UniProt's own taxid when the entry is active
                        "uniprot_status": arch[a]["uniprot_status"], "uniparc": upi,
                        "records": "; ".join(kinds) or "(none found)",
                        "evidence_type": "; ".join(sorted({k.split(":", 1)[1] for k in kinds})) or "unknown",
                        "n_exons": n_ex, "intron_positions_protein": exon_txt,
                        "adometdc_end": ad_end, "plp_start": plp_start, "junction_vs_introns": junction})
        # C. linker
        pl = plddt(RUN / "models" / f"{a}.pdb")
        link = seqs[a][ad_end:plp_start - 1]
        ln_rows.append({"accession": a, "organism": master[a]["organism"], "adometdc_end": ad_end,
                        "plp_start": plp_start, "linker_len": len(link), "linker_seq": link,
                        "linker_plddt_mean": round(np.mean([pl[n] for n in range(ad_end + 1, plp_start) if n in pl]), 1) if link else "",
                        "frac_plddt_below_50": round(np.mean([pl[n] < 50 for n in range(ad_end + 1, plp_start) if n in pl]), 2) if link else "",
                        "frac_D_E_K_R_S": round(sum(c in "DEKRS" for c in link) / len(link), 2) if link else "",
                        "frac_N_Q": round(sum(c in "NQ" for c in link) / len(link), 2) if link else ""})

    # B. same species (UniProt, active entries only)
    done = set()
    for r in gm_rows:
        t = r["taxid"]
        if t in done:
            continue
        done.add(t)
        q = (f"taxonomy_id:{t} AND (xref:pfam-PF01536 OR xref:pfam-PF02675 OR xref:pfam-PF02784 OR xref:pfam-PF00278)")
        st, _, body = get("https://rest.uniprot.org/uniprotkb/search?query=" + re.sub(" ", "%20", q)
                          + "&format=tsv&fields=accession,length,xref_pfam,protein_name,organism_name&size=200")
        lines = body.strip().splitlines()[1:] if st == 200 else []
        found = []
        for l in lines:
            acc, L, pf, name, org = (l.split("\t") + [""] * 5)[:5]
            pfs = set(pf.strip(";").split(";")) - {""}
            cls = ("fusion" if pfs & ADOMETDC and pfs & PLP else "AdoMetDC only" if pfs & ADOMETDC
                   else "decarboxylase only")
            found.append(f"{acc} ({L} aa, {cls})")
        cl = [f.split(", ")[-1].rstrip(")") for f in found]
        sp_rows.append({"taxid": t, "organism": r["organism"],
                        "n_fusion": cl.count("fusion"), "n_standalone_adometdc": cl.count("AdoMetDC only"),
                        "n_standalone_decarboxylase": cl.count("decarboxylase only"),
                        "entries": "; ".join(found) or "(none active in UniProt)"})

    # Plasmodium linker for comparison
    st, _, body = get("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/Q8IJ77/")
    pf_rows = []
    for x in json.loads(body).get("results", []) if st == 200 and body.strip() else []:
        for prot in x["proteins"]:
            for loc in prot["entry_protein_locations"]:
                for fr in loc["fragments"]:
                    pf_rows.append((x["metadata"]["accession"], fr["start"], fr["end"]))
    if pf_rows:
        ae = max(e for pf, s, e in pf_rows if pf in ADOMETDC)
        ps = min(s for pf, s, e in pf_rows if pf in PLP and s > ae)
        ln_rows.append({"accession": "Q8IJ77 (reference)", "organism": "Plasmodium falciparum 3D7", "adometdc_end": ae,
                        "plp_start": ps, "linker_len": ps - ae - 1, "linker_seq": refseq["Q8IJ77"][ae:ps - 1],
                        "linker_plddt_mean": "", "frac_plddt_below_50": "",
                        "frac_D_E_K_R_S": round(sum(c in "DEKRS" for c in refseq["Q8IJ77"][ae:ps - 1]) / max(1, ps - ae - 1), 2),
                        "frac_N_Q": round(sum(c in "NQ" for c in refseq["Q8IJ77"][ae:ps - 1]) / max(1, ps - ae - 1), 2)})

    # D. integrity
    odc_s = load_structure(str(af_model("P11926")))
    odc = sub_structure(odc_s, 1, ODC_ORDERED)
    al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1, end_gap_score=0)

    def seq_sites(q):
        a_ = al.align(refseq["P11926"], q)[0]
        m = {}
        for (rs, re_), (qs, qe) in zip(*a_.aligned):
            for i, j in zip(range(rs, re_), range(qs, qe)):
                m[i + 1] = q[j]
        ident = sum(refseq["P11926"][i] == q[j] for (rs, re_), (qs, qe) in zip(*a_.aligned)
                    for i, j in zip(range(rs, re_), range(qs, qe)))
        return m, ident / ODC_ORDERED

    targets = [(a, RUN / "models" / f"{a}.pdb", plddt_start - 25 if (plddt_start := next(int(x["plp_start"]) for x in gm_rows if x["accession"] == a)) else 1, 10 ** 6, master[a]["organism"])
               for a in p6]
    targets += [(c, af_model(c), 1, 10 ** 6, lab) for c, lab in CONTROLS.items()]
    for a, path, lo, hi, lab in targets:
        tm, n_al, sites, maxgap = structural(odc, path, lo, hi)
        s = load_structure(str(path))
        q = "".join(c for c, n in zip(s.seq, s.resnums) if lo <= n <= hi)
        m, ident = seq_sites(q)
        pl = plddt(path)
        kept = [p for p in ODC_SITES if sites[p][0] == refseq["P11926"][p - 1] and sites[p][2] is not None and sites[p][2] <= 4]
        it_rows.append({"accession": a, "organism": lab, "partner_from": lo, "partner_len": len(q),
                        "TM_to_human_ODC": round(tm, 3), "n_ODC_positions_aligned": f"{n_al}/{ODC_ORDERED}",
                        "longest_unaligned_ODC_stretch": maxgap,
                        "seq_identity_to_ODC": round(ident, 3),
                        "partner_plddt_mean": round(np.mean([v for n, v in pl.items() if lo <= n <= hi]), 1),
                        "functional_residues_identical_in_place": f"{len(kept)}/{len(ODC_SITES)}",
                        **{ODC_SITES[p].split()[0]:
                           (f"{sites[p][0]}{sites[p][1]} ({sites[p][2]:.1f}A)" if sites[p][2] is not None else "-")
                           for p in ODC_SITES}})
    for name, rows in (("p7_gene_models.tsv", gm_rows), ("p7_same_species.tsv", sp_rows),
                       ("p7_linkers.tsv", ln_rows), ("p7_partner_integrity.tsv", it_rows)):
        write(name, rows)
        print(f"\n===== {name}")
        for r in rows:
            print("  " + " | ".join(f"{k}={str(v)[:70]}" for k, v in r.items()))


if __name__ == "__main__":
    main()
