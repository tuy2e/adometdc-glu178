"""p24: the fusion as its own reference - internal AdoMetDC / ODC-domain architecture of the AFDB fusion models.

Proteins: every row of results/p10_all_fusions.tsv whose model is "AlphaFold DB" (the fusions with a model in
data/p10_models/). The current AFDB model and its PAE file are downloaded together into data/p24/ (a PAE belongs to
one model version); the model's sequence must equal our UniParc sequence.

Domains (per protein, Pfam spans from p10_fusion_search.tsv): AdoMetDC = first..last PF01536/PF02675 residue;
ODC = PLP-decarboxylase part, first PLP Pfam start to the C-terminus when AdoMetDC comes first, N-terminus to the
last PLP Pfam end when the PLP domain comes first (as p20).

1. GATE. Inter-domain PAE = mean of PAE[AdoMetDC rows, ODC columns] and PAE[ODC rows, AdoMetDC columns] (Angstrom).
   Also given over residues with pLDDT >= 70 only. Pass = mean inter-domain PAE <= PAE_GATE (10 A, the cutoff used in
   p13). Only passing models are interpreted below; numbers for failing models are written to the TSV but marked.
   Contacts: heavy atoms within 4.5 A between the two domains; buried surface = SASA(domain alone) summed for both
   domains minus SASA(both together), Shrake-Rupley.
2. Conservation (passing models only): for each pair, TM-align AdoMetDC domain of A onto that of B, apply the same
   rotation to A's ODC domain, then TM-align the two ODC domains only to get residue correspondence and report the
   CA RMSD of those corresponding ODC residues WITHOUT refitting, plus the ODC centroid displacement.
3. Dimer surfaces: ODC-domain contact residues mapped to human ODC numbering (TM-align onto the AlphaFold model of
   human ODC, residues 1-422, CA within 4 A) and compared with the 1D7K interface (p12_interface_residues.tsv,
   contact 4.5 A); AdoMetDC contact residues mapped onto 1MSV chain A and compared with the 1MSV interface.
   Expected overlap by chance = share of mapped positions that are interface positions.
4. Degraded vs intact: ODC-half status from p21b_rescore_p10.tsv (core rule).

Writes results/p24_architecture.tsv, results/p24_pairs.tsv (step 2), results/p24_report.md
"""
import csv
import itertools
import json
import statistics as st
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

import numpy as np
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from tmtools import tm_align

import p7_fusion_validation as p7
from foldtrace.io import load_structure
from p2_fetch_annotations import DATA, RESULTS
from p12_interface_vs_catalytic import pdb as rcsb_pdb
from p20_dnds_ci import regions_any

sys.stdout.reconfigure(encoding="utf8")
P24 = DATA / "p24"
P24.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "adometdc-dataset-survey/1.0"}
PAE_GATE, CONTACT, OFFSET, PLDDT_OK = 10.0, 4.5, 4.0, 70
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def get(url, path, binary=False):
    if not path.exists():
        for k in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                    data = r.read()
                break
            except (urllib.error.URLError, OSError):
                time.sleep(2 ** k)
        else:
            raise RuntimeError(url)
        path.write_bytes(data)
        time.sleep(0.3)
    return path


def afdb(acc):
    meta = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", P24 / f"{acc}.api.json").read_text())[0]
    pdbp = get(meta["pdbUrl"], P24 / f"AF-{acc}-v{meta['latestVersion']}.pdb")
    paep = get(meta["paeDocUrl"], P24 / f"AF-{acc}-v{meta['latestVersion']}.pae.json")
    pae = json.loads(paep.read_text())
    pae = pae[0] if isinstance(pae, list) else pae
    return meta, pdbp, np.asarray(pae["predicted_aligned_error"], float)


def residues_in(chain, lo, hi):
    return [r for r in chain if r.id[0] == " " and lo <= r.id[1] <= hi and "CA" in r]


def ca(res):
    return np.array([r["CA"].coord for r in res], float), "".join(THREE.get(r.get_resname(), "X") for r in res)


def sasa_sum(res_lists):
    s, m, c = Structure("x"), Model(0), Chain("A")
    s.add(m)
    m.add(c)
    for rl in res_lists:
        for r in rl:
            c.add(r.copy())
    ShrakeRupley().compute(s, level="S")
    return s.sasa


def mapping(res, ref_ca, ref_seq, ref_nums):
    """{residue number: reference position} via TM-align, CA within OFFSET after superposition."""
    x, s = ca(res)
    r = tm_align(ref_ca, x, ref_seq, s)
    on = ref_ca @ np.asarray(r.u).T + np.asarray(r.t)
    m, i, j = {}, 0, 0
    for a, b in zip(r.seqxA, r.seqyA):
        if a != "-" and b != "-" and np.linalg.norm(on[i] - x[j]) <= OFFSET:
            m[res[j].id[1]] = ref_nums[i]
        i += a != "-"
        j += b != "-"
    return m, float(r.tm_norm_chain1)


def main():
    search = {}
    for r in read("p10_fusion_search.tsv"):
        for k in [r["upi"]] + [a.split(".")[0] for a in r["uniprot_accessions"].split(",") if a]:
            search[k] = r
    ents = json.loads(next((DATA / "p10").glob("uniparc__feature_id_PF01536*.json")).read_text(encoding="utf8"))
    seq_of = {e["uniParcId"]: e["sequence"]["value"] for e in ents}
    core = {r["accession"].split(".")[0]: r["call_core"] for r in read("p21b_rescore_p10.tsv")}
    iface = {"ODC": set(), "AdoMetDC": set()}
    for r in read("p12_interface_residues.tsv"):
        if r["contact_4.5"] == "True":
            iface[r["enzyme"]].add(int(r["pos"]))
    odc_s = load_structure(str(p7.af_model("P11926")))
    hodc = p7.sub_structure(odc_s, 1, p7.ODC_ORDERED)
    msv = rcsb_pdb("1MSV")[0]["A"]
    msv_res = [r for r in msv if r.id[0] == " " and "CA" in r]
    msv_ca, msv_seq = ca(msv_res)
    msv_nums = [r.id[1] for r in msv_res]
    parser = PDBParser(QUIET=True)

    rows, keep = [], {}
    fus = [f for f in read("p10_all_fusions.tsv") if f["model"] == "AlphaFold DB"]
    for f in fus:
        acc = f["accession"].split(".")[0]
        s = search[acc]
        seq = seq_of[s["upi"]]
        meta, pdbp, pae = afdb(acc)
        ch = parser.get_structure(acc, str(pdbp))[0]["A"]
        mseq = "".join(THREE.get(r.get_resname(), "X") for r in ch if r.id[0] == " ")
        reg = regions_any(s["pfam_architecture"], len(seq))
        A = residues_in(ch, *reg["AdoMetDC"])
        O = residues_in(ch, *reg["ODC"])
        pl = {r.id[1]: r["CA"].get_bfactor() for r in ch if "CA" in r}
        ia = [r.id[1] - 1 for r in A]
        io = [r.id[1] - 1 for r in O]
        inter = (pae[np.ix_(ia, io)].mean() + pae[np.ix_(io, ia)].mean()) / 2
        ia70 = [i for i in ia if pl[i + 1] >= PLDDT_OK]
        io70 = [i for i in io if pl[i + 1] >= PLDDT_OK]
        inter70 = ((pae[np.ix_(ia70, io70)].mean() + pae[np.ix_(io70, ia70)].mean()) / 2) if ia70 and io70 else None
        intraA = pae[np.ix_(ia, ia)].mean()
        intraO = pae[np.ix_(io, io)].mean()
        ns = NeighborSearch([a for r in O for a in r if a.element != "H"])
        pairs = set()
        for r in A:
            for a in r:
                for b in ns.search(a.coord, CONTACT):
                    pairs.add((r.id[1], b.get_parent().id[1]))
        cA = sorted({a for a, _ in pairs})
        cO = sorted({b for _, b in pairs})
        bsa = sasa_sum([A]) + sasa_sum([O]) - sasa_sum([A, O]) if pairs else 0.0
        ok = inter <= PAE_GATE
        row = {"accession": f["accession"], "organism": f["organism"], "group": f["group"],
               "odc_core_call": core.get(acc, ""), "length": len(seq), "afdb_version": meta["latestVersion"],
               "model_matches_sequence": mseq == seq, "adometdc_domain": "%d-%d" % reg["AdoMetDC"],
               "odc_domain": "%d-%d" % reg["ODC"], "linker": "%d-%d" % reg["linker"],
               "pae_interdomain_mean": round(inter, 2),
               "pae_interdomain_plddt70": round(inter70, 2) if inter70 is not None else "",
               "pae_intra_adometdc": round(intraA, 2), "pae_intra_odc": round(intraO, 2),
               "gate": "pass" if ok else "FAIL",
               "contact_residues_adometdc": len(cA), "contact_residues_odc": len(cO),
               "contact_plddt_mean": round(st.mean([pl[x] for x in cA + cO]), 1) if pairs else "",
               "bsa_between_domains_A2": round(bsa)}
        if pairs:
            mo, tmo = mapping(O, *hodc)
            ma, tma = mapping(A, msv_ca, msv_seq, msv_nums)
            ho = {mo[x] for x in cO if x in mo}
            ha = {ma[x] for x in cA if x in ma}
            row.update({"odc_contacts_mapped_to_human": len(ho), "odc_contacts_at_1D7K_interface": len(ho & iface["ODC"]),
                        "odc_chance_share": round(len(set(mo.values()) & iface["ODC"]) / max(1, len(set(mo.values()))), 3),
                        "adometdc_contacts_mapped_to_1MSV": len(ha), "adometdc_contacts_at_1MSV_interface": len(ha & iface["AdoMetDC"]),
                        "adometdc_chance_share": round(len(set(ma.values()) & iface["AdoMetDC"]) / max(1, len(set(ma.values()))), 3),
                        "odc_contact_positions_human": " ".join(map(str, sorted(ho))),
                        "adometdc_contact_positions_1MSV": " ".join(map(str, sorted(ha)))})
        rows.append(row)
        if ok:
            keep[acc] = (A, O, row)
        print(f"{acc:12s} {f['organism'][:30]:30s} PAE {inter:5.1f} ({row['gate']}) contacts {len(cA)}/{len(cO)}", flush=True)
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(RESULTS / "p24_architecture.tsv", "w", encoding="utf8", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)

    prs = []
    for a, b in itertools.combinations(keep, 2):
        Aa, Oa, ra = keep[a]
        Ab, Ob, rb = keep[b]
        xa, sa = ca(Aa)
        xb, sb = ca(Ab)
        t = tm_align(xa, xb, sa, sb)
        U, T = np.asarray(t.u), np.asarray(t.t)
        oa, soa = ca(Oa)
        ob, sob = ca(Ob)
        oa_on = oa @ U.T + T
        t2 = tm_align(oa, ob, soa, sob)
        i = j = 0
        d = []
        for x, y in zip(t2.seqxA, t2.seqyA):
            if x != "-" and y != "-":
                d.append(np.linalg.norm(oa_on[i] - ob[j]))
            i += x != "-"
            j += y != "-"
        prs.append({"a": a, "group_a": ra["group"], "b": b, "group_b": rb["group"],
                    "adometdc_tm": round(float(t.tm_norm_chain1), 3), "odc_tm": round(float(t2.tm_norm_chain1), 3),
                    "odc_rmsd_after_adometdc_superposition": round(float(np.sqrt(np.mean(np.square(d)))), 1) if d else "",
                    "odc_centroid_shift_A": round(float(np.linalg.norm(oa_on.mean(0) - ob.mean(0))), 1)})
    if prs:
        with open(RESULTS / "p24_pairs.tsv", "w", encoding="utf8", newline="") as fh:
            w = csv.DictWriter(fh, list(prs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(prs)
    report(rows, prs)


def report(rows, prs):
    n = len(rows)
    passed = [r for r in rows if r["gate"] == "pass"]
    pv = sorted(r["pae_interdomain_mean"] for r in rows)
    L = ["# p24: the fusion as its own reference - internal domain architecture", "",
         "Generated by `scripts/p24_architecture.py` (method in its docstring).", "",
         "## 1. The gate: inter-domain PAE", "",
         f"**{len(passed)} of {n} AFDB fusion models pass** (mean inter-domain PAE <= {PAE_GATE:g} A). Inter-domain PAE "
         f"across all {n}: median {st.median(pv):.1f} A, range {pv[0]:.1f}-{pv[-1]:.1f} A. For comparison, PAE within a "
         f"domain: AdoMetDC median {st.median(r['pae_intra_adometdc'] for r in rows):.1f} A, ODC median "
         f"{st.median(r['pae_intra_odc'] for r in rows):.1f} A.", "",
         "| PAE between the two domains | models |", "|---|---|"]
    for lo, hi in ((0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 40)):
        L.append(f"| {lo}-{hi} A | {sum(lo <= x < hi for x in pv)} |")
    L += ["", "| group | models | pass | inter-domain PAE, median (range) | with pLDDT >= 70 only, median |", "|---|---|---|---|---|"]
    for g, c in Counter(r["group"] for r in rows).most_common():
        rs = [r for r in rows if r["group"] == g]
        p70 = [r["pae_interdomain_plddt70"] for r in rs if r["pae_interdomain_plddt70"] != ""]
        L.append(f"| {g} | {len(rs)} | {sum(r['gate'] == 'pass' for r in rs)} | "
                 f"{st.median(r['pae_interdomain_mean'] for r in rs):.1f} ({min(r['pae_interdomain_mean'] for r in rs):.1f}-"
                 f"{max(r['pae_interdomain_mean'] for r in rs):.1f}) | {st.median(p70):.1f} |" if p70 else "| - |")
    if not passed:
        L += ["", "**No model passes, so the internal arrangement is not predicted confidently by these AFDB models and "
              "is not interpreted.** Steps 2-4 are therefore not reported. The contact and buried-surface columns in "
              "the TSV are what each model happens to show and should not be read as an architecture.", ""]
    else:
        L += ["", "## 1b. Contacts in the passing models", "",
              "| accession | organism | group | ODC call | PAE | contact residues AdoMetDC / ODC | buried surface (A2) | ODC contacts at 1D7K interface / mapped (chance share) | AdoMetDC contacts at 1MSV interface / mapped (chance share) |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in passed:
            L.append(f"| {r['accession']} | *{r['organism']}* | {r['group']} | {r['odc_core_call']} | {r['pae_interdomain_mean']} | "
                     f"{r['contact_residues_adometdc']} / {r['contact_residues_odc']} | {r['bsa_between_domains_A2']} | "
                     f"{r.get('odc_contacts_at_1D7K_interface', '-')} / {r.get('odc_contacts_mapped_to_human', '-')} ({r.get('odc_chance_share', '-')}) | "
                     f"{r.get('adometdc_contacts_at_1MSV_interface', '-')} / {r.get('adometdc_contacts_mapped_to_1MSV', '-')} ({r.get('adometdc_chance_share', '-')}) |")
        L.append("")
        if prs:
            L += ["## 2. Is the arrangement conserved?", "",
                  "AdoMetDC domains superposed; ODC-domain CA RMSD without refitting (low = same placement).", "",
                  "| comparison | pairs | ODC RMSD after AdoMetDC superposition, median (range) |", "|---|---|---|"]
            def grp(p):
                return " vs ".join(sorted([p["group_a"], p["group_b"]]))
            for g, c in Counter(map(grp, prs)).most_common():
                v = [p["odc_rmsd_after_adometdc_superposition"] for p in prs if grp(p) == g and p["odc_rmsd_after_adometdc_superposition"] != ""]
                L.append(f"| {g} | {c} | {st.median(v):.1f} ({min(v):.1f}-{max(v):.1f}) |" if v else f"| {g} | {c} | - |")
            L.append("")
    (RESULTS / "p24_report.md").write_text("\n".join(L) + "\n", encoding="utf8")


if __name__ == "__main__":
    main()
