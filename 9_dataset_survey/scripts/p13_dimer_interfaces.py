"""Track B, the dimer question: AlphaFold Server homodimer predictions for the p11 jobs.

Reads data/p11_af_results/fold_*.zip directly (no unzipping). Per job, the top-ranked model (highest
ranking_score of the five; checked, not assumed to be model_0), its full_data PAE matrix and pLDDT.

1. TRUST GATE first. Mean inter-chain PAE = mean of PAE[A,B] and PAE[B,A] over all residue pairs.
   PAE_GATE (10 A) or below = the arrangement may be interpreted; above = unsupported, no geometry is
   interpreted for that job. The cutoff is a constant below; raw values are always written.
   The same mean is also given per domain block (ODC(A)xODC(B), AdoMetDC(A)xAdoMetDC(B)) and over
   residues with pLDDT >= 70 only, but the gate is the whole-chain mean as specified.
2. Domains are found structurally: TM-align of chain A to 1D7K chain A (ODC) and to 1MSV chain A
   (AdoMetDC); a domain = the model residues aligned within 5 A after superposition, span min..max.
   Everything outside both spans (linker, termini, parasite insertions) is "other".
3. Buried surface: Shrake-Rupley SASA of each chain alone minus in the dimer, summed over both chains
   (as p12). Each buried residue is assigned to the region of its nearest partner residue (heavy atoms,
   within 8 A), giving BSA per region pair: ODC-ODC, AdoMetDC-AdoMetDC, AdoMetDC-ODC (either way round),
   other. Contacts = heavy atoms within 4.5 A (p12 cutoff).
4. Real interface? The model's ODC(A)+ODC(B) pair (chains concatenated A then B) is TM-aligned onto
   1D7K A+B, TM-score normalised by the crystal dimer; the monomer TM-score (ODC(A) onto 1D7K A) is given
   beside it, so a dimer score near the monomer score means chain B sits where the crystal partner sits.
   Contact agreement: the model's ODC-ODC contact residues are put on human numbering through the monomer
   alignment and compared with the 1D7K contact set (same 4.5 A rule): recall = share of crystal interface
   positions also in contact in the model; precision = share of model contacts that are crystal interface
   positions. K69-C360: CA-CA distance of the model residues aligned to human K69 of one chain and C360 of
   the other, both ways round, against the same distance in 1D7K. AdoMetDC the same against 1MSV A+B.

Writes results/p13_dimer_interfaces.tsv and results/p13_dimer_interfaces.md
Run:  py -3.12 p13_dimer_interfaces.py                 (gate 10 A)
      py -3.12 p13_dimer_interfaces.py --gate 15       (a different cutoff, if you decide on one)
      --out DIR writes the two files to DIR instead of results/
"""
import io
import json
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.SeqUtils import seq1
from tmtools import tm_align

from p2_fetch_annotations import DATA, RESULTS
from p12_interface_vs_catalytic import pdb, protein_only

sys.stdout.reconfigure(encoding="utf8")
AF = DATA / "p11_af_results"
PAE_GATE = float(sys.argv[sys.argv.index("--gate") + 1]) if "--gate" in sys.argv else 10.0
OUT = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else None
CONTACT, ASSIGN, PLDDT_OK, ALIGN_D = 4.5, 8.0, 70, 5.0
JOBS = {  # job: (accession, organism, role, expected domains)
    "01": ("A0A1Y1HLG4", "Klebsormidium nitens", "fusion, ODC half intact", "both"),
    "02": ("A0A2V3IZY9", "Gracilariopsis chorda", "fusion, ODC half degraded", "both"),
    "03": ("A0A9C7PQ51", "Galdieria partita", "fusion, ODC half degraded", "both"),
    "04": ("C1E5E7", "Micromonas commoda", "fusion (Green 2011)", "both"),
    "05": ("Q8IJ77", "Plasmodium falciparum", "known bifunctional enzyme", "both"),
    "06": ("A0AAV8UQI0", "Rhodosorus marinus", "control: standalone ODC", "odc"),
    "07": ("A0A9D4U7K2", "Adiantum capillus-veneris", "control: standalone AdoMetDC", "ado"),
}


# ---------------------------------------------------------------- structure helpers
# 1D7K stores Lys69 as LLP (lysine carrying PLP, a HETATM record); it is kept as K, lysine atoms only
MODIFIED = {"LLP": "K"}
LYS_ATOMS = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"}


def residues(chain):
    return [r for r in chain if (r.id[0] == " " or r.get_resname() in MODIFIED) and "CA" in r]


def aa(r):
    return MODIFIED.get(r.get_resname()) or seq1(r.get_resname())


def ca_seq(res):
    return np.array([r["CA"].coord for r in res], float), "".join(aa(r) for r in res)


def heavy(res):
    return [a for r in res for a in r if a.element != "H"
            and (r.get_resname() not in MODIFIED or a.get_id() in LYS_ATOMS)]


def align_map(ref_res, mod_res):
    """TM-align ref onto model; return {ref resnum: model resnum} for pairs within ALIGN_D, and TM (ref-normalised)."""
    rx, rs = ca_seq(ref_res)
    mx, ms = ca_seq(mod_res)
    r = tm_align(rx, mx, rs, ms)
    on = rx @ np.asarray(r.u).T + np.asarray(r.t)
    m, i, j = {}, 0, 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-" and np.linalg.norm(on[i] - mx[j]) <= ALIGN_D:
            m[ref_res[i].id[1]] = mod_res[j].id[1]
        i += x != "-"
        j += y != "-"
    return m, float(r.tm_norm_chain1)


def dimer_tm(refA, refB, modA, modB):
    rx, rs = ca_seq(refA + refB)
    mx, ms = ca_seq(modA + modB)
    r = tm_align(rx, mx, rs, ms)
    # which model chain each crystal chain landed on (A->A, B->B expected; A->B means swapped, fine for a homodimer)
    nA = len(refA)
    land = defaultdict(int)
    i = j = 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-":
            land[("A" if i < nA else "B") + ("A" if j < len(modA) else "B")] += 1
        i += x != "-"
        j += y != "-"
    return float(r.tm_norm_chain1), dict(land)


def contacts(resA, resB, cut=CONTACT):
    """Residue pairs (a resnum, b resnum) with any heavy atoms within cut."""
    ns = NeighborSearch(heavy(resB))
    out = set()
    for r in resA:
        for a in r:
            if a.element == "H":
                continue
            for b in ns.search(a.coord, cut):
                out.add((r.id[1], b.get_parent().id[1]))
    return out


def bsa_per_residue(chA, chB):
    sr = ShrakeRupley()
    alone = {}
    for c in (chA, chB):
        s, m = Structure("x"), Model(0)
        s.add(m)
        cc = c.copy()
        m.add(cc)
        sr.compute(s, level="R")
        alone[c.id] = {r.id[1]: r.sasa for r in cc}
    s, m = Structure("d"), Model(0)
    s.add(m)
    a2, b2 = chA.copy(), chB.copy()
    m.add(a2)
    m.add(b2)
    sr.compute(s, level="R")
    return {c.id: {r.id[1]: alone[c.id][r.id[1]] - r.sasa for r in c} for c in (a2, b2)}


def nearest_partner(res, partner_res):
    ns = NeighborSearch(heavy(partner_res))
    best = {}
    for r in res:
        d0, p0 = None, None
        for a in r:
            for b in ns.search(a.coord, ASSIGN):
                d = float(np.linalg.norm(a.coord - b.coord))
                if d0 is None or d < d0:
                    d0, p0 = d, b.get_parent().id[1]
        best[r.id[1]] = p0
    return best


def by_num(chain, n):
    return next((r for r in residues(chain) if r.id[1] == n), None)


def ca_dist(chX, nx, chY, ny):
    x, y = by_num(chX, nx), by_num(chY, ny)
    if x is None or y is None:
        return None
    return round(float(np.linalg.norm(x["CA"].coord - y["CA"].coord)), 1)


# ---------------------------------------------------------------- crystal references
def reference(pid):
    st = pdb(pid)
    A, B = st[0]["A"], st[0]["B"]
    ra, rb = residues(A), residues(B)
    iface = {a for a, _ in contacts(ra, rb)} | {a for a, _ in contacts(rb, ra)}
    return {"A": A, "B": B, "ra": ra, "rb": rb, "iface": iface}


# ---------------------------------------------------------------- one job
def top_model(zf, stem):
    scores = []
    for i in range(5):
        s = json.loads(zf.read(f"{stem}_summary_confidences_{i}.json"))
        scores.append((s["ranking_score"], i, s))
    scores.sort(reverse=True)
    return scores[0][1], scores[0][2], [x[2] for x in scores]


def block_mean(pae, idx_r, idx_c):
    if not len(idx_r) or not len(idx_c):
        return None
    return float(pae[np.ix_(idx_r, idx_c)].mean())


def r1(x):
    return "" if x is None else round(x, 1)


def r2(x):
    return "" if x is None else round(x, 2)


def job(path, refs):
    zf = zipfile.ZipFile(path)
    stem = path.stem
    n = re.match(r"fold_(\d\d)", stem).group(1)
    acc, org, role, dom = JOBS[n]
    k, summ, allsum = top_model(zf, stem)
    st = MMCIFParser(QUIET=True).get_structure(stem, io.StringIO(zf.read(f"{stem}_model_{k}.cif").decode("utf8")))
    A, B = st[0]["A"], st[0]["B"]
    ra, rb = residues(A), residues(B)
    plddt = {c.id: {r.id[1]: float(np.mean([a.get_bfactor() for a in r])) for r in residues(c)} for c in (A, B)}

    full = json.loads(zf.read(f"{stem}_full_data_{k}.json"))
    pae = np.asarray(full["pae"], float)
    tch, tres = full["token_chain_ids"], full["token_res_ids"]
    row = {"job": n, "accession": acc, "organism": org, "role": role, "top_model": k,
           "length_per_chain": len(ra), "iptm": summ["iptm"], "ptm": summ["ptm"],
           "iptm_range_5_models": f"{min(s['iptm'] for s in allsum)}-{max(s['iptm'] for s in allsum)}",
           "pae_min_interchain": round(min(summ["chain_pair_pae_min"][0][1], summ["chain_pair_pae_min"][1][0]), 2)}

    # domains, from chain A
    # AdoMetDC first; ODC is then searched only outside the AdoMetDC span (TM-align of a whole fusion chain
    # onto ODC can otherwise pick up AdoMetDC residues)
    spans = {}
    maps = {}
    for key, pid in (("ado", "1MSV"), ("odc", "1D7K")):
        if dom not in (key, "both"):
            continue
        pool = [r for r in ra if not ("ado" in spans and spans["ado"][0] <= r.id[1] <= spans["ado"][1])]
        m, tm = align_map(refs[pid]["ra"], pool)
        if len(m) >= 30:
            spans[key] = (min(m.values()), max(m.values()))
            maps[key] = (m, tm)
    region = lambda num: next((k_ for k_, (lo, hi) in spans.items() if lo <= num <= hi), "other")
    for key in ("odc", "ado"):
        row[f"{key}_span"] = f"{spans[key][0]}-{spans[key][1]}" if key in spans else ""
        row[f"{key}_monomer_tm"] = round(maps[key][1], 3) if key in maps else ""

    # 1. trust
    idx = {c: [i for i, x in enumerate(tch) if x == c] for c in ("A", "B")}
    row["pae_mean_interchain"] = r2((block_mean(pae, idx["A"], idx["B"]) + block_mean(pae, idx["B"], idx["A"])) / 2)
    ordered = {c: [i for i in idx[c] if plddt[c].get(tres[i], 0) >= PLDDT_OK] for c in ("A", "B")}
    row["pae_mean_interchain_plddt70"] = r2((block_mean(pae, ordered["A"], ordered["B"])
                                              + block_mean(pae, ordered["B"], ordered["A"])) / 2)
    for key in ("odc", "ado"):
        if key in spans:
            lo, hi = spans[key]
            ia = [i for i in idx["A"] if lo <= tres[i] <= hi]
            ib = [i for i in idx["B"] if lo <= tres[i] <= hi]
            row[f"pae_{key}{key}_block"] = r2((block_mean(pae, ia, ib) + block_mean(pae, ib, ia)) / 2)
        else:
            row[f"pae_{key}{key}_block"] = ""
    ok = row["pae_mean_interchain"] <= PAE_GATE
    row["trust_gate"] = "pass" if ok else "FAIL: unsupported, geometry not interpreted"
    if not ok:
        return row

    # 2. contacts and buried surface by region pair
    bsa = bsa_per_residue(protein_only(A), protein_only(B))
    pair_bsa = defaultdict(float)
    for c, own, other in (("A", ra, rb), ("B", rb, ra)):
        near = nearest_partner(own, other)
        for num, b in bsa[c].items():
            if b <= 0:
                continue
            p = near.get(num)
            if p is None:
                pair_bsa["unassigned"] += b
                continue
            pr = tuple(sorted((region(num), region(p))))
            pair_bsa[{("odc", "odc"): "odc_odc", ("ado", "ado"): "ado_ado",
                      ("ado", "odc"): "ado_odc"}.get(pr, "other")] += b
    row["bsa_total_A2"] = round(sum(sum(v.values()) for v in bsa.values()))
    for kk in ("odc_odc", "ado_ado", "ado_odc", "other", "unassigned"):
        row[f"bsa_{kk}_A2"] = round(pair_bsa.get(kk, 0.0))
    cab = contacts(ra, rb)
    for kk, test in (("odc_odc", ("odc", "odc")), ("ado_ado", ("ado", "ado"))):
        cc = [(a, b) for a, b in cab if (region(a), region(b)) == test]
        resset = {("A", a) for a, _ in cc} | {("B", b) for _, b in cc}
        row[f"contact_residues_{kk}"] = len(resset)
        row[f"contact_plddt_{kk}"] = r1(np.mean([plddt[c][x] for c, x in resset])) if resset else ""

    # 3. real interface?
    for key, pid in (("odc", "1D7K"), ("ado", "1MSV")):
        cols = [f"{key}_dimer_tm", f"{key}_dimer_chain_landing", f"{key}_iface_recall", f"{key}_iface_precision"]
        if key not in spans:
            row.update({c: "" for c in cols})
            continue
        ref = refs[pid]
        lo, hi = spans[key]
        ma = [r for r in ra if lo <= r.id[1] <= hi]
        mb = [r for r in rb if lo <= r.id[1] <= hi]
        tm, land = dimer_tm(ref["ra"], ref["rb"], ma, mb)
        m = maps[key][0]
        back = {v: k_ for k_, v in m.items()}
        mod_iface = {back[a] for a, _ in contacts(ma, mb) if a in back} | {back[a] for a, _ in contacts(mb, ma) if a in back}
        comparable = {p for p in ref["iface"] if p in m}
        row[f"{key}_dimer_tm"] = round(tm, 3)
        row[f"{key}_dimer_chain_landing"] = " ".join(f"{k_}:{v}" for k_, v in sorted(land.items()))
        row[f"{key}_iface_recall"] = f"{len(mod_iface & comparable)}/{len(comparable)}"
        row[f"{key}_iface_precision"] = f"{len(mod_iface & ref['iface'])}/{len(mod_iface)}"
    if "odc" in spans:
        m = maps["odc"][0]
        k69, c360 = m.get(69), m.get(360)
        row["model_res_at_K69"] = f"{aa(by_num(A, k69))}{k69}" if k69 else "not aligned"
        row["model_res_at_C360"] = f"{aa(by_num(A, c360))}{c360}" if c360 else "not aligned"
        row["K69A_C360B_CA_A"] = ca_dist(A, k69, B, c360) if k69 and c360 else ""
        row["K69B_C360A_CA_A"] = ca_dist(B, k69, A, c360) if k69 and c360 else ""
    return row


def clean(chain):
    """Copy of a chain with only the residues() set (LLP reduced to its lysine atoms), for crystal buried surface."""
    c = chain.copy()
    keep = {r.id for r in residues(chain)}
    for r in list(c):
        if r.id not in keep:
            c.detach_child(r.id)
        elif r.get_resname() in MODIFIED:
            for a in list(r):
                if a.get_id() not in LYS_ATOMS:
                    r.detach_child(a.get_id())
    return c


def md(rows, crys, path):
    passed = [r for r in rows if r["trust_gate"] == "pass"]
    failed = [r for r in rows if r["trust_gate"] != "pass"]
    L = ["# p13: AlphaFold Server homodimers of the p11 jobs",
         "",
         "Generated by `scripts/p13_dimer_interfaces.py`. One seed per job (AlphaFold Server, pilot); the top-ranked "
         "of five models is analysed. Method notes are in the script docstring.",
         "",
         f"## 1. Trust gate (mean inter-chain PAE, cutoff {PAE_GATE:g} A)",
         "",
         "Whole-chain = mean of PAE over every residue of chain A against every residue of chain B, both ways. "
         "The block columns restrict the same mean to one domain pair; the pLDDT >= 70 column drops disordered residues. "
         "The gate is applied to the whole-chain column only.",
         "",
         "| job | protein | role | ipTM (5 models) | min inter-chain PAE | mean inter-chain PAE | mean, pLDDT>=70 | ODC x ODC block | AdoMetDC x AdoMetDC block | gate |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['job']} | {r['accession']} *{r['organism']}* | {r['role']} | {r['iptm']} ({r['iptm_range_5_models']}) | "
                 f"{r['pae_min_interchain']} | {r['pae_mean_interchain']} | {r['pae_mean_interchain_plddt70']} | "
                 f"{r['pae_odcodc_block'] or '-'} | {r['pae_adoado_block'] or '-'} | "
                 f"{'pass' if r['trust_gate'] == 'pass' else '**unsupported**'} |")
    L += ["", f"**Cannot be interpreted (mean inter-chain PAE above {PAE_GATE:g} A):** "
          + (", ".join(f"{r['job']} *{r['organism']}* ({r['pae_mean_interchain']} A)" for r in failed) or "none")
          + ". For these, no contact, buried-surface or superposition result is reported.", ""]
    L += ["## 2-4. Interfaces (jobs that passed the gate only)", ""]
    if not passed:
        L += ["No job passed the gate, so no geometry is reported.", ""]
    else:
        L += ["| job | BSA total | ODC-ODC | AdoMetDC-AdoMetDC | AdoMetDC-ODC | other | ODC dimer TM (monomer TM) | ODC interface recall / precision | K69(A)-C360(B), K69(B)-C360(A) CA | AdoMetDC dimer TM (monomer TM) | AdoMetDC interface recall / precision |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in passed:
            g = lambda k: r.get(k, "")
            L.append(f"| {r['job']} *{r['organism']}* | {g('bsa_total_A2')} | {g('bsa_odc_odc_A2')} | {g('bsa_ado_ado_A2')} | "
                     f"{g('bsa_ado_odc_A2')} | {g('bsa_other_A2')} | {g('odc_dimer_tm') or '-'} ({g('odc_monomer_tm') or '-'}) | "
                     f"{g('odc_iface_recall') or '-'} / {g('odc_iface_precision') or '-'} | "
                     f"{g('K69A_C360B_CA_A') or '-'}, {g('K69B_C360A_CA_A') or '-'} | {g('ado_dimer_tm') or '-'} ({g('ado_monomer_tm') or '-'}) | "
                     f"{g('ado_iface_recall') or '-'} / {g('ado_iface_precision') or '-'} |")
        L.append("")
    L += ["## Crystal references (same code, same 4.5 A contact rule)", "",
          f"- 1D7K (human ODC) chains A+B: {crys['odc_iface']} interface residues, buried surface {crys['odc_bsa']} A2 "
          f"(both chains); K69(A)-C360(B) CA {crys['K69A_C360B']} A, K69(B)-C360(A) CA {crys['K69B_C360A']} A. "
          "Lys69 is the PLP-bound LLP residue in this entry and is counted as lysine (PLP atoms excluded).",
          f"- 1MSV (human AdoMetDC S68A) chains A+B: {crys['ado_iface']} interface residues, buried surface {crys['ado_bsa']} A2.",
          ""]
    path.write_text("\n".join(L), encoding="utf8")


def main():
    refs = {pid: reference(pid) for pid in ("1D7K", "1MSV")}
    crys = {"K69A_C360B": ca_dist(refs["1D7K"]["A"], 69, refs["1D7K"]["B"], 360),
            "K69B_C360A": ca_dist(refs["1D7K"]["B"], 69, refs["1D7K"]["A"], 360),
            "odc_iface": len(refs["1D7K"]["iface"]), "ado_iface": len(refs["1MSV"]["iface"])}
    for key, pid in (("odc_bsa", "1D7K"), ("ado_bsa", "1MSV")):
        b = bsa_per_residue(clean(refs[pid]["A"]), clean(refs[pid]["B"]))
        crys[key] = round(sum(sum(v.values()) for v in b.values()))
    print("crystal:", crys)
    rows = []
    for p in sorted(AF.glob("fold_*.zip")):
        r = job(p, refs)
        print(r["job"], r["organism"], r["trust_gate"], r["pae_mean_interchain"])
        rows.append(r)
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    out = OUT or RESULTS
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "p13_dimer_interfaces.tsv", "w", encoding="utf8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    md(rows, crys, out / "p13_dimer_interfaces.md")


if __name__ == "__main__":
    main()
