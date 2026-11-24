"""Track B mechanism: are degraded ODC halves kept for structure (the dimer interface) rather than chemistry?

1. Interfaces derived from crystal structures, not assumed:
     ODC       human ODC 1D7K (2.1 A; chains A+B = the PDB/PISA biological dimer, assembly 1)
     AdoMetDC  human AdoMetDC S68A proenzyme 1MSV (1.75 A; chains A+B = assembly 1, PISA dimer)
   Protein atoms only (waters and ligands removed). Per residue: buried surface on dimerisation
   (Shrake-Rupley SASA of the chain alone minus SASA in the dimer), and contact = any heavy atom within
   4.5 A of the other chain (counts at 4.0 and 5.0 A also written). Interface set = residues contacting the
   partner at 4.5 A in either chain (numbers are the same in both chains of a homodimer; the union is used).
   ODC numbering is checked against UniProt P11926, AdoMetDC numbering against P17707.
2. Catalytic sets: ODC = the 13 human-ODC functional residues used in p7/p9; AdoMetDC = Marc's five sites.
   Residues in both sets are reported, and for the per-protein tests they are counted as catalytic only.
3. Mapping with the p9 gate: kept = identical, CA within 4 A after superposition, pLDDT >= 70; lost =
   substituted with the same two conditions; unresolved = otherwise. Fraction conserved = kept / (kept +
   lost), unresolved left out. ODC: TM-align of the partner region onto the AlphaFold model of human ODC
   (1-422), as p7; interface positions beyond 422 are dropped. AdoMetDC: p8.map_sites (1MSV chain A).
   Proteins: ODC side - our 23 fusions, the p10 fusions with an AlphaFold model, standalone ODC baseline
   (reviewed UniProt EC 4.1.1.17 with PF02784 and an AlphaFold model, plus the red-algal standalone
   ODC-family proteins from p10 with a model). AdoMetDC side - our 23, p10 fusions with a model, and the
   590 standalone AdoMetDCs.
4. Tests (no test on fewer than 3 proteins per group):
   per protein   Fisher exact, interface kept/lost vs catalytic kept/lost
   within group  exact sign test on d = interface fraction - catalytic fraction (one value per protein,
                 so no lineage assumption across proteins is needed)
   between groups  difference in mean d (or mean interface fraction), labels shuffled within taxon
                 stratum (p5 control; 10,000 permutations, seed 1); strata = p10 groups (Rhodophyta,
                 Viridiplantae, Alveolata, ...) on the ODC side, p5 phylum on the AdoMetDC side. The number of
                 strata holding both groups is reported; with none, no lineage-controlled p is possible.
   Test group on the ODC side = the degraded fusions of Rhodophyta and Heterolobosea (p10 status).

Reads Marc's models folder and data/p10_models read-only. Downloads to data/p12/.
Writes results/p12_interface_residues.tsv, p12_odc_per_protein.tsv, p12_adometdc_per_protein.tsv,
       p12_tests.tsv, p12_report.md
Run:  py -3.12 p12_interface_vs_catalytic.py           (AdoMetDC side reuses the p3/p8 alignment cache)
      py -3.12 p12_interface_vs_catalytic.py --smoke   (interfaces + our 23 only, no downloads of baselines)
"""
import copy
import csv
import json
import math
import random
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

import numpy as np
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from tmtools import tm_align
from foldtrace.io import load_structure

import p7_fusion_validation as p7
import p8_putrescine_site as p8
from p2_fetch_annotations import DATA, RESULTS, get
from p5_loss_vs_fusion import fisher
from p3_residue_grid import REF_PDB, REF_CHAIN, PLDDT_OK, RUN, human_sequence, plddt_by_resnum
from p9_self_processing import ODC_SITE_LABELS, marc_sites

sys.stdout.reconfigure(encoding="utf8")
SMOKE = "--smoke" in sys.argv
P12 = DATA / "p12"
MODELS = P12 / "models"
MODELS.mkdir(parents=True, exist_ok=True)
P10_MODELS = DATA / "p10_models"
CONTACT, N_PERM = 4.5, 10000
ODC_CAT = [int(l[1:]) for l in ODC_SITE_LABELS]
FOCUS = ["Rhodophyta", "Glaucocystophyceae", "Cryptophyceae", "Heterolobosea", "Haptophyta", "Stramenopiles",
         "Alveolata", "Rhizaria", "Viridiplantae", "Metazoa", "Fungi", "Euglenozoa", "Bacteria", "Archaea"]


# ---------------------------------------------------------------- 1. interfaces
def pdb(pid):
    p = P12 / f"{pid}.pdb"
    if not p.exists():
        st, _, body = get(f"https://files.rcsb.org/download/{pid}.pdb")
        p.write_text(body, encoding="utf8")
    return PDBParser(QUIET=True).get_structure(pid, str(p))


def protein_only(chain):
    # 1D7K stores Lys69 as LLP (lysine carrying PLP, HETATM); keep it as lysine atoms, drop the PLP atoms
    c = copy.deepcopy(chain)
    for r in list(c):
        if r.get_resname() == "LLP":
            for a in list(r):
                if a.get_id() not in ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"):
                    r.detach_child(a.get_id())
        elif r.id[0] != " ":
            c.detach_child(r.id)
    return c


def lone(chain):
    s, m = Structure("x"), Model(0)
    s.add(m)
    m.add(chain)
    return s


def interface(pid, ca, cb, uniprot_seq):
    st = pdb(pid)
    A, B = protein_only(st[0][ca]), protein_only(st[0][cb])
    sr = ShrakeRupley()
    alone = {}
    for c in (A, B):
        cc = copy.deepcopy(c)
        sr.compute(lone(cc), level="R")
        alone[c.id] = {r.id[1]: r.sasa for r in cc}
    dimer = Structure("d")
    m = Model(0)
    dimer.add(m)
    a2, b2 = copy.deepcopy(A), copy.deepcopy(B)
    m.add(a2)
    m.add(b2)
    sr.compute(dimer, level="R")
    inside = {c.id: {r.id[1]: r.sasa for r in c} for c in (a2, b2)}
    mism = [f"{c.id}{r.id[1]}{r.get_resname()}" for c in (A, B) for r in c
            if r.id[1] <= len(uniprot_seq) and seq1(r.get_resname()) not in (uniprot_seq[r.id[1] - 1], "X")]
    rows = defaultdict(dict)
    for c, other in ((A, B), (B, A)):
        ns = NeighborSearch([a for a in other.get_atoms() if a.element != "H"])
        for r in c:
            n = r.id[1]
            dmin = min((np.linalg.norm(a.coord - b.coord) for a in r if a.element != "H"
                        for b in ns.search(a.coord, 5.0)), default=None)
            rows[n][c.id] = {"resname": r.get_resname(), "bsa": alone[c.id][n] - inside[c.id][n], "dmin": dmin}
    out = []
    for n in sorted(rows):
        d = rows[n]
        bsa = [v["bsa"] for v in d.values()]
        dm = [v["dmin"] for v in d.values() if v["dmin"] is not None]
        out.append({"pos": n, "uniprot_aa": uniprot_seq[n - 1] if n <= len(uniprot_seq) else "",
                    "chains_present": "".join(sorted(d)), "bsa_mean_A2": round(float(np.mean(bsa)), 1),
                    "min_dist_to_partner_A": round(min(dm), 2) if dm else "",
                    "contact_4.0": any(x <= 4.0 for x in dm), "contact_4.5": any(x <= CONTACT for x in dm),
                    "contact_5.0": any(x <= 5.0 for x in dm)})
    total_bsa = sum(alone[c.id][n] - inside[c.id][n] for c in (A, B) for n in alone[c.id])
    return out, total_bsa, mism


def seq1(three):
    from Bio.SeqUtils import seq1 as s1
    return s1(three)


# ---------------------------------------------------------------- 2. mapping
def gate(aa, human_aa, off, p):
    if aa in (None, "-", ""):
        return "unresolved"
    if off > p8.OFFSET or p is None or p < PLDDT_OK:
        return "unresolved"
    return "kept" if aa == human_aa else "lost"


def map_odc(odc, odc_seq, path, lo, hi, positions):
    s = load_structure(str(path))
    ca, seq, nums = p7.sub_structure(s, lo, hi)
    oca, oseq, onums = odc
    r = tm_align(oca, ca, oseq, seq)
    corr, i, j = {}, 0, 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-":
            corr[onums[i]] = j
        i += x != "-"
        j += y != "-"
    on = oca @ np.asarray(r.u).T + np.asarray(r.t)
    pl = plddt_by_resnum(path)
    out = {}
    for pos in positions:
        if pos in corr:
            jj = corr[pos]
            off = float(np.linalg.norm(on[onums.index(pos)] - ca[jj]))
            out[pos] = gate(seq[jj], odc_seq[pos - 1], off, pl.get(nums[jj]))
        else:
            out[pos] = "unresolved"
    return out, float(r.tm_norm_chain1)


def af_model(acc):
    for d in (P10_MODELS, MODELS):
        p = d / f"AF-{acc}.pdb"
        if p.exists():
            return p
    miss = MODELS / f"AF-{acc}.none"
    if miss.exists():
        return None
    try:
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60) as r:
            ent = json.loads(r.read().decode("utf8"))[0]
    except urllib.error.HTTPError as e:
        miss.write_text(f"HTTP {e.code}")
        return None
    st, _, body = get(ent["pdbUrl"])
    p = MODELS / f"AF-{acc}.pdb"
    p.write_text(body, encoding="utf8")
    return p


def group_of(lineage_names):
    names = set(lineage_names)
    return next((g for g in FOCUS if g in names), "other")


def frac(states, positions):
    c = Counter(states[p] for p in positions)
    k, l = c["kept"], c["lost"]
    return (k / (k + l) if k + l else None), k, l, c["unresolved"]


# ---------------------------------------------------------------- 4. statistics
def sign_test_greater(k, n):
    """Exact one-sided binomial P(X >= k | n, 0.5)."""
    return sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n


def perm_strata(rows, key, g1, g2, stratum):
    """Mean(key | g1) - mean(key | g2); labels shuffled within stratum (p5 control)."""
    rows = [r for r in rows if r[key] is not None and r["group"] in (g1, g2)]
    strata = defaultdict(list)
    for r in rows:
        strata[r[stratum]].append((r[key], r["group"]))
    both = sum(1 for s in strata.values() if {g for _, g in s} == {g1, g2})

    def stat(ss):
        a = [v for s in ss for v, g in s if g == g1]
        b = [v for s in ss for v, g in s if g == g2]
        return np.mean(a) - np.mean(b)
    groups = list(strata.values())
    n1 = sum(1 for s in groups for _, g in s if g == g1)
    n2 = sum(1 for s in groups for _, g in s if g == g2)
    if n1 < 3 or n2 < 3:
        return None, None, both, n1, n2
    obs = stat(groups)
    if both == 0:
        return obs, None, both, n1, n2
    rng, hits = random.Random(1), 0
    for _ in range(N_PERM):
        perm = []
        for s in groups:
            labels = [g for _, g in s]
            rng.shuffle(labels)
            perm.append([(v, g) for (v, _), g in zip(s, labels)])
        hits += abs(stat(perm)) >= abs(obs) - 1e-12
    return obs, (hits + 1) / (N_PERM + 1), both, n1, n2


def f3(x):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.3g}"


# ---------------------------------------------------------------- main
def main():
    from p6_decarboxylase_partner import references
    refs = references()
    odc_seq = refs["P11926"]
    ado_seq = human_sequence()
    rep = ["# Track B mechanism: interface vs catalytic conservation (p12, generated)", "",
           "Generated by `scripts/p12_interface_vs_catalytic.py`; definitions in its docstring.", ""]

    # 1. interfaces
    odc_if, odc_bsa, odc_mm = interface("1D7K", "A", "B", odc_seq)
    ado_if, ado_bsa, ado_mm = interface("1MSV", "A", "B", ado_seq)
    cat_ado = [s[1] for s in marc_sites()]
    alt = {s[1]: s[4] for s in marc_sites()}
    with open(RESULTS / "p12_interface_residues.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, ["enzyme", "structure", "pos", "uniprot_aa", "chains_present", "bsa_mean_A2",
                                "min_dist_to_partner_A", "contact_4.0", "contact_4.5", "contact_5.0", "catalytic"],
                           delimiter="\t")
        w.writeheader()
        for enz, pid, rows, cat in (("ODC", "1D7K", odc_if, ODC_CAT), ("AdoMetDC", "1MSV", ado_if, cat_ado)):
            for r in rows:
                w.writerow({"enzyme": enz, "structure": pid, **r, "catalytic": r["pos"] in cat})
    I_odc = [r["pos"] for r in odc_if if r["contact_4.5"]]
    I_ado = [r["pos"] for r in ado_if if r["contact_4.5"]]
    I_odc_m = [p for p in I_odc if p <= p7.ODC_ORDERED]
    rep += ["## 1. Interfaces derived from crystal structures", "",
            "| enzyme | structure | buried surface, both chains (A^2) | interface residues (4.0 / 4.5 / 5.0 A) | "
            "residues with any buried surface | in the catalytic set too | numbering mismatches vs UniProt |",
            "|---|---|---|---|---|---|---|"]
    for enz, pid, rows, bsa, cat, mm in (("ODC", "1D7K", odc_if, odc_bsa, ODC_CAT, odc_mm),
                                         ("AdoMetDC", "1MSV", ado_if, ado_bsa, cat_ado, ado_mm)):
        I = [r["pos"] for r in rows if r["contact_4.5"]]
        rep.append(f"| {enz} | {pid} A+B | {bsa:.0f} | {sum(r['contact_4.0'] for r in rows)} / {len(I)} / "
                   f"{sum(r['contact_5.0'] for r in rows)} | {sum(r['bsa_mean_A2'] > 1 for r in rows)} | "
                   f"{', '.join(f'{(odc_seq if enz == 'ODC' else ado_seq)[p - 1]}{p}' for p in cat if p in I) or 'none'} | "
                   f"{len(mm)}{' (' + ', '.join(mm[:8]) + ')' if mm else ''} |")
    rep += ["", f"ODC interface positions beyond the ordered human-ODC model (> {p7.ODC_ORDERED}), not mappable: "
            f"{', '.join(str(p) for p in I_odc if p > p7.ODC_ORDERED) or 'none'}.", ""]
    I_odc_only = [p for p in I_odc_m if p not in ODC_CAT]
    I_ado_only = [p for p in I_ado if p not in cat_ado]

    # 2-3. ODC side
    integ = {r["accession"]: r for r in p7.read("p7_partner_integrity.tsv")}
    p10f = {r["accession"]: r for r in p7.read("p10_all_fusions.tsv")}
    p10s = {r["uniprot_accessions"].split(",")[0].split(".")[0]: r for r in p7.read("p10_fusion_search.tsv") if r["uniprot_accessions"]}
    odc_s = load_structure(str(p7.af_model("P11926")))
    odc = p7.sub_structure(odc_s, 1, p7.ODC_ORDERED)
    gated = defaultdict(dict)
    for r in p7.read("p9_odc_sites_gated.tsv"):
        gated[r["accession"]][int(r["human_odc_pos"])] = r["state"]
    targets = []   # (acc, kind, path, lo, hi, group, status, organism)
    for a, r in integ.items():
        if a in p10f and p10f[a]["set"] == "our 23":
            targets.append((a, "fusion (our 23)", RUN / "models" / f"{a}.pdb", int(r["partner_from"]), 10 ** 6,
                            p10f[a]["group"], p10f[a]["odc_status"], r["organism"]))
    if not SMOKE:
        for a, r in p10f.items():
            if r["set"] != "new" or not r["odc_status"] or not (P10_MODELS / f"AF-{a}.pdb").exists():
                continue
            s = p10s.get(a)
            if not s:
                continue
            lo = max(1, int(s["plp_start"]) - 25)
            hi = 10 ** 6 if s["adometdc_then_plp"] == "True" else int(s["plp_end"]) + 25
            targets.append((a, "fusion (new, p10)", P10_MODELS / f"AF-{a}.pdb", lo, hi, r["group"], r["odc_status"],
                            r["organism"]))
        base = p10.uniprot_tsv_like("reviewed:true AND ec:4.1.1.17 AND xref:pfam-PF02784", "accession,organism_name,lineage")
        for r in base:
            path = af_model(r["Entry"])
            if path:
                targets.append((r["Entry"], "standalone ODC (reviewed)", path, 1, 10 ** 6,
                                group_of([x.split(" (")[0] for x in r["Taxonomic lineage"].split(", ")]), "standalone",
                                r["Organism"]))
        for r in p7.read("p10_species_enzymes.tsv"):
            if r["class"].startswith("ODC family") and r["structure"]:
                acc = r["structure"].split(":")[0]
                path = af_model(acc)
                if path:
                    targets.append((acc, "standalone ODC-family (red alga / Percolomonas, p10)", path, 1, 10 ** 6,
                                    "Rhodophyta" if r["species"] != "Percolomonas cosmopolitus" else "Heterolobosea",
                                    "standalone", r["species"]))
    odc_rows, mismatch = [], []
    seen = set()
    for a, kind, path, lo, hi, grp, status, org in targets:
        if a in seen:
            continue
        seen.add(a)
        st, tm = map_odc(odc, odc_seq, path, lo, hi, sorted(set(ODC_CAT) | set(I_odc_m)))
        if a in gated:
            mismatch += [f"{a}:{p}" for p in ODC_CAT if gated[a][p] != st[p]]
        fc, kc, lc, uc = frac(st, ODC_CAT)
        fi, ki, li, ui = frac(st, I_odc_only)
        p_f = fisher(ki, li, kc, lc) if (ki + li) and (kc + lc) else None
        if kind.startswith("standalone"):
            g = "standalone ODC"
        elif status == "degraded" and grp in ("Rhodophyta", "Heterolobosea"):
            g = "degraded (red + Percolomonas)"
        elif status == "intact":
            g = "intact fusion"
        else:
            g = f"other fusion ({status})"
        odc_rows.append({"accession": a, "organism": org, "kind": kind, "stratum": grp, "odc_status": status,
                         "group": g, "TM_to_human_ODC": round(tm, 3),
                         "catalytic_kept": kc, "catalytic_lost": lc, "catalytic_unresolved": uc, "catalytic_frac": fc,
                         "interface_kept": ki, "interface_lost": li, "interface_unresolved": ui, "interface_frac": fi,
                         "d_interface_minus_catalytic": None if fi is None or fc is None else fi - fc,
                         "fisher_p_interface_vs_catalytic": p_f})
    write("p12_odc_per_protein.tsv", odc_rows)
    rep += ["## 2-3. ODC half", "",
            f"Catalytic set: the 13 residues. Interface set: {len(I_odc_m)} mappable 1D7K contact residues, "
            f"{len(I_odc_only)} after removing those also catalytic. Check against p9: "
            f"{len(mismatch)} catalytic calls differ from `p9_odc_sites_gated.tsv`"
            + (f" ({', '.join(mismatch[:10])})" if mismatch else "") + ".", ""]
    rep += summary_table(odc_rows)
    tests = []
    deg = "degraded (red + Percolomonas)"
    rep += ["", "Per degraded protein:", "", "| accession | organism | catalytic kept/lost/unres | interface kept/lost/unres | "
            "catalytic frac | interface frac | Fisher p |", "|---|---|---|---|---|---|---|"]
    for r in odc_rows:
        if r["group"] == deg:
            rep.append(f"| {r['accession']} | *{r['organism']}* | {r['catalytic_kept']}/{r['catalytic_lost']}/{r['catalytic_unresolved']} | "
                       f"{r['interface_kept']}/{r['interface_lost']}/{r['interface_unresolved']} | {f3(r['catalytic_frac'])} | "
                       f"{f3(r['interface_frac'])} | {f3(r['fisher_p_interface_vs_catalytic'])} |")
    rep += ["", "Tests:", ""]
    rep += group_tests(odc_rows, "ODC", deg, ["intact fusion", "standalone ODC"], "stratum", tests)

    # 4. AdoMetDC side
    human = human_sequence()
    ref = load_structure(str(REF_PDB), chain=REF_CHAIN)
    ref_idx = {n: k for k, n in enumerate(ref.resnums)}
    I_ado_map = [p for p in I_ado if p in ref_idx]
    p8.POSITIONS = sorted(set(cat_ado) | set(I_ado_map))
    p5 = {r["accession"]: r for r in p7.read("p5_loss_vs_fusion_proteins.tsv")}
    status = {r["accession"]: r for r in odc_rows}
    ado_targets = [(a, None, "fusion (our 23)") for a in integ if a in p10f and p10f[a]["set"] == "our 23"]
    if not SMOKE:
        ado_targets += [(a, P10_MODELS / f"AF-{a}.pdb", "fusion (new, p10)") for a, r in p10f.items()
                        if r["set"] == "new" and r["adometdc_half"] and (P10_MODELS / f"AF-{a}.pdb").exists()
                        and r["adometdc_then_plp"] == "True"]
        ado_targets += [(a, None, "standalone AdoMetDC") for a, r in p5.items() if r["class"] == "standalone"]
    ado_rows = []
    for k, (a, path, kind) in enumerate(ado_targets, 1):
        st_, _ = p8.map_sites(f"p10_{a}" if path else a, ref, ref_idx, human, mpath=path)
        st = {p: ("kept" if v["state"] == "lost" and v["cand_aa"] in alt.get(p, ()) else v["state"]) for p, v in st_.items()}
        fc, kc, lc, uc = frac(st, cat_ado)
        fi, ki, li, ui = frac(st, I_ado_only)
        p_f = fisher(ki, li, kc, lc) if (ki + li) and (kc + lc) else None
        odcst = status.get(a, {}).get("odc_status", "")
        g = "standalone AdoMetDC" if kind.startswith("standalone") else (
            "fusion, ODC degraded" if odcst == "degraded" else "fusion, ODC intact" if odcst == "intact" else "fusion, other")
        ado_rows.append({"accession": a, "organism": (p5.get(a) or p10f.get(a) or {}).get("organism", ""), "kind": kind,
                         "stratum": (p5[a]["tax_phylum"] or "(none)") if a in p5 else (p10f.get(a, {}).get("group", "")),
                         "odc_status": odcst, "group": g,
                         "catalytic_kept": kc, "catalytic_lost": lc, "catalytic_unresolved": uc, "catalytic_frac": fc,
                         "interface_kept": ki, "interface_lost": li, "interface_unresolved": ui, "interface_frac": fi,
                         "d_interface_minus_catalytic": None if fi is None or fc is None else fi - fc,
                         "fisher_p_interface_vs_catalytic": p_f})
        if k % 100 == 0:
            print(f"  AdoMetDC {k}/{len(ado_targets)}", flush=True)
    write("p12_adometdc_per_protein.tsv", ado_rows)
    rep += ["", "## 4. AdoMetDC half", "",
            f"Catalytic set: Marc's five sites (S229T and E247Q count as kept, as in p9). Interface set: "
            f"{len(I_ado_map)} 1MSV contact residues with coordinates in chain A, {len(I_ado_only)} after removing "
            "catalytic ones.", ""]
    rep += summary_table(ado_rows)
    rep += ["", "Tests:", ""]
    for r in ado_rows:
        r["group2"] = "standalone AdoMetDC" if r["kind"] == "standalone AdoMetDC" else "fusion (all)"
    ado2 = [dict(r, group=r["group2"]) for r in ado_rows]
    rep += group_tests(ado2, "AdoMetDC", "fusion (all)", ["standalone AdoMetDC"], "stratum", tests, key="interface_frac")
    rep += group_tests(ado_rows, "AdoMetDC", "fusion, ODC degraded", ["fusion, ODC intact", "standalone AdoMetDC"],
                       "stratum", tests)
    write("p12_tests.tsv", tests)
    (RESULTS / "p12_report.md").write_text("\n".join(rep) + "\n", encoding="utf8")
    print("\n".join(rep))


def write(name, rows):
    if not rows:
        return
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(RESULTS / name, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f3(v) if isinstance(v, float) else ("" if v is None else v)) for k, v in r.items()})


def summary_table(rows):
    out = ["| group | n proteins | catalytic frac, median (range) | interface frac, median (range) | "
           "d = interface - catalytic, median | proteins with d > 0 / d < 0 / d = 0 |", "|---|---|---|---|---|---|"]
    for g in sorted({r["group"] for r in rows}):
        rs = [r for r in rows if r["group"] == g]

        def med(k):
            v = [r[k] for r in rs if r[k] is not None]
            return f"{np.median(v):.2f} ({min(v):.2f}-{max(v):.2f})" if v else "NA"
        d = [r["d_interface_minus_catalytic"] for r in rs if r["d_interface_minus_catalytic"] is not None]
        out.append(f"| {g} | {len(rs)} | {med('catalytic_frac')} | {med('interface_frac')} | "
                   f"{f'{np.median(d):.2f}' if d else 'NA'} | {sum(x > 0 for x in d)} / {sum(x < 0 for x in d)} / {sum(x == 0 for x in d)} |")
    return out


def group_tests(rows, enz, g1, others, stratum, tests, key="d_interface_minus_catalytic"):
    out = []
    d = [r["d_interface_minus_catalytic"] for r in rows if r["group"] == g1 and r["d_interface_minus_catalytic"] is not None]
    pos, neg = sum(x > 0 for x in d), sum(x < 0 for x in d)
    if key == "d_interface_minus_catalytic":
        if pos + neg >= 3:
            p = sign_test_greater(pos, pos + neg)
            out.append(f"- {enz}, {g1}: interface > catalytic in {pos} of {pos + neg} proteins with a difference "
                       f"(exact sign test, one-sided p = {p:.3g}).")
            tests.append({"enzyme": enz, "test": "sign test d > 0", "group": g1, "vs": "", "n1": pos + neg, "n2": "",
                          "statistic": pos, "p": f"{p:.3g}", "strata_with_both": ""})
        else:
            out.append(f"- {enz}, {g1}: {pos + neg} proteins with a difference; too few for a sign test.")
    for g2 in others:
        obs, p, both, n1, n2 = perm_strata(rows, key, g1, g2, stratum)
        if obs is None:
            out.append(f"- {enz}, {g1} vs {g2} on {key}: n = {n1} / {n2}; too few to test.")
            continue
        out.append(f"- {enz}, {g1} vs {g2} on {key}: difference in means {obs:+.3f} (n = {n1} / {n2}); "
                   + (f"within-{stratum} permutation p = {p:.3g} ({both} strata hold both groups)." if p is not None else
                      "no stratum holds both groups, so no lineage-controlled p-value is possible."))
        tests.append({"enzyme": enz, "test": f"mean {key}, within-{stratum} permutation", "group": g1, "vs": g2,
                      "n1": n1, "n2": n2, "statistic": f"{obs:.3f}", "p": "" if p is None else f"{p:.3g}",
                      "strata_with_both": both})
    return out


class p10:   # minimal UniProt TSV fetch, cached in data/p12 (does not touch p10's cache)
    @staticmethod
    def uniprot_tsv_like(query, fields):
        c = P12 / ("uniprot_" + re.sub(r"\W+", "_", query)[:100] + ".tsv")
        if not c.exists():
            st, _, body = get("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
                {"query": query, "format": "tsv", "fields": fields, "size": 500}))
            c.write_text(body, encoding="utf8")
        return list(csv.DictReader(c.read_text(encoding="utf8").splitlines(), delimiter="\t"))


if __name__ == "__main__":
    sys.exit(main())
