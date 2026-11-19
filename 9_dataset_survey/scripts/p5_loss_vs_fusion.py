"""Track A: is loss of FOLDTRACE sites less common in fused AdoMetDCs, once relatedness is accounted for?

Inputs (all already in results/): master_candidates.tsv (FOLDTRACE site calls, fold gate),
architecture_per_protein.tsv (Pfam architecture call). Taxonomic lineage for every protein comes from
the UniProt taxonomy API by taxid (ranked lineage), cached in data/taxonomy.json, so active and
deleted accessions get the same kind of lineage.

Definitions (fixed here, before looking at the result):
  set         fold-passers, excluding A0A8C5BDJ7 (sequence/model mismatch)
  standalone  architecture "AdoMetDC only" (one PF01536/PF02675 match, no other Pfam domain)
  fused       AdoMetDC Pfam plus at least one other Pfam domain (AdoMetDC first / last / middle)
  repeated    two or more AdoMetDC Pfam matches and nothing else: kept separate (not fused, not
              clearly standalone); sensitivity analysis counts them as standalone
  loss        at least one of the five sites called "lost" by FOLDTRACE
  denominators  A = all proteins in the class; B = only proteins with all five sites resolved
Relatedness controls
  - number of distinct taxa (genus ... phylum) among proteins with a loss
  - collapsed test: one unit per (taxon, class); a unit has a loss if any member does; Fisher on units
  - within-lineage permutation: fusion labels shuffled only within each taxon (class, order),
    10,000 permutations; plus the Mantel-Haenszel odds ratio across those strata
Also: pLDDT (whole model, AdoMetDC domain, the five sites), taxonomic makeup, each site separately,
and every result repeated for the InterPro-only set (annotation source = interpro).

Writes results/p5_loss_vs_fusion_report.txt and results/p5_loss_vs_fusion_proteins.tsv
Run (about a minute the first time, for taxonomy and pLDDT):  py -3.12 p5_loss_vs_fusion.py
"""
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

from p2_fetch_annotations import DATA, RESULTS, get, save_json

RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
SITES = ["Ser68", "Ser229", "His243", "Cys82", "Glu247"]
EXCLUDE = {"A0A8C5BDJ7"}
RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]
ADOMETDC = {"PF01536", "PF02675"}
LEADER = "PF08132"
N_PERM = 10000
random.seed(1)


def read(name):
    with open(RESULTS / name, encoding="utf8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ statistics
def lf(n):
    return math.lgamma(n + 1)


def fisher(a, b, c, d):
    """Two-sided Fisher exact test for [[a, b], [c, d]] (sum of tables no more likely than observed)."""
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d

    def p(x):
        return math.exp(lf(r1) + lf(r2) + lf(c1) + lf(n - c1) - lf(n) - lf(x) - lf(r1 - x) - lf(c1 - x) - lf(r2 - c1 + x))
    obs = p(a)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(p(x) for x in range(lo, hi + 1) if p(x) <= obs * (1 + 1e-9)))


def mann_whitney(x, y):
    """Two-sided Mann-Whitney U, normal approximation with tie correction. Returns (U, p)."""
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks, i = [0.0] * len(allv), 0
    ties = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[k] = (i + j) / 2 + 1
        t = j - i + 1
        ties += t ** 3 - t
        i = j + 1
    n1, n2 = len(x), len(y)
    r1 = sum(r for r, (_, g) in zip(ranks, allv) if g == 0)
    u = r1 - n1 * (n1 + 1) / 2
    n = n1 + n2
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1))))
    z = (u - n1 * n2 / 2) / sd if sd else 0
    return u, math.erfc(abs(z) / math.sqrt(2))


def mantel_haenszel(strata):
    """strata: list of (a, b, c, d) with a = standalone loss, b = standalone no loss, c = fused loss, d = fused no loss."""
    num = sum(a * d / (a + b + c + d) for a, b, c, d in strata if a + b and c + d)
    den = sum(b * c / (a + b + c + d) for a, b, c, d in strata if a + b and c + d)
    return num / den if den else float("inf") if num else float("nan")


def rate(k, n):
    return f"{k}/{n} ({100 * k / n:.1f}%)" if n else "0/0"


# ------------------------------------------------------------------ data
def taxonomy(taxids):
    p = DATA / "taxonomy.json"
    cache = json.loads(p.read_text(encoding="utf8")) if p.exists() else {}
    todo = sorted({t for t in taxids if t and t not in cache})
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        q = "%20OR%20".join(f"tax_id:{t}" for t in chunk)
        st, _, body = get(f"https://rest.uniprot.org/taxonomy/search?query={q}&format=json"
                          f"&fields=id,scientific_name,rank,lineage&size=500")
        for r in json.loads(body).get("results", []):
            lin = {x.get("rank"): x.get("scientificName", "") for x in r.get("lineage", [])}
            lin[r.get("rank")] = r.get("scientificName", "")
            cache[str(r["taxonId"])] = {"name": r.get("scientificName", ""), "ranks": lin,
                                        "path": [x.get("scientificName", "") for x in r.get("lineage", [])],
                                        "inactive": r.get("inactiveReason")}
        for t in chunk:
            cache.setdefault(t, {"name": "", "ranks": {}, "path": []})
        save_json(p, cache)
    return cache


def ranked(tax, rank):
    """Name at a rank; if the taxon has no such rank, fall back to the nearest named higher rank."""
    r = tax["ranks"]
    if r.get(rank):
        return r[rank]
    for up in reversed(RANKS[:RANKS.index(rank)]):
        if r.get(up):
            return f"{r[up]} (no {rank})"
    return "(unknown)"


def domain_plddt(acc, spans):
    pl = {}
    with open(RUN / "models" / f"{acc}.pdb", encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("ATOM") and ln[12:16] == " CA ":
                pl[int(ln[22:26])] = float(ln[60:66])
    vals = [v for n, v in pl.items() if any(s <= n <= e for s, e in spans)]
    return sum(vals) / len(vals) if vals else None


def load():
    master = {r["accession"]: r for r in read("master_candidates.tsv")}
    arch = {r["accession"]: r for r in read("architecture_per_protein.tsv")}
    tax = taxonomy({m["taxid"] for m in master.values()})
    rows = []
    for a, m in master.items():
        if m["fold_ok"] != "true" or a in EXCLUDE:
            continue
        ar = arch[a]
        order = ar["adometdc_order"]
        cls = {"AdoMetDC only": "standalone", "AdoMetDC repeated": "repeated"}.get(
            order, "fused" if order.startswith("AdoMetDC") else "other")
        partners = {x.split()[0] for x in ar["partner_pfam"].split("; ") if x}
        states = {s: m[f"{s}_state"] for s in SITES}
        spans = [(int(p.split()[1].split("-")[0]), int(p.split()[1].split("-")[1]))
                 for p in ar["domains_with_positions"].split("; ") if p and p.split()[0] in ADOMETDC]
        site_pl = [float(m[f"{s}_plddt"]) for s in SITES if m[f"{s}_plddt"]]
        t = tax.get(m["taxid"], {"ranks": {}, "path": []})
        row = {"accession": a, "organism": m["organism"], "taxid": m["taxid"], "source": ar["source"],
               "architecture": ar["arch_pfam_ids"], "class": cls,
               "leader_peptide_only": cls == "fused" and partners == {LEADER},
               "loss": any(v == "lost" for v in states.values()),
               "all_resolved": all(v != "unresolved" for v in states.values()),
               "lost_sites": ",".join(s for s in SITES if states[s] == "lost"),
               "plddt_model": float(m["plddt_mean"]),
               "plddt_adometdc_domain": domain_plddt(a, spans) if spans else None,
               "plddt_sites": sum(site_pl) / len(site_pl) if site_pl else None,
               "lost_site_plddt": ",".join(m[f"{s}_plddt"] for s in SITES if states[s] == "lost"),
               **{f"{s}_state": states[s] for s in SITES},
               **{f"tax_{r}": ranked(t, r) for r in RANKS}}
        rows.append(row)
    return rows


# ------------------------------------------------------------------ analysis
def compare(rows, label, out):
    S = [r for r in rows if r["class"] == "standalone"]
    F = [r for r in rows if r["class"] == "fused"]
    out.append(f"\n{'=' * 100}\n{label}\n{'=' * 100}")
    out.append(f"classes: {dict(Counter(r['class'] for r in rows))}")

    def test(s, f, name):
        a, c = sum(r["loss"] for r in s), sum(r["loss"] for r in f)
        p = fisher(a, len(s) - a, c, len(f) - c) if s and f else float("nan")
        out.append(f"  {name:<62} standalone {rate(a, len(s)):<16} fused {rate(c, len(f)):<14} Fisher p = {p:.3g}")
        return a, c, p

    out.append("\n[1] Loss of >= 1 site, standalone vs fused (architecture-based)")
    test(S, F, "A: all proteins")
    test([r for r in S if r["all_resolved"]], [r for r in F if r["all_resolved"]], "B: all five sites resolved")
    out.append("  sensitivity:")
    test(S + [r for r in rows if r["class"] == "repeated"], F, "repeated PF01536 counted as standalone")
    test(S, [r for r in F if not r["leader_peptide_only"]], "fused excluding leader-peptide-only (PF08132)")
    test(S, [r for r in F if r["leader_peptide_only"]], "leader-peptide-only fusions vs standalone")
    out.append(f"  repeated PF01536 class on its own: {rate(sum(r['loss'] for r in rows if r['class'] == 'repeated'), sum(r['class'] == 'repeated' for r in rows))}")

    out.append("\n[5] Each site separately (resolved calls only: lost vs retained+altered)")
    for s in SITES:
        ss = [r for r in S if r[f"{s}_state"] != "unresolved"]
        ff = [r for r in F if r[f"{s}_state"] != "unresolved"]
        a = sum(r[f"{s}_state"] == "lost" for r in ss)
        c = sum(r[f"{s}_state"] == "lost" for r in ff)
        p = fisher(a, len(ss) - a, c, len(ff) - c) if ss and ff else float("nan")
        out.append(f"  {s:<8} standalone {rate(a, len(ss)):<16} fused {rate(c, len(ff)):<14} Fisher p = {p:.3g}")

    out.append("\n[2] Relatedness")
    L = [r for r in S + F if r["loss"]]
    for rk in ("genus", "family", "order", "class", "phylum"):
        out.append(f"  proteins with a loss ({len(L)}) fall in {len({r["tax_" + rk] for r in L})} distinct {rk}")
    out.append("  losses by order: " + "; ".join(f"{k}: {v}" for k, v in Counter(r["tax_order"] for r in L).most_common()))
    out.append("  proteins with a loss: " + "; ".join(f"{r['accession']} {r['organism']} [{r['class']}] lost {r['lost_sites']}"
                                                     f" (pLDDT {r['lost_site_plddt']})" for r in L))
    for rk in ("genus", "family", "order", "class"):
        units = defaultdict(list)
        for r in S + F:
            units[(r["tax_" + rk], r["class"])].append(r["loss"])
        su = [any(v) for (t, c), v in units.items() if c == "standalone"]
        fu = [any(v) for (t, c), v in units.items() if c == "fused"]
        p = fisher(sum(su), len(su) - sum(su), sum(fu), len(fu) - sum(fu)) if su and fu else float("nan")
        out.append(f"  collapsed to one unit per {rk:<7}: standalone {rate(sum(su), len(su)):<16} "
                   f"fused {rate(sum(fu), len(fu)):<14} Fisher p = {p:.3g}")
    for rk in ("phylum", "class", "order"):
        strata = defaultdict(list)
        for r in S + F:
            strata[r["tax_" + rk]].append(r)
        tabs, informative = [], []
        for t, rs in strata.items():
            a = sum(r["loss"] for r in rs if r["class"] == "standalone")
            b = sum(not r["loss"] for r in rs if r["class"] == "standalone")
            c = sum(r["loss"] for r in rs if r["class"] == "fused")
            d = sum(not r["loss"] for r in rs if r["class"] == "fused")
            tabs.append((a, b, c, d))
            if a + b and c + d:
                informative.append((t, a, b, c, d))

        def stat(groups):
            s_n = sum(1 for g in groups for r in g if r[1] == "standalone")
            f_n = sum(1 for g in groups for r in g if r[1] == "fused")
            s_l = sum(1 for g in groups for r in g if r[1] == "standalone" and r[0])
            f_l = sum(1 for g in groups for r in g if r[1] == "fused" and r[0])
            return s_l / s_n - f_l / f_n
        groups = [[(r["loss"], r["class"]) for r in rs] for rs in strata.values()]
        obs = stat(groups)
        hits = 0
        for _ in range(N_PERM):
            perm = []
            for g in groups:
                labels = [c for _, c in g]
                random.shuffle(labels)
                perm.append([(l, c) for (l, _), c in zip(g, labels)])
            hits += abs(stat(perm)) >= abs(obs) - 1e-12
        mh = mantel_haenszel(tabs)
        inf_loss = sum(a + c for _, a, b, c, d in informative)
        out.append(f"  within-{rk:<6} permutation: observed difference {100 * obs:.2f} points, p = {(hits + 1) / (N_PERM + 1):.3g}; "
                   f"MH odds ratio {mh:.3g}; {len(informative)} {rk} contain both classes "
                   f"({sum(a + b + c + d for _, a, b, c, d in informative)} proteins, {inf_loss} with a loss)")

    out.append("\n[3] pLDDT, standalone vs fused (Mann-Whitney, two-sided)")
    for key in ("plddt_model", "plddt_adometdc_domain", "plddt_sites"):
        x = [r[key] for r in S if r[key] is not None]
        y = [r[key] for r in F if r[key] is not None]
        _, p = mann_whitney(x, y)
        md = lambda v: sorted(v)[len(v) // 2]
        out.append(f"  {key:<24} standalone mean {sum(x) / len(x):.1f} median {md(x):.1f} (n {len(x)}); "
                   f"fused mean {sum(y) / len(y):.1f} median {md(y):.1f} (n {len(y)}); p = {p:.3g}")
    lost_pl = [float(v) for r in L for v in r["lost_site_plddt"].split(",") if v]
    out.append(f"  pLDDT at the lost sites themselves: {len(lost_pl)} calls, {sum(v < 70 for v in lost_pl)} below 70, "
               f"median {sorted(lost_pl)[len(lost_pl) // 2] if lost_pl else 'NA'}")
    hi = [r for r in S + F if r["plddt_sites"] is not None and r["plddt_sites"] >= 70]
    test([r for r in hi if r["class"] == "standalone"], [r for r in hi if r["class"] == "fused"],
         "loss test restricted to mean site pLDDT >= 70")

    out.append("\n[4] Taxonomic makeup")
    for rk in ("domain", "kingdom", "phylum"):
        out.append(f"  by {rk}:")
        cs, cf = Counter(r["tax_" + rk] for r in S), Counter(r["tax_" + rk] for r in F)
        for t in sorted(set(cs) | set(cf), key=lambda t: -(cs[t] + cf[t])):
            ls = sum(r["loss"] for r in S if r["tax_" + rk] == t)
            lfu = sum(r["loss"] for r in F if r["tax_" + rk] == t)
            out.append(f"    {t:<45} standalone {cs[t]:>4} (loss {ls:>2})   fused {cf[t]:>3} (loss {lfu})")


def main():
    rows = load()
    out = [f"Track A: site loss vs fusion. {len(rows)} fold-passers (A0A8C5BDJ7 excluded). "
           f"Taxonomy: UniProt taxonomy API, cached in data/taxonomy.json. Permutations: {N_PERM}, seed 1."]
    compare(rows, "ALL SOURCES (interpro + replacement_record + uniparc_features)", out)
    compare([r for r in rows if r["source"] == "interpro"], "INTERPRO-ONLY SET", out)
    txt = "\n".join(out) + "\n"
    (RESULTS / "p5_loss_vs_fusion_report.txt").write_text(txt, encoding="utf8")
    with open(RESULTS / "p5_loss_vs_fusion_proteins.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(txt)


if __name__ == "__main__":
    main()
