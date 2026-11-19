"""Phase 4a: where do the 773 fold-passing proteins differ from human, and how confident is each difference?

No hypothesis is prioritised: every human position is treated the same way. The five prespecified
FOLDTRACE sites and the residues within 6 A of putrescine in 1MSV are only LABELLED, for orientation.

pLDDT RULE (applies everywhere): a difference is "high confidence" only if the pLDDT of that region in
that model is >= 70 (residue pLDDT for substitutions/displacements, mean over inserted residues for
insertions/replacements, mean of the two flanking residues for deletions). Low-confidence differences
are counted and reported in their own columns but never used in rankings, groups or "unique" lists.
Events at residues 22-27 (no coordinates in 1MSV) are also kept out of headline results.

Needs p3_residue_grid.py to have been run. Uses results/annotation_uniprot.tsv (Phase 2) for
taxonomic lineage if it exists; otherwise only the organism name is available.

Writes (results/p4_*):
  p4_positions_ranked.tsv        1. every human position: conservation, how many differ, mean pLDDT
  p4_position_differences.tsv    1. which proteins differ at which position (long list)
  p4_regions_ranked.tsv          4a. 9-residue windows that vary most
  p4_indel_hotspots.tsv          2. insertions/deletions grouped by where they sit in the human frame
  p4_proteins_ranked.tsv         4b. proteins that differ most, with the confidence of those regions
  p4_protein_change_profiles.tsv     one line per protein listing every high-confidence change
  p4_unique_changes.tsv          changes seen in only 1-2 proteins at positions conserved in >=90%
  p4_co_changing_groups.tsv      5. groups of proteins that changed at the same positions together
  p4_co_changing_members.tsv     5. members of each group
figures/p4_change_map.png        tracks along human positions
figures/p4_heatmap.png           proteins x human positions, coloured by call
Run (about a minute; figures need matplotlib: py -3.12 -m pip install matplotlib):
    py -3.12 p4_variation.py
"""
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser

BASE = Path(__file__).resolve().parent.parent
RES, FIG = BASE / "results", BASE / "figures"
RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
PLDDT_OK = 70.0
REF_GAP = range(22, 28)             # 1MSV has no coordinates here
SITES = {68: "Ser68", 82: "Cys82", 229: "Ser229", 243: "His243", 247: "Glu247"}
CONSERVED = 0.90                    # "conserved position" for the unique-change list
GROUP_MIN, GROUP_JACCARD, CORE_FRAC = 3, 0.7, 0.7
DIFF = ("substituted", "displaced", "deleted")
EXCLUDE_POSITIONS = {"A0A8C5BDJ7"}   # UniProt sequence 227 aa != AlphaFold model 245 aa


def read(name):
    with open(RES / name, encoding="utf8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, rows, fields=None):
    with open(RES / name, "w", newline="", encoding="utf8") as fh:
        if not rows:
            fh.write("(none)\n")
            return
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fnum(x):
    return float(x) if x not in ("", None) else None


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 1) if xs else ""


def pocket_positions():
    s = PDBParser(QUIET=True).get_structure("r", str(RUN / "1MSV_A.pdb"))
    atoms = [a for a in s.get_atoms() if a.get_parent().id[0] == " "]
    lig = [a for a in s.get_atoms() if a.get_parent().get_resname() == "PUT"]
    ns = NeighborSearch(atoms)
    return sorted({a.get_parent().id[1] for l in lig for a in ns.search(l.coord, 6.0)})


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    grid = [g for g in read("p3_residue_grid.tsv") if g["accession"] not in EXCLUDE_POSITIONS]
    events = [e for e in read("p3_indel_events.tsv") if e["accession"] not in EXCLUDE_POSITIONS]
    source = {}
    if (RES / "architecture_per_protein.tsv").exists():
        source = {r["accession"]: r["source"] for r in read("architecture_per_protein.tsv")}
    summ = {r["accession"]: r for r in read("p3_protein_summary.tsv")}
    master = {r["accession"]: r for r in read("master_candidates.tsv")}
    lineage = {}
    if (RES / "annotation_uniprot.tsv").exists():
        lineage = {r["accession"]: r["lineage"] for r in read("annotation_uniprot.tsv")}
    have_lineage = bool(lineage)

    def org(a):
        return master[a]["organism"]

    def group_of(a, depth=3):
        lin = lineage.get(a, "")
        return ";".join(lin.split(";")[:depth]) if lin else ""

    pocket = set(pocket_positions())

    def label(p):
        out = []
        if p in SITES:
            out.append(f"FOLDTRACE site {SITES[p]}")
        if p in pocket:
            out.append("within 6 A of putrescine (1MSV)")
        if p in REF_GAP:
            out.append("no coordinates in 1MSV")
        return "; ".join(out)

    prot = sorted({g["accession"] for g in grid})
    N = len(prot)
    positions = sorted({int(g["ref_pos"]) for g in grid})
    human = {}
    by_pos = defaultdict(list)
    by_prot = defaultdict(dict)
    for g in grid:
        p = int(g["ref_pos"])
        human[p] = g["human_aa"]
        g["pl"] = fnum(g["plddt"])
        g["hc"] = g["pl"] is not None and g["pl"] >= PLDDT_OK
        by_pos[p].append(g)
        by_prot[g["accession"]][p] = g

    # ---------- events, with high-confidence flag
    for e in events:
        q = fnum(e["plddt_mean"]) if e["cand_len"] not in ("0", "") else fnum(e["plddt_flank_mean"])
        e["hc"] = q is not None and q >= PLDDT_OK
        e["q"] = q
        e["headline"] = e["hc"] and e["at_reference_gap"] != "True" and e["event"] in ("insertion", "deletion", "replacement")
    ins_after = Counter()
    for e in events:
        if e["headline"] and e["event"] == "insertion":
            ins_after[int(e["after_ref_pos"])] += 1

    # ---------- 1. per-position conservation
    pos_rows, diff_rows = [], []
    for p in positions:
        rows = by_pos[p]
        c = Counter((g["call"], g["hc"]) for g in rows)
        inplace_hc = [g["cand_aa"] for g in rows if g["hc"] and g["call"] in ("identical", "substituted")]
        aa = Counter(inplace_hc)
        H = -sum(v / len(inplace_hc) * math.log2(v / len(inplace_hc)) for v in aa.values()) if inplace_hc else ""
        n_hc_called = sum(1 for g in rows if g["hc"] and g["call"] != "not_covered")
        n_diff_hc = sum(c[(d, True)] for d in DIFF)
        n_diff_lc = sum(c[(d, False)] for d in DIFF)
        pos_rows.append({
            "ref_pos": p, "human_aa": human[p], "label": label(p),
            "n_proteins": N, "n_not_covered": c[("not_covered", False)] + c[("not_covered", True)],
            "n_high_conf_calls": n_hc_called,
            "pct_identical_high_conf": round(100 * c[("identical", True)] / n_hc_called, 1) if n_hc_called else "",
            "n_differ_high_conf": n_diff_hc,
            "n_substituted_hc": c[("substituted", True)], "n_displaced_hc": c[("displaced", True)],
            "n_deleted_hc": c[("deleted", True)],
            "n_differ_low_conf_excluded": n_diff_lc,
            "n_insertions_after_hc": ins_after[p],
            "entropy_bits_hc": round(H, 3) if H != "" else "",
            "residues_seen_hc": ",".join(f"{k}:{v}" for k, v in aa.most_common(6)),
            "mean_plddt": mean([g["pl"] for g in rows if g["call"] != "deleted"]),
            "mean_plddt_where_different": mean([g["pl"] for g in rows if g["call"] in DIFF]),
        })
        for g in rows:
            if g["call"] in DIFF:
                diff_rows.append({"ref_pos": p, "human_aa": human[p], "accession": g["accession"],
                                  "organism": org(g["accession"]), "group": group_of(g["accession"]),
                                  "call": g["call"], "cand_aa": g["cand_aa"], "cand_resnum": g["cand_resnum"],
                                  "ca_offset": g["ca_offset"], "plddt": g["plddt"],
                                  "confidence": "high" if g["hc"] else "low", "label": label(p)})
    pos_by = {r["ref_pos"]: r for r in pos_rows}
    ranked_pos = sorted((r for r in pos_rows if r["ref_pos"] not in REF_GAP),
                        key=lambda r: (-r["n_differ_high_conf"], -(r["entropy_bits_hc"] or 0)))
    for i, r in enumerate(ranked_pos, 1):
        r["rank"] = i
    write("p4_positions_ranked.tsv", ranked_pos, ["rank"] + list(pos_rows[0]))
    write("p4_position_differences.tsv", diff_rows)

    # ---------- 4a. regions (9-residue windows, non-overlapping greedy)
    W = 9
    pos_hc_prots = defaultdict(set)
    for d in diff_rows:
        if d["confidence"] == "high":
            pos_hc_prots[d["ref_pos"]].add(d["accession"])
    ev_by_prot = defaultdict(list)
    for e in events:
        ev_by_prot[e["accession"]].append(e)
    win = []
    for k in range(len(positions) - W + 1):
        ps = positions[k:k + W]
        if ps[-1] - ps[0] != W - 1 or any(p in REF_GAP for p in ps):
            continue
        prots_hc = set().union(*(pos_hc_prots[p] for p in ps))
        prots_hc |= {e["accession"] for e in events if e["headline"] and
                     int(e["after_ref_pos"] or 0) >= ps[0] - 1 and int(e["before_ref_pos"] or 0) <= ps[-1] + 1}
        win.append({"start": ps[0], "end": ps[-1], "human_seq": "".join(human[p] for p in ps),
                    "n_proteins_differ_hc": len(prots_hc),
                    "mean_pct_identical_hc": round(sum(pos_by[p]["pct_identical_high_conf"] or 0 for p in ps) / W, 1),
                    "n_indel_events_hc": sum(1 for e in events if e["headline"] and ps[0] - 1 <= int(e["after_ref_pos"] or 0) <= ps[-1]),
                    "mean_plddt": round(sum(pos_by[p]["mean_plddt"] or 0 for p in ps) / W, 1),
                    "labels": "; ".join(sorted({f"{p}:{label(p)}" for p in ps if label(p)}))})
    win.sort(key=lambda r: (-r["n_proteins_differ_hc"], r["mean_pct_identical_hc"]))
    chosen, used = [], set()
    for r in win:
        span = set(range(r["start"], r["end"] + 1))
        if span & used:
            continue
        used |= span
        chosen.append(r)
    for i, r in enumerate(chosen, 1):
        r["rank"] = i
    write("p4_regions_ranked.tsv", chosen, ["rank"] + [k for k in chosen[0] if k != "rank"] if chosen else None)

    # ---------- 2. indel hotspots (same event type, anchors within 2 positions)
    hot = []
    for ev in ("insertion", "deletion", "replacement"):
        es = sorted((e for e in events if e["event"] == ev), key=lambda e: int(e["after_ref_pos"]))
        cl = []
        for e in es:
            if cl and int(e["after_ref_pos"]) - int(cl[-1][-1]["after_ref_pos"]) <= 2:
                cl[-1].append(e)
            else:
                cl.append([e])
        for c in cl:
            hc = [e for e in c if e["hc"]]
            lens = [int(e["cand_len"]) if ev != "deletion" else int(e["n_ref_missing"]) for e in c]
            hot.append({"event": ev,
                        "human_region": f"{min(int(e['after_ref_pos']) for e in c)}-{max(int(e['before_ref_pos']) for e in c)}",
                        "at_reference_gap": any(e["at_reference_gap"] == "True" for e in c),
                        "n_proteins_high_conf": len({e["accession"] for e in hc}),
                        "n_proteins_low_conf_excluded": len({e["accession"] for e in c} - {e["accession"] for e in hc}),
                        "length_min": min(lens), "length_median": sorted(lens)[len(lens) // 2], "length_max": max(lens),
                        "mean_plddt_high_conf": mean([e["q"] for e in hc]),
                        "mean_plddt_all": mean([e["q"] for e in c]),
                        "n_groups_high_conf": len({group_of(e["accession"]) for e in hc}) if have_lineage else "",
                        "examples_high_conf": "; ".join(f"{e['accession']} ({org(e['accession'])}, {e['cand_len'] if ev != 'deletion' else e['n_ref_missing']} res, pLDDT {e['q']})" for e in hc[:6])})
    hot.sort(key=lambda r: (r["at_reference_gap"], -r["n_proteins_high_conf"]))
    write("p4_indel_hotspots.tsv", hot)

    # ---------- 4b. proteins ranked + change profiles
    prow, profiles = [], []
    for a in prot:
        cells = by_prot[a]
        hc_d = sorted(p for p, g in cells.items() if g["call"] in DIFF and g["hc"] and p not in REF_GAP)
        lc_d = sorted(p for p, g in cells.items() if g["call"] in DIFF and not g["hc"])
        evs = ev_by_prot[a]
        hc_e = [e for e in evs if e["headline"]]
        lc_e = [e for e in evs if e["event"] in ("insertion", "deletion", "replacement") and not e["headline"]]
        # merge differing positions into regions (allow gaps of <= 2)
        regs = []
        for p in hc_d:
            if regs and p - regs[-1][-1] <= 3:
                regs[-1].append(p)
            else:
                regs.append([p])
        reg_txt = [f"{r[0]}-{r[-1]} ({len(r)} diff, pLDDT {mean([cells[p]['pl'] for p in r])})"
                   for r in sorted(regs, key=len, reverse=True)[:8]]
        s = summ[a]
        changes = []
        for p in hc_d:
            g = cells[p]
            changes.append(f"{human[p]}{p}{g['cand_aa'] if g['call'] == 'substituted' else ('-' if g['call'] == 'deleted' else g['cand_aa'].lower() + '*')}")
        for e in hc_e:
            if e["event"] == "insertion":
                changes.append(f"ins{e['after_ref_pos']}^{e['before_ref_pos']}+{e['cand_len']}")
            elif e["event"] == "deletion":
                changes.append(f"del{e['ref_positions_missing']}")
            else:
                changes.append(f"repl{e['ref_positions_missing']}:{e['cand_len']}res")
        n_hc = len(hc_d) + sum(int(e["cand_len"]) if e["event"] != "deletion" else 0 for e in hc_e)
        prow.append({"accession": a, "organism": org(a), "group": group_of(a), "annotation_source": source.get(a, ""),
                     "length": s["cand_length"],
                     "tm_norm_ref": s["tm_norm_ref"], "seq_identity_in_place": s["seq_identity_in_place"],
                     "n_diff_positions_high_conf": len(hc_d), "n_indel_events_high_conf": len(hc_e),
                     "inserted_residues_high_conf": sum(int(e["cand_len"]) for e in hc_e if e["event"] != "deletion"),
                     "difference_score_high_conf": n_hc,
                     "n_diff_positions_low_conf_excluded": len(lc_d), "n_indel_events_low_conf_excluded": len(lc_e),
                     "mean_plddt_at_hc_differences": mean([cells[p]["pl"] for p in hc_d]),
                     "plddt_mean_aligned": s["plddt_mean_aligned"],
                     "main_difference_regions": "; ".join(reg_txt),
                     "N_ext_len": s["N_ext_len"], "C_ext_len": s["C_ext_len"],
                     "fragment_in_name": master[a]["fragment_in_name"]})
        profiles.append({"accession": a, "organism": org(a), "group": group_of(a),
                         "n_high_conf_changes": len(changes), "changes_high_conf": " ".join(changes)})
    prow.sort(key=lambda r: -r["difference_score_high_conf"])
    for i, r in enumerate(prow, 1):
        r["rank"] = i
    write("p4_proteins_ranked.tsv", prow, ["rank"] + [k for k in prow[0] if k != "rank"])
    write("p4_protein_change_profiles.tsv", profiles)

    # ---------- unique changes at conserved positions
    uniq = []
    sub_counts = Counter((d["ref_pos"], d["call"], d["cand_aa"]) for d in diff_rows if d["confidence"] == "high")
    for d in diff_rows:
        if d["confidence"] != "high" or d["ref_pos"] in REF_GAP:
            continue
        cons = pos_by[d["ref_pos"]]["pct_identical_high_conf"]
        n = sub_counts[(d["ref_pos"], d["call"], d["cand_aa"])]
        if cons != "" and cons >= 100 * CONSERVED and n <= 2:
            uniq.append({**d, "pct_identical_at_position": cons, "n_proteins_with_this_change": n})
    for h in hot:
        if h["n_proteins_high_conf"] in (1, 2) and not h["at_reference_gap"]:
            uniq.append({"ref_pos": h["human_region"], "call": h["event"], "accession": h["examples_high_conf"],
                         "n_proteins_with_this_change": h["n_proteins_high_conf"]})
    uniq.sort(key=lambda r: -(r.get("pct_identical_at_position") or 0))
    write("p4_unique_changes.tsv", uniq, ["ref_pos", "human_aa", "call", "cand_aa", "accession", "organism", "group",
                                          "plddt", "pct_identical_at_position", "n_proteins_with_this_change", "label",
                                          "cand_resnum", "ca_offset"])

    # ---------- 5. groups that changed together
    feats = defaultdict(set)
    for d in diff_rows:
        if d["confidence"] == "high" and d["ref_pos"] not in REF_GAP:
            feats[(d["ref_pos"], f"{d['human_aa']}{d['ref_pos']}{d['cand_aa'] or '-'}" if d["call"] != "displaced" else f"{d['human_aa']}{d['ref_pos']}displaced")].add(d["accession"])
    for e in events:
        if e["headline"]:
            key = (int(e["after_ref_pos"]), f"{e['event']}{e['after_ref_pos']}^{e['before_ref_pos']}")
            feats[key].add(e["accession"])
    idx = {a: i for i, a in enumerate(prot)}
    F = [(k, v) for k, v in feats.items() if GROUP_MIN <= len(v) <= N // 2]
    bits = [sum(1 << idx[a] for a in v) for _, v in F]
    sizes = [len(v) for _, v in F]
    parent = list(range(len(F)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(F)):
        for j in range(i + 1, len(F)):
            if min(sizes[i], sizes[j]) < GROUP_JACCARD * max(sizes[i], sizes[j]):
                continue
            inter = (bits[i] & bits[j]).bit_count()
            if inter / (sizes[i] + sizes[j] - inter) >= GROUP_JACCARD:
                parent[find(i)] = find(j)
    blocks = defaultdict(list)
    for i in range(len(F)):
        blocks[find(i)].append(i)
    grows, gmem = [], []
    for members in blocks.values():
        if len({F[i][0][0] for i in members}) < 2:
            continue
        cnt = Counter(a for i in members for a in F[i][1])
        core = sorted(a for a, n in cnt.items() if n >= CORE_FRAC * len(members))
        if len(core) < GROUP_MIN:
            continue
        groups = Counter(group_of(a) for a in core) if have_lineage else Counter()
        pls = [by_prot[a][F[i][0][0]]["pl"] for i in members for a in F[i][1]
               if a in core and F[i][0][0] in by_prot[a]]
        grows.append({"n_features": len(members), "n_positions": len({F[i][0][0] for i in members}),
                      "n_core_proteins": len(core),
                      "n_core_interpro_source": sum(source.get(a) == "interpro" for a in core),
                      "holds_in_interpro_only_set": sum(source.get(a) == "interpro" for a in core) >= GROUP_MIN,
                      "features": " ".join(F[i][0][1] for i in sorted(members, key=lambda i: F[i][0][0])),
                      "mean_plddt_core": mean(pls),
                      "n_lineage_groups": len(groups) if have_lineage else "",
                      "lineage_groups": "; ".join(f"{g}:{n}" for g, n in groups.most_common()) if have_lineage else "",
                      "spans_distinct_lineages": (len(groups) > 1) if have_lineage else "lineage needs Phase 2",
                      "organisms": "; ".join(sorted({org(a) for a in core}))[:600],
                      "_core": core})
    grows.sort(key=lambda r: (-r["n_positions"], -r["n_core_proteins"]))
    for i, r in enumerate(grows, 1):
        r["group_id"] = f"G{i:03d}"
        for a in r["_core"]:
            gmem.append({"group_id": r["group_id"], "accession": a, "organism": org(a), "group": group_of(a),
                         "tm_norm_ref": summ[a]["tm_norm_ref"]})
    write("p4_co_changing_groups.tsv", grows, ["group_id", "n_positions", "n_features", "n_core_proteins",
                                               "n_core_interpro_source", "holds_in_interpro_only_set",
                                               "spans_distinct_lineages", "n_lineage_groups", "lineage_groups",
                                               "mean_plddt_core", "features", "organisms"])
    write("p4_co_changing_members.tsv", gmem)

    print(f"{N} fold-passing proteins x {len(positions)} human positions")
    print(f"lineage from Phase 2: {'yes' if have_lineage else 'NO (run p2_fetch_annotations.py first for lineage columns)'}")
    print(f"positions ranked: {len(ranked_pos)}; regions: {len(chosen)}; indel hotspots: {len(hot)}; "
          f"unique changes: {len(uniq)}; co-changing groups: {len(grows)}")
    print("\nTop 10 most variable positions (high-confidence differences only):")
    for r in ranked_pos[:10]:
        print(f"  {r['human_aa']}{r['ref_pos']:<4} differ {r['n_differ_high_conf']:>3}  "
              f"(low-conf excluded {r['n_differ_low_conf_excluded']:>3})  identical {r['pct_identical_high_conf']}%  "
              f"pLDDT {r['mean_plddt']}  {r['label']}")
    print("\nTop 10 indel hotspots:")
    for h in hot[:10]:
        print(f"  {h['event']:<11} {h['human_region']:<9} {h['n_proteins_high_conf']:>3} proteins (hc), "
              f"len {h['length_min']}-{h['length_max']}, pLDDT {h['mean_plddt_high_conf']}"
              f"{'  [reference gap]' if h['at_reference_gap'] else ''}")
    plots(positions, pos_by, prot, by_prot, group_of, summ, ins_after, pocket)


def plots(positions, pos_by, prot, by_prot, group_of, summ, ins_after, pocket):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
    except ImportError:
        print("\nmatplotlib not installed: tables written, figures skipped. "
              "Install with: py -3.12 -m pip install matplotlib")
        return
    import numpy as np
    x = positions
    fig, ax = plt.subplots(4, 1, figsize=(16, 10), sharex=True)
    ax[0].plot(x, [pos_by[p]["pct_identical_high_conf"] or np.nan for p in x], lw=1, color="#1f5f8b")
    ax[0].set_ylabel("% identical\nto human (hc)")
    bottom = np.zeros(len(x))
    for col, c, lab in (("n_substituted_hc", "#3b7dd8", "substituted"), ("n_displaced_hc", "#e08a1e", "displaced >4 A"),
                        ("n_deleted_hc", "#222222", "deleted"), ("n_differ_low_conf_excluded", "#bbbbbb", "low pLDDT (excluded)")):
        v = np.array([pos_by[p][col] for p in x], float)
        ax[1].bar(x, v, bottom=bottom, width=1, color=c, label=lab)
        bottom += v
    ax[1].set_ylabel("proteins\ndiffering")
    ax[1].legend(fontsize=8, ncol=4, loc="upper right")
    ax[2].bar(x, [ins_after[p] for p in x], width=1, color="#c0392b")
    ax[2].set_ylabel("insertions\nafter (hc)")
    ax[3].plot(x, [pos_by[p]["mean_plddt"] or np.nan for p in x], lw=1, color="#555")
    ax[3].axhline(70, ls="--", lw=0.8, color="#c0392b")
    ax[3].set_ylabel("mean pLDDT")
    ax[3].set_xlabel("human AdoMetDC position (P17707 numbering; 1MSV frame)")
    for a in ax:
        for p in SITES:
            a.axvline(p, color="#8e44ad", lw=0.6, alpha=0.6)
        for p in pocket:
            a.axvline(p, color="#27ae60", lw=0.4, alpha=0.4)
        a.axvspan(21.5, 27.5, color="#eeeeee", zorder=0)
    ax[0].set_title("Where the 773 fold-passing proteins differ from human  "
                    "(purple = FOLDTRACE sites, green = within 6 A of putrescine, grey band = no 1MSV coordinates)",
                    fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "p4_change_map.png", dpi=150)
    plt.close(fig)

    code = {"identical": 0, "substituted": 1, "displaced": 2, "deleted": 3, "not_covered": 4}
    order = sorted(prot, key=lambda a: (group_of(a, 4), -float(summ[a]["tm_norm_ref"])))
    M = np.full((len(order), len(x)), 4)
    for r, a in enumerate(order):
        for c, p in enumerate(x):
            g = by_prot[a].get(p)
            if g is None:
                continue
            M[r, c] = 5 if (g["call"] in DIFF and not g["hc"]) else code[g["call"]]
    cmap = ListedColormap(["#f2f2f2", "#3b7dd8", "#e08a1e", "#222222", "#ffffff", "#c8c8c8"])
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.imshow(M, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=5,
              extent=(x[0] - 0.5, x[-1] + 0.5, len(order), 0))
    ax.set_xlabel("human position")
    ax.set_ylabel(f"{len(order)} proteins, sorted by taxonomic group (Phase 2 lineage if available), then TM-score")
    ax.set_title("Each row = one protein. light = same as human, blue = substituted, orange = displaced >4 A, "
                 "black = deleted, grey = difference at pLDDT < 70, white = not covered", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "p4_heatmap.png", dpi=150)
    plt.close(fig)
    print(f"\nfigures: {FIG / 'p4_change_map.png'}\n         {FIG / 'p4_heatmap.png'}")


if __name__ == "__main__":
    main()
