"""p21 C: rebuild the ODC-half tree inputs and test sets, run IQ-TREE if it is installed, and score support.

IQ-TREE is required for the tree; there is no NJ fallback. If no iqtree2 / iqtree / iqtree3 binary is on PATH
(or given with --iqtree PATH), the script writes the alignment, the test sets and the exact command, and stops.

Sequences
  - every row of results/p10_fusion_search.tsv (UniParc sequences from the p10 cache), ODC half = first PLP Pfam
    start to the C-terminus when AdoMetDC comes first, N-terminus to the last PLP Pfam end when the PLP domain comes
    first (as p20);
  - the p6/p10 reference decarboxylases (Q8IJ77 is already a fusion row), whole sequence;
  - every standalone PLP-decarboxylase of the four p10 species (data/p10 species caches), whole sequence.
Coverage filter (replaces p10's 250-residue Pfam-span filter): the segment must cover >= MIN_COV of human ODC
  residues 1-422, counted as aligned positions - the p10 structural count (odc_positions_aligned) where p10 has one,
  otherwise a local sequence alignment to P11926 (BLOSUM62, 11/1). Short Pfam spans no longer remove a protein.
Alignment: FAMSA (pyfamsa), columns with > 50 % gaps removed (as p10).
Test sets, from p21 B's core-site calls (results/p21b_rescore.tsv for the 23, p21b_rescore_p10.tsv for the rest):
  red-algal fusions; red-algal fusions + Percolomonas; degraded ODC halves (core rule); degraded ODC halves,
  eukaryotes only; green-lineage fusions; Alveolata fusions. Only sequences in the tree count.
Tree: iqtree2 -s <aln> -m MFP -bb 1000 -seed 1 -nt AUTO. Scoring (unrooted): a set is monophyletic if it is one
  side of a branch of the ML tree; its UFBoot support is then that branch's support. Otherwise the smallest side of
  any branch that contains the whole set is reported, with the extra members, and no support is given for the set.

Writes results/p21c_alignment.fasta, results/p21c_sets.tsv, results/p21c_report.md and, after a tree,
results/p21c_tree.treefile (IQ-TREE output prefix results/p21c_tree) and results/p21c_tree_support.tsv
Run:  py -3.12 p21c_tree.py                    (build + run IQ-TREE if found)
      py -3.12 p21c_tree.py --iqtree C:\\path\\iqtree2.exe
"""
import csv
import json
import shutil
import subprocess
import sys
from collections import defaultdict

from Bio import Phylo
from pyfamsa import Aligner, Sequence

import p6_decarboxylase_partner as p6
import p10_red_algal_pattern as p10
from p2_fetch_annotations import DATA, RESULTS
from p20_dnds_ci import regions_any

sys.stdout.reconfigure(encoding="utf8")
MIN_COV = 0.5            # share of human ODC 1-422
HUMAN_LEN = 422
PREFIX = RESULTS / "p21c_tree"
ALN = RESULTS / "p21c_alignment.fasta"
STANDALONE_TAXA = {83374: "Galdieria partita", 101924: "Rhodosorus marinus", 448386: "Gracilariopsis chorda",
                   63605: "Percolomonas cosmopolitus"}


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def seq_cov(seg, odc_seq):
    a = p10.LOCAL.align(seg, odc_seq[:HUMAN_LEN])[0]
    return sum(te - ts for ts, te in a.aligned[1])


def build():
    refs = p10.references()
    odc_seq = refs["P11926"]
    ents = json.loads(next((DATA / "p10").glob("uniparc__feature_id_PF01536*.json")).read_text(encoding="utf8"))
    seq_of = {e["uniParcId"]: e["sequence"]["value"] for e in ents}
    allf = {r["accession"].split(".")[0]: r for r in read("p10_all_fusions.tsv")}
    calls = {r["accession"]: r["call_core"] for r in read("p21b_rescore.tsv")}
    calls.update({r["accession"].split(".")[0]: r["call_core"] for r in read("p21b_rescore_p10.tsv")})
    rows, seqs = [], {}
    for r in read("p10_fusion_search.tsv"):
        acc = (r["in_our_23"] or (r["uniprot_accessions"].split(",")[0] if r["uniprot_accessions"] else r["upi"])).split(".")[0]
        s = seq_of[r["upi"]]
        try:
            lo, hi = regions_any(r["pfam_architecture"], len(s))["ODC"]
        except ValueError:
            rows.append({"id": acc, "kind": "fusion", "organism": r["organism"], "group": r["group"],
                         "in_tree": False, "why": "ODC half not definable from Pfam"})
            continue
        seg = s[lo - 1:hi]
        f = allf.get(acc, {})
        struct = f.get("odc_positions_aligned", "")
        cov = int(struct.split("/")[0]) if struct else seq_cov(seg, odc_seq)
        rows.append({"id": acc, "kind": "fusion", "organism": r["organism"], "group": r["group"],
                     "segment": f"{lo}-{hi}", "coverage_of_human_1_422": cov,
                     "coverage_method": "structure (p10)" if struct else "sequence",
                     "core_call": calls.get(acc, "no structural read"), "in_tree": cov >= MIN_COV * HUMAN_LEN,
                     "why": "" if cov >= MIN_COV * HUMAN_LEN else f"covers {cov} < {MIN_COV * HUMAN_LEN:.0f}"})
        if rows[-1]["in_tree"]:
            seqs[acc] = seg
    for acc, (label, _) in {**p6.REFS, **p10.EXTRA_REFS}.items():
        if acc in seqs or acc == "Q8IJ77":
            continue
        cov = seq_cov(refs[acc], odc_seq)
        rows.append({"id": acc, "kind": "reference", "organism": label, "group": "reference", "segment": "whole",
                     "coverage_of_human_1_422": cov, "coverage_method": "sequence",
                     "in_tree": cov >= MIN_COV * HUMAN_LEN, "why": "" if cov >= MIN_COV * HUMAN_LEN else "low coverage"})
        if rows[-1]["in_tree"]:
            seqs[acc] = refs[acc]
    for taxid, name in STANDALONE_TAXA.items():
        f = next((DATA / "p10").glob(f"uniparc_organism_id_{taxid}_AND_*.json"))
        for e in json.loads(f.read_text(encoding="utf8")):
            ids = {x[0] for x in p10.pfam_spans(e)}
            if not ids & p10.PLP or ids & p10.ADOMETDC:
                continue
            s = e["sequence"]["value"]
            cov = seq_cov(s, odc_seq)
            uid = e["uniParcId"]
            rows.append({"id": uid, "kind": "standalone", "organism": name, "group": "standalone PLP decarboxylase",
                         "segment": "whole", "coverage_of_human_1_422": cov, "coverage_method": "sequence",
                         "in_tree": cov >= MIN_COV * HUMAN_LEN,
                         "why": "" if cov >= MIN_COV * HUMAN_LEN else f"covers {cov} < {MIN_COV * HUMAN_LEN:.0f}"})
            if rows[-1]["in_tree"]:
                seqs[uid] = s
    ids = list(seqs)
    msa = Aligner().align([Sequence(i.encode(), seqs[i].encode()) for i in ids])
    aln = {s.id.decode(): s.sequence.decode() for s in msa}
    L = len(aln[ids[0]])
    keep = [k for k in range(L) if sum(aln[i][k] == "-" for i in ids) <= 0.5 * len(ids)]
    with open(ALN, "w", encoding="utf8") as fh:
        for i in ids:
            fh.write(f">{i}\n{''.join(aln[i][k] for k in keep)}\n")
    return rows, len(ids), len(keep)


def test_sets(rows):
    t = [r for r in rows if r["in_tree"] and r["kind"] == "fusion"]
    return {
        "red-algal fusions": [r["id"] for r in t if r["group"] == "Rhodophyta"],
        "red-algal fusions + Percolomonas": [r["id"] for r in t if r["group"] in ("Rhodophyta", "Heterolobosea")],
        "degraded ODC halves (core rule)": [r["id"] for r in t if r.get("core_call") == "degraded"],
        "degraded ODC halves, eukaryotes only": [r["id"] for r in t if r.get("core_call") == "degraded" and r["group"] != "Bacteria"],
        "green-lineage fusions": [r["id"] for r in t if r["group"] == "Viridiplantae"],
        "Alveolata fusions": [r["id"] for r in t if r["group"] == "Alveolata"],
    }


def find_iqtree():
    if "--iqtree" in sys.argv:
        return sys.argv[sys.argv.index("--iqtree") + 1]
    return next((shutil.which(b) for b in ("iqtree2", "iqtree3", "iqtree") if shutil.which(b)), None)


def best_model():
    for ln in open(str(PREFIX) + ".iqtree", encoding="utf8", errors="replace"):
        if ln.startswith("Best-fit model according to"):
            return ln.split(":", 1)[1].strip()
    return "?"


def score(sets, labels):
    """Unrooted: every internal branch splits the leaves in two; a 'clade' here is either side of a branch."""
    tree = Phylo.read(str(PREFIX) + ".treefile", "newick")
    leaves = {t.name for t in tree.get_terminals()}
    splits = []
    for cl in tree.get_nonterminals():
        side = frozenset(t.name for t in cl.get_terminals())
        if 1 < len(side) < len(leaves) - 1:
            sup = cl.confidence
            if sup is None and cl.name:
                try:
                    sup = float(cl.name.split("/")[-1])
                except ValueError:
                    sup = None
            splits.append((side, sup))
    sides = [(s, sup) for s, sup in splits] + [(frozenset(leaves - s), sup) for s, sup in splits]
    lab = lambda x: f"{x} ({labels.get(x, '?')})"
    out = []
    for name, mem in sets.items():
        m = set(mem) & leaves
        hit = next((sup for s, sup in sides if s == m), None) if len(m) > 1 else None
        if len(m) > 1 and any(s == m for s, _ in sides):
            out.append({"set": name, "n_in_tree": len(m), "monophyletic": True, "ufboot_of_set_node": hit,
                        "pieces": 1, "pure_subclades": "", "intruders_in_smallest_enclosing_clade": "",
                        "smallest_enclosing_clade_size": len(m), "smallest_enclosing_clade_ufboot": hit})
            continue
        # maximal pure pieces: largest existing clades made only of set members (singletons where none)
        pure = sorted({s for s, _ in sides if s <= m and len(s) > 1}, key=len, reverse=True)
        pieces, covered = [], set()
        for s in pure:
            if not (s & covered):
                pieces.append(s)
                covered |= s
        pieces += [frozenset([x]) for x in sorted(m - covered)]
        sup_of = {s: sup for s, sup in sides}
        enc = min((s for s, _ in sides if m <= s), key=len, default=None)
        out.append({"set": name, "n_in_tree": len(m), "monophyletic": False, "ufboot_of_set_node": "",
                    "pieces": len(pieces),
                    "pure_subclades": " | ".join(
                        (f"{len(p)} seqs, UFBoot {sup_of.get(p)}: " if len(p) > 1 else "singleton: ")
                        + ", ".join(lab(x) for x in sorted(p)) for p in pieces),
                    "intruders_in_smallest_enclosing_clade": "; ".join(lab(x) for x in sorted(enc - m)) if enc else "",
                    "smallest_enclosing_clade_size": len(enc) if enc else "",
                    "smallest_enclosing_clade_ufboot": sup_of.get(enc) if enc else ""})
    return out


def main():
    rows, n, L = build()
    sets = test_sets(rows)
    with open(RESULTS / "p21c_sets.tsv", "w", encoding="utf8", newline="") as f:
        cols = list(dict.fromkeys(k for r in rows for k in r))
        w = csv.DictWriter(f, cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    iq = find_iqtree()
    cmd = [iq or "iqtree2", "-s", str(ALN), "-m", "MFP", "-bb", "1000", "-seed", "1", "-nt", "AUTO",
           "-pre", str(PREFIX), "-redo"]
    support = None
    if iq or "--score-only" in sys.argv:
        if "--score-only" not in sys.argv:
            print("running:", " ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)
        labels = {r["id"]: f"{r['organism']}, {r['group']}" for r in rows}
        support = score(sets, labels)
        with open(RESULTS / "p21c_tree_support.tsv", "w", encoding="utf8", newline="") as f:
            cols = list(dict.fromkeys(k for r in support for k in r))
            w = csv.DictWriter(f, cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(support)
    report(rows, sets, n, L, iq, cmd, support)


def report(rows, sets, n, L, iq, cmd, support):
    fus = [r for r in rows if r["kind"] == "fusion"]
    out = [r for r in fus if not r["in_tree"]]
    g = {r["id"]: r for r in rows}
    R = ["# p21 C: ODC-half tree, rebuilt", "", "Generated by `scripts/p21c_tree.py` (method in its docstring).", ""]
    if not iq:
        R += ["**IQ-TREE is not installed on this machine, so no tree was built.** No NJ or other fallback was used. "
              "The alignment and test sets below are ready; once IQ-TREE 2 is installed, rerun the script (it finds "
              "`iqtree2` on PATH, or take `--iqtree <path>`), or run this command and then rerun the script to score it:", "",
              "```", " ".join(cmd), "```", "",
              "`p21c_tree_support.tsv` and the tree file are written only after that run.", ""]
    else:
        R += [f"IQ-TREE 2: `{' '.join(cmd)}`. ML tree with UFBoot support: `p21c_tree.treefile`; full log "
              "`p21c_tree.iqtree`.", "",
              f"**Model selected by ModelFinder (BIC): {best_model()}.**", "",
              "The tree is unrooted, so a \"clade\" here means one side of a branch. Monophyly and support are reported "
              "separately: a support value is given only for a node that exists in the ML tree. A set that is not a "
              "clade gets no support value of its own.", "",
              "## 1. Is each test set a clade?", "",
              "| set | sequences in tree | monophyletic | split into how many pure pieces | smallest clade holding the whole set (size) | sequences that break it up |",
              "|---|---|---|---|---|---|"]
        for s in support:
            R.append(f"| {s['set']} | {s['n_in_tree']} | **{'yes' if s['monophyletic'] else 'no'}** | {s['pieces']} | "
                     f"{s['smallest_enclosing_clade_size']} | "
                     f"{len(s['intruders_in_smallest_enclosing_clade'].split('; ')) if s['intruders_in_smallest_enclosing_clade'] else 0} intruders |")
        R += ["", "## 2. Support for nodes that exist in the ML tree", "",
              "| set | node | UFBoot |", "|---|---|---|"]
        for s in support:
            if s["monophyletic"]:
                R.append(f"| {s['set']} | the set itself | {s['ufboot_of_set_node']} |")
            else:
                for p in s["pure_subclades"].split(" | "):
                    if p.startswith("singleton"):
                        continue
                    head, mem = p.split(": ", 1)
                    R.append(f"| {s['set']} | pure sub-clade of {head.split(',')[0]} | {head.split('UFBoot ')[1]} |")
                R.append(f"| {s['set']} | smallest clade holding the whole set ({s['smallest_enclosing_clade_size']} "
                         f"sequences, includes intruders) | {s['smallest_enclosing_clade_ufboot']} |")
        R += ["", "UFBoot >= 95 is the usual threshold for a well-supported node in IQ-TREE's ultrafast bootstrap.", "",
              "## 3. Sets that are not clades: how they break up", ""]
        for s in support:
            if s["monophyletic"]:
                continue
            R += [f"### {s['set']}", "", "Pure pieces (largest clades made only of set members; singletons sit apart):", ""]
            R += [f"- {p}" for p in s["pure_subclades"].split(" | ")]
            intr = s["intruders_in_smallest_enclosing_clade"]
            R += ["", f"Smallest clade holding all of them: {s['smallest_enclosing_clade_size']} sequences "
                      f"(UFBoot {s['smallest_enclosing_clade_ufboot']}). Sequences inside it that are not in the set "
                      f"({len(intr.split('; ')) if intr else 0}):", ""]
            ilist = intr.split("; ") if intr else []
            R += [f"- {x}" for x in ilist[:60]]
            if len(ilist) > 60:
                R.append(f"- ... and {len(ilist) - 60} more (full list in p21c_tree_support.tsv)")
            R.append("")
    R += [f"Alignment: {n} sequences, {L} columns after removing columns with > 50 % gaps (`p21c_alignment.fasta`).", "",
          f"Coverage filter: segment must cover >= {MIN_COV:.0%} of human ODC 1-422 ({MIN_COV * HUMAN_LEN:.0f} positions). "
          f"Fusions removed: {len(out)} of {len(fus)}"
          + (": " + ", ".join(f"{r['id']} ({r['organism']}, {r['why']})" for r in out) if out else "") + ".", "",
          "Proteins the old 250-residue filter dropped, now: "
          + "; ".join(f"{a} {g[a]['organism']}: coverage {g[a].get('coverage_of_human_1_422')} "
                      f"({g[a].get('coverage_method')}), {'in the tree' if g[a]['in_tree'] else 'still out'}"
                      for a in ("A0A9C7PQ51", "A0A9C7PSW2", "M2YAF1", "A0A1X6PGY4") if a in g) + ".", "",
          "## Test sets (from p21 B's core-site calls)", ""]
    for k, v in sets.items():
        R.append(f"- **{k}** ({len(v)}): " + ", ".join(f"{a} ({g[a]['organism']})" for a in v))
    u = g.get("UPI0017429AE9")
    R += ["", "## UPI0017429AE9", "",
          "Not added to any fusion set. It is a *Rhodosorus marinus* standalone PLP decarboxylase (355 aa, "
          "PF02784 5-77 > PF00278 221-316, no AdoMetDC domain; p10: nearest AZIN1, 28 % identity to human ODC, 3/13 "
          "sites), so it cannot be a member of a set of fusions. It "
          + ("is in the tree as a standalone" if u and u["in_tree"] else
             f"is not in the tree (coverage {u.get('coverage_of_human_1_422') if u else '?'} of human ODC 1-422)")
          + "; in p10 it was the one extra member of the smallest clade holding the red-algal fusions.", "",
          "## Disagreements with the task's premises", "",
          "- Q8IJ77 leaves the degraded set, as the task expected: its only substitution is Y331F, which p21 B drops "
          "from the core.",
          "- A0A2E3ETT7 (*Tistrella*, a bacterium) stays degraded under the core list: it has lost H197 (H197N) as well "
          "as Y331. Keeping K69 and C360 is not enough by the rule. Because it is a bacterium with a different "
          "architecture (PF02675 before the PLP domain), a eukaryote-only degraded set is also scored.", ""]
    (RESULTS / "p21c_report.md").write_text("\n".join(R) + "\n", encoding="utf8")


if __name__ == "__main__":
    main()
