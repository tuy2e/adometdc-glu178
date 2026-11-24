"""p23 part 2: re-check AlphaFold DB for the p10 fusions that had no model, download any that now exist, and list
what is still missing for Colab, red algae and other non-green lineages first.

Rows: every row of results/p10_all_fusions.tsv with no AF-<acc>.pdb in data/p10_models/ and not one of our 23
(those have Marc's models). For a UniParc-only row (UPI...), UniParc is re-queried for UniProtKB accessions assigned
since p10 (AFDB takes only UniProt accessions; a UPI is rejected with HTTP 400). Each accession is tried at
https://alphafold.ebi.ac.uk/api/prediction/<acc>; a model is kept only if AFDB's sequence equals our UniParc
sequence. New models go to data/p23_models/ (data/p10_models is left as it was). API answers cached in data/p23/afdb/.

Writes results/p23_unmodelled.tsv, results/p23_unmodelled.md
Run:  py -3.12 p23_afdb_recheck.py
"""
import csv
import json
import sys
import time
import urllib.error
import urllib.request

from p2_fetch_annotations import DATA, RESULTS

sys.stdout.reconfigure(encoding="utf8")
UA = {"User-Agent": "adometdc-dataset-survey/1.0"}
CACHE = DATA / "p23" / "afdb"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = DATA / "p23_models"
OUT.mkdir(parents=True, exist_ok=True)
P10 = DATA / "p10_models"
PRIORITY = ["Rhodophyta", "Heterolobosea", "Glaucocystophyceae", "Stramenopiles", "other", "Alveolata", "Bacteria",
            "Metazoa", "Viridiplantae"]


def fetch(url, cache_name):
    p = CACHE / cache_name
    if p.exists():
        code, body = p.read_text(encoding="utf8").split("\n", 1)
        return int(code), body
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
            code, body = r.status, r.read().decode("utf8")
    except urllib.error.HTTPError as e:
        code, body = e.code, e.read().decode("utf8", "replace")
    p.write_text(f"{code}\n{body}", encoding="utf8")
    time.sleep(0.3)
    return code, body


def read(name):
    with open(RESULTS / name, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    allf = read("p10_all_fusions.tsv")
    search = {}
    for r in read("p10_fusion_search.tsv"):
        for k in [r["upi"]] + [a.split(".")[0] for a in r["uniprot_accessions"].split(",") if a]:
            search[k] = r
    ents = json.loads(next((DATA / "p10").glob("uniparc__feature_id_PF01536*.json")).read_text(encoding="utf8"))
    seq_of = {e["uniParcId"]: e["sequence"]["value"] for e in ents}
    rows = []
    for f in allf:
        if f["set"] == "our 23":
            continue
        accs = [a.split(".")[0] for a in f["accession"].split(",")]
        if any((P10 / f"AF-{a}.pdb").exists() for a in accs):
            continue
        s = search[accs[0]]
        upi, seq = s["upi"], seq_of[s["upi"]]
        tried = [a for a in accs if not a.startswith("UPI")]
        if not tried:   # UniParc-only in p10: any UniProt accession since?
            code, body = fetch(f"https://rest.uniprot.org/uniparc/{upi}?format=json", f"{upi}.uniparc.json")
            if code == 200:
                tried = [a.split(".")[0] for a in json.loads(body).get("uniProtKBAccessions", []) or []]
        status, model_acc, note = "no model", "", ""
        for a in tried:
            code, body = fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{a}", f"{a}.afdb.json")
            if code != 200:
                note += f"{a}: HTTP {code}; "
                continue
            ent = json.loads(body)[0]
            mseq = ent.get("uniprotSequence") or ent.get("sequence")
            if mseq and mseq != seq:
                note += f"{a}: AFDB model sequence differs from this UniParc sequence; "
                continue
            pdb = OUT / f"AF-{a}.pdb"
            if not pdb.exists():
                with urllib.request.urlopen(urllib.request.Request(ent["pdbUrl"], headers=UA), timeout=120) as r:
                    pdb.write_bytes(r.read())
            status, model_acc = "NEW model downloaded", a
            break
        if not tried:
            note = "UniParc-only (no UniProt accession), AFDB cannot be queried"
        rows.append({"accession": f["accession"], "upi": upi, "organism": f["organism"], "group": f["group"],
                     "length": len(seq), "adometdc_then_plp": s["adometdc_then_plp"], "accessions_tried": " ".join(tried),
                     "afdb_status": status, "model": f"data/p23_models/AF-{model_acc}.pdb" if model_acc else "",
                     "note": note.strip("; "), "colab_dimer_residues": 2 * len(seq)})
        print(f"{f['accession'][:20]:20s} {f['organism'][:30]:30s} {status}", flush=True)
    rank = {g: k for k, g in enumerate(PRIORITY)}
    rows.sort(key=lambda r: (r["afdb_status"] != "no model", rank.get(r["group"], len(PRIORITY)), r["length"]))
    with open(RESULTS / "p23_unmodelled.tsv", "w", encoding="utf8", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    new = [r for r in rows if r["afdb_status"] != "no model"]
    miss = [r for r in rows if r["afdb_status"] == "no model"]
    from collections import Counter
    L = ["# p23: fusions still without a structure model", "",
         "Generated by `scripts/p23_afdb_recheck.py` (method in its docstring). Our 23 have Marc's models and are not "
         "listed; the 56 p10 fusions already in data/p10_models are not listed either.", "",
         f"Checked: {len(rows)} fusions without a model. **New in AlphaFold DB since p10: {len(new)}** "
         f"(downloaded to `data/p23_models/`). **Still missing: {len(miss)}.**", ""]
    if new:
        L += ["## Newly available", "", "| accession | organism | group | length | model |", "|---|---|---|---|---|"]
        L += [f"| {r['accession']} | *{r['organism']}* | {r['group']} | {r['length']} | `{r['model']}` |" for r in new]
        L.append("")
    reasons = Counter("UniParc-only, no UniProt accession" if r["note"].startswith("UniParc-only") else
                      "sequence differs" if "differs" in r["note"] else "AFDB 404 (not in AFDB)" for r in miss)
    L += ["## Still missing, by priority", "",
          "Why: " + "; ".join(f"{v} {k}" for k, v in reasons.most_common()) + ".", "",
          "Priority order: " + " > ".join(PRIORITY) + "; within a group, shortest first. Monomer length is given; a "
          "homodimer job is twice that. Colab's free GPU limits are not fixed; very long jobs (several thousand "
          "residues as a dimer) may need a paid tier or a monomer-only run.", "",
          "| # | accession | organism | group | length | domain order AdoMetDC first | why no model |",
          "|---|---|---|---|---|---|---|"]
    for k, r in enumerate(miss, 1):
        L.append(f"| {k} | {r['accession']} | *{r['organism']}* | {r['group']} | {r['length']} | {r['adometdc_then_plp']} | "
                 f"{r['note'] or '-'} |")
    L += ["", "By group: " + ", ".join(f"{g} {n}" for g, n in Counter(r["group"] for r in miss).most_common()) + ".",
          "", "The Metazoa row (*Piliocolobus*) is parasite contamination (p21d) and does not need folding as a "
          "metazoan protein.", ""]
    (RESULTS / "p23_unmodelled.md").write_text("\n".join(L), encoding="utf8")


if __name__ == "__main__":
    main()
