"""Phase 2b: explain and, where possible, recover the UniProt accessions that are no longer active.

For every accession that p2_fetch_annotations.py cached as inactive:
  1. UniProt's own reason: inactiveReasonType (DELETED / MERGED / DEMERGED), deletedReason, and
     mergeDemergeTo (the replacement accession(s) UniProt names for merged/demerged entries).
  2. UniParc record of the same sequence (UPI from the inactive record), fetched in batches from
     /uniparc/upis. It gives:
       - every UniProtKB accession that carries the identical sequence (a deleted entry often has an
         active twin, e.g. after a proteome change); these are "same-sequence" successors, which is
         NOT the same as a UniProt merge and is labelled separately;
       - the sequence itself (checked against the AlphaFold model);
       - sequence features (Pfam, Gene3D, SUPERFAMILY, PANTHER, NCBIfam, PROSITE, ...) with positions.
  3. Replacement records: UniProt-named replacements (merge/demerge) and same-sequence active
     accessions are fetched from UniProt in batches, and their InterPro matches per accession.

Nothing already downloaded is changed or dropped. RESUMABLE: everything is cached under data/
(uniparc/, uniprot_successors/, interpro_successors/); re-runs skip cached items. Retries with
backoff come from p2_fetch_annotations.get().

Writes
    results/inactive_resolution.tsv   one row per inactive accession
    results/domain_hits_uniparc.tsv   UniParc sequence features for the inactive accessions
    results/p2b_manifest.json         counts, releases, dates
Run:
    py -3.12 p2b_resolve_inactive.py              (full; ~5 min, mostly InterPro for successors)
    py -3.12 p2b_resolve_inactive.py --no-interpro (UniProt/UniParc only, ~1 min)
"""
import argparse
import csv
import json
import time
from collections import Counter
from pathlib import Path

from p2_fetch_annotations import (DATA, IP_DIR, INTERPRO_ONE, POLITE_DELAY, RESULTS, UNIPROT_BATCH, UP_DIR,
                                  get, now, read_candidates, save_json)

UPARC_DIR, SUCC_DIR, SUCC_IP_DIR = DATA / "uniparc", DATA / "uniprot_successors", DATA / "interpro_successors"
UNIPARC_BATCH = "https://rest.uniprot.org/uniparc/upis?upis={}&format=json"
BATCH = 100


def versions(h):
    return {"release": h.get("X-UniProt-Release"), "release_date": h.get("X-UniProt-Release-Date")}


def fetch_uniparc(upis):
    todo = sorted(u for u in upis if not (UPARC_DIR / f"{u}.json").exists())
    print(f"UniParc: {len(upis) - len(todo)} cached, {len(todo)} to fetch")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        st, h, body = get(UNIPARC_BATCH.format(",".join(chunk)))
        data = json.loads(body) if st == 200 else []
        if isinstance(data, dict):          # the API wraps the list when asked for application/json
            data = data.get("results", [])
        got = {r["uniParcId"]: r for r in data}
        for u in chunk:
            save_json(UPARC_DIR / f"{u}.json", {"uniparc_id": u, "downloaded_utc": now(), "versions": versions(h),
                                                "found": u in got, "entry": got.get(u)})
        print(f"  batch {i // BATCH + 1}: {len(got)}/{len(chunk)} found", flush=True)


def fetch_successors(accs):
    todo = sorted(a for a in accs if not (SUCC_DIR / f"{a}.json").exists())
    print(f"Successor UniProt records: {len(accs) - len(todo)} cached, {len(todo)} to fetch")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        st, h, body = get(UNIPROT_BATCH.format(",".join(chunk)))
        got = {e["primaryAccession"]: e for e in (json.loads(body).get("results", []) if st == 200 else [])}
        for a in chunk:
            save_json(SUCC_DIR / f"{a}.json", {"accession": a, "downloaded_utc": now(), "versions": versions(h),
                                               "status": "active" if a in got else "not_returned", "entry": got.get(a)})


def fetch_successor_interpro(accs):
    todo = sorted(a for a in accs if not (SUCC_IP_DIR / f"{a}.json").exists())
    print(f"Successor InterPro: {len(accs) - len(todo)} cached, {len(todo)} to fetch")
    for n, a in enumerate(todo, 1):
        url, results, ver, status = INTERPRO_ONE.format(a), [], None, "ok"
        try:
            while url:
                st, h, body = get(url)
                ver = ver or {"interpro_version": h.get("InterPro-Version")}
                if st == 204 or not body.strip():
                    status = "no_matches"
                    break
                if st in (404, 410):
                    status = f"not_found_http_{st}"
                    break
                page = json.loads(body)
                results.extend(page.get("results", []))
                url = page.get("next")
                time.sleep(POLITE_DELAY)
        except RuntimeError as e:
            print(f"  {a}: {e}")
            continue
        save_json(SUCC_IP_DIR / f"{a}.json", {"accession": a, "downloaded_utc": now(), "versions": ver,
                                              "status": status, "results": results})
        if n % 25 == 0 or n == len(todo):
            print(f"  {n}/{len(todo)}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-interpro", action="store_true")
    args = ap.parse_args()
    for d in (UPARC_DIR, SUCC_DIR, SUCC_IP_DIR):
        d.mkdir(parents=True, exist_ok=True)
    seqs = read_candidates()
    fold_ok = {}
    with open(RESULTS / "master_candidates.tsv", encoding="utf8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fold_ok[r["accession"]] = r["fold_ok"]

    inactive = {}
    for a in seqs:
        c = json.loads((UP_DIR / f"{a}.json").read_text(encoding="utf8"))
        if c["status"] != "active":
            inactive[a] = c
    print(f"{len(inactive)} inactive accessions")

    upi = {a: ((c.get("entry") or {}).get("extraAttributes") or {}).get("uniParcId") for a, c in inactive.items()}
    fetch_uniparc({u for u in upi.values() if u})

    rows, feat_rows, succ_needed = [], [], {}
    for a, c in inactive.items():
        e = c.get("entry") or {}
        r = e.get("inactiveReason") or {}
        named = r.get("mergeDemergeTo") or []
        up = json.loads((UPARC_DIR / f"{upi[a]}.json").read_text(encoding="utf8")) if upi[a] else None
        ue = (up or {}).get("entry") or {}
        same_seq = [x for x in ue.get("uniProtKBAccessions", []) if "." not in x and x != a]
        succ_needed[a] = named or same_seq
        rows.append({"accession": a, "fold_ok": fold_ok[a], "uniprot_status": c["status"],
                     "reason_type": r.get("inactiveReasonType", ""), "deleted_reason": r.get("deletedReason", ""),
                     "uniprot_named_replacement": ";".join(named), "uniparc_id": upi[a] or "",
                     "uniparc_found": bool(ue),
                     "uniparc_seq_equals_model": (ue.get("sequence", {}).get("value") == seqs[a]) if ue else "",
                     "same_sequence_uniprot_accessions": ";".join(same_seq),
                     "common_taxon": ";".join(t.get("commonTaxon", "") for t in ue.get("commonTaxons", [])),
                     "n_uniparc_features": len(ue.get("sequenceFeatures", [])),
                     "uniprot_release": (c.get("versions") or {}).get("release", ""),
                     "uniparc_release": ((up or {}).get("versions") or {}).get("release", ""),
                     "downloaded_utc": (up or {}).get("downloaded_utc", "")})
        for f in ue.get("sequenceFeatures", []):
            for loc in f.get("locations", []):
                feat_rows.append({"accession": a, "uniparc_id": upi[a], "source_database": f.get("database", ""),
                                  "entry_accession": f.get("databaseId", ""),
                                  "interpro_group": (f.get("interproGroup") or {}).get("id", ""),
                                  "interpro_name": (f.get("interproGroup") or {}).get("name", ""),
                                  "start": loc.get("start"), "end": loc.get("end"),
                                  "uniparc_release": ((up or {}).get("versions") or {}).get("release", "")})

    all_succ = sorted({s for v in succ_needed.values() for s in v})
    fetch_successors(all_succ)
    succ = {s: json.loads((SUCC_DIR / f"{s}.json").read_text(encoding="utf8")) for s in all_succ}
    if not args.no_interpro:
        fetch_successor_interpro([s for s in all_succ if succ[s]["status"] == "active"])

    for row in rows:
        cands = [s for s in succ_needed[row["accession"]] if succ[s]["status"] == "active"]
        chosen = cands[0] if cands else ""
        se = succ[chosen]["entry"] if chosen else {}
        desc = (se or {}).get("proteinDescription", {})
        rec = desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
        ipf = SUCC_IP_DIR / f"{chosen}.json"
        row["replacement_used"] = chosen
        row["replacement_kind"] = ("uniprot_merge_or_demerge" if chosen and row["uniprot_named_replacement"]
                                   else "same_sequence_active_entry" if chosen else "")
        row["replacement_seq_equals_model"] = ((se or {}).get("sequence", {}).get("value") ==
                                               seqs[row["accession"]]) if chosen else ""
        row["replacement_protein_name"] = (rec.get("fullName") or {}).get("value", "") if chosen else ""
        row["replacement_organism"] = (se or {}).get("organism", {}).get("scientificName", "") if chosen else ""
        row["replacement_interpro"] = (json.loads(ipf.read_text(encoding="utf8"))["status"]
                                       if chosen and ipf.exists() else "")
        row["annotation_recovered"] = ("full (replacement UniProt + InterPro)" if row["replacement_interpro"] == "ok"
                                       else "replacement UniProt record" if chosen
                                       else "domains only (UniParc features)" if row["n_uniparc_features"]
                                       else "none")

    fields = ["accession", "fold_ok", "uniprot_status", "reason_type", "deleted_reason", "uniprot_named_replacement",
              "same_sequence_uniprot_accessions", "replacement_used", "replacement_kind",
              "replacement_seq_equals_model", "replacement_protein_name", "replacement_organism",
              "replacement_interpro", "uniparc_id", "uniparc_found", "uniparc_seq_equals_model", "common_taxon",
              "n_uniparc_features", "annotation_recovered", "uniprot_release", "uniparc_release", "downloaded_utc"]
    with open(RESULTS / "inactive_resolution.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if feat_rows:
        with open(RESULTS / "domain_hits_uniparc.tsv", "w", newline="", encoding="utf8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(feat_rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(feat_rows)

    def tab(key, subset=rows):
        return dict(Counter(r[key] for r in subset))
    fo = [r for r in rows if r["fold_ok"] == "true"]
    man = {"run_finished_utc": now(), "n_inactive": len(rows), "n_inactive_fold_ok": len(fo),
           "reason_type": tab("reason_type"), "deleted_reason": tab("deleted_reason"),
           "deleted_reason_fold_ok": tab("deleted_reason", fo),
           "uniparc_seq_equals_model": tab("uniparc_seq_equals_model"),
           "replacement_kind": tab("replacement_kind"), "replacement_kind_fold_ok": tab("replacement_kind", fo),
           "annotation_recovered": tab("annotation_recovered"),
           "annotation_recovered_fold_ok": tab("annotation_recovered", fo),
           "uniprot_releases": tab("uniparc_release"), "interpro_skipped": args.no_interpro}
    save_json(RESULTS / "p2b_manifest.json", man)
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
