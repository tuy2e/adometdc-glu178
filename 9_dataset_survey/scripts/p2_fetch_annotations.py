"""Phase 2, step 1: download UniProt and InterPro annotation for all 1,000 candidates.

RESUMABLE. Every response is cached under 9_dataset_survey/data/ as soon as it arrives
(written to a temporary file, then renamed, so an interrupted run never leaves a half-written
cache file). On a re-run, accessions that are already cached are skipped. Failed requests are
retried with exponential backoff, honouring the server's Retry-After header. Anything that still
fails is listed at the end; just run the script again to pick it up.

Sources
    UniProtKB REST  batch endpoint /uniprotkb/accessions (100 accessions per request).
                    Accessions missing from a batch reply (deleted, merged, demerged) are then
                    fetched one by one, which returns the inactive-entry record and its reason.
    InterPro API    /api/entry/all/protein/uniprot/<acc>/ : InterPro entries AND every member
                    database (Pfam, CDD, SMART, CATH-Gene3D, SUPERFAMILY, PANTHER, ...) with
                    residue positions. InterPro has no multi-accession endpoint, so this is one
                    request per accession (paginated if needed).

Versions and dates
    Every cache file stores the UTC download time and the release reported by the server
    (UniProt: X-UniProt-Release / X-UniProt-Release-Date headers; InterPro: InterPro-Version
    header). The output tables carry these per row, and results/p2_download_manifest.json
    summarises them, flagging if a resumed run mixed releases.

Reads (never writes) the FOLDTRACE run folder, for the list of accessions and the model sequences.

Writes
    data/uniprot/<acc>.json        cached UniProt entry (or inactive record / not-found marker)
    data/interpro/<acc>.json       cached InterPro matches (all pages merged)
    results/annotation_uniprot.tsv one row per candidate
    results/domain_hits.tsv        one row per domain/family/site match fragment
    results/p2_download_manifest.json

Run (expect roughly 15-30 minutes, mostly the 1,000 InterPro requests):
    py -3.12 p2_fetch_annotations.py
Options:
    --only uniprot | interpro      fetch one source only
    --limit N                      only the first N accessions (for a quick test)
"""
import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

RUN = Path(r"C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc")
BASE = Path(__file__).resolve().parent.parent
DATA, RESULTS = BASE / "data", BASE / "results"
UP_DIR, IP_DIR = DATA / "uniprot", DATA / "interpro"

UNIPROT_BATCH = "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&format=json&size=500"
UNIPROT_ONE = "https://rest.uniprot.org/uniprotkb/{}?format=json"
INTERPRO_ONE = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{}/?page_size=200"
BATCH = 100
MAX_TRIES = 6
POLITE_DELAY = 0.2       # seconds between InterPro requests
USER_AGENT = "adometdc-dataset-survey/1.0 (Python urllib; academic research)"


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def save_json(path: Path, obj) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj), encoding="utf8")
    os.replace(tmp, path)


def get(url: str):
    """GET with retries. Returns (status, headers, body_text). Status 404/204 are returned, not retried."""
    for attempt in range(1, MAX_TRIES + 1):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.status, dict(r.headers), r.read().decode("utf8")
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                return e.code, dict(e.headers), e.read().decode("utf8", "replace")
            wait = 2 ** attempt
            ra = e.headers.get("Retry-After")
            if ra and ra.isdigit():
                wait = max(wait, int(ra))
            reason = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            wait, reason = 2 ** attempt, type(e).__name__
        if attempt < MAX_TRIES:
            print(f"    {reason}; retry {attempt}/{MAX_TRIES - 1} in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"gave up after {MAX_TRIES} tries: {url}")


def read_candidates() -> dict:
    """accession -> model sequence, in Foldseek order."""
    seqs = {}
    with open(RUN / "foldseek_hits.m8", encoding="utf8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) == 21:
                seqs[re.search(r"AF-([A-Z0-9]+)-F\d+", row[1]).group(1)] = row[18]
    return seqs


# ------------------------------------------------------------------ UniProt
def fetch_uniprot(accs: list, failed: list) -> None:
    todo = [a for a in accs if not (UP_DIR / f"{a}.json").exists()]
    print(f"UniProt: {len(accs) - len(todo)} cached, {len(todo)} to fetch")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        print(f"  batch {i // BATCH + 1}/{(len(todo) + BATCH - 1) // BATCH} ({len(chunk)} accessions)", flush=True)
        try:
            status, hdr, body = get(UNIPROT_BATCH.format(",".join(chunk)))
        except RuntimeError as e:
            print(f"  {e}")
            failed.extend(("uniprot", a) for a in chunk)
            continue
        ver = {"release": hdr.get("X-UniProt-Release"), "release_date": hdr.get("X-UniProt-Release-Date"),
               "api_deployment_date": hdr.get("X-API-Deployment-Date")}
        got = {}
        if status == 200:
            for e in json.loads(body).get("results", []):
                got[e["primaryAccession"]] = e
        stamp = now()
        for a in chunk:
            if a in got:
                save_json(UP_DIR / f"{a}.json", {"accession": a, "downloaded_utc": stamp, "versions": ver,
                                                 "status": "active", "url": "batch", "entry": got[a]})
        # Missing from the batch: fetch singly to get the inactive record (or a 404)
        for a in (x for x in chunk if x not in got):
            try:
                st, h1, b1 = get(UNIPROT_ONE.format(a))
            except RuntimeError as e:
                print(f"    {e}")
                failed.append(("uniprot", a))
                continue
            v1 = {"release": h1.get("X-UniProt-Release"), "release_date": h1.get("X-UniProt-Release-Date"),
                  "api_deployment_date": h1.get("X-API-Deployment-Date")}
            entry = json.loads(b1) if st == 200 and b1.strip() else None
            if entry and entry.get("primaryAccession") not in (None, a):
                status_s = f"redirected_to_{entry['primaryAccession']}"
            elif entry and entry.get("entryType", "").lower().startswith("inactive"):
                status_s = "inactive"
            elif entry:
                status_s = "active"
            else:
                status_s = f"not_found_http_{st}"
            save_json(UP_DIR / f"{a}.json", {"accession": a, "downloaded_utc": now(), "versions": v1,
                                             "status": status_s, "url": "single", "entry": entry})
            print(f"    {a}: {status_s}")


# ------------------------------------------------------------------ InterPro
def fetch_interpro(accs: list, failed: list) -> None:
    todo = [a for a in accs if not (IP_DIR / f"{a}.json").exists()]
    print(f"InterPro: {len(accs) - len(todo)} cached, {len(todo)} to fetch")
    t0 = time.time()
    for n, a in enumerate(todo, 1):
        url, results, version, status_s = INTERPRO_ONE.format(a), [], None, "ok"
        try:
            while url:
                st, hdr, body = get(url)
                version = version or {"interpro_version": hdr.get("InterPro-Version"),
                                      "interpro_version_minor": hdr.get("InterPro-Version-Minor")}
                if st == 204 or not body.strip():
                    status_s = "no_matches"
                    break
                if st in (404, 410):
                    status_s = f"not_found_http_{st}"
                    break
                page = json.loads(body)
                results.extend(page.get("results", []))
                url = page.get("next")
                time.sleep(POLITE_DELAY)
        except RuntimeError as e:
            print(f"  {a}: {e}")
            failed.append(("interpro", a))
            continue
        save_json(IP_DIR / f"{a}.json", {"accession": a, "downloaded_utc": now(), "versions": version,
                                         "status": status_s, "results": results})
        if n % 25 == 0 or n == len(todo):
            rate = (time.time() - t0) / n
            print(f"  {n}/{len(todo)} done, about {rate * (len(todo) - n) / 60:.0f} min left", flush=True)


# ------------------------------------------------------------------ tables
def build_tables(accs: list, seqs: dict) -> dict:
    up_rows, dom_rows = [], []
    up_rel, ip_rel, statuses = Counter(), Counter(), Counter()
    for a in accs:
        p = UP_DIR / f"{a}.json"
        c = json.loads(p.read_text(encoding="utf8")) if p.exists() else None
        e = (c or {}).get("entry") or {}
        active = c is not None and c["status"] == "active"
        statuses[c["status"] if c else "not_downloaded"] += 1
        if c:
            up_rel[(c["versions"] or {}).get("release")] += 1
        desc = e.get("proteinDescription", {})
        rec = desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
        seq = (e.get("sequence") or {}).get("value", "")
        org = e.get("organism", {})
        xref = e.get("uniProtKBCrossReferences", [])
        up_rows.append({
            "accession": a,
            "uniprot_status": c["status"] if c else "not_downloaded",
            "inactive_reason": json.dumps(e.get("inactiveReason")) if e.get("inactiveReason") else "",
            "entry_type": e.get("entryType", ""),
            "entry_name": e.get("uniProtkbId", ""),
            "protein_name": (rec.get("fullName") or {}).get("value", ""),
            "fragment_flag": (desc.get("flag") or ""),
            "gene": ";".join(g.get("geneName", {}).get("value", "") for g in e.get("genes", []) if g.get("geneName")),
            "organism": org.get("scientificName", ""),
            "taxid": org.get("taxonId", ""),
            "lineage": ";".join(org.get("lineage", [])),
            "protein_existence": e.get("proteinExistence", ""),
            "annotation_score": e.get("annotationScore", ""),
            "uniprot_length": len(seq) if seq else "",
            "model_length": len(seqs[a]),
            "sequence_version": (e.get("entryAudit") or {}).get("sequenceVersion", ""),
            "last_sequence_update": (e.get("entryAudit") or {}).get("lastSequenceUpdateDate", ""),
            "uniprot_seq_equals_model": (seq == seqs[a]) if seq else "",
            "ec_numbers": ";".join(x.get("value", "") for x in (rec.get("ecNumbers") or [])),
            "pfam_ids": ";".join(x["id"] for x in xref if x["database"] == "Pfam"),
            "interpro_ids": ";".join(x["id"] for x in xref if x["database"] == "InterPro"),
            "keywords": ";".join(k.get("name", "") for k in e.get("keywords", [])),
            "uniprot_release": (c["versions"] or {}).get("release", "") if c else "",
            "uniprot_downloaded_utc": c["downloaded_utc"] if c else "",
        })
        q = IP_DIR / f"{a}.json"
        ip = json.loads(q.read_text(encoding="utf8")) if q.exists() else None
        up_rows[-1]["interpro_status"] = ip["status"] if ip else "not_downloaded"
        if not ip:
            continue
        ipv = (ip["versions"] or {}).get("interpro_version", "")
        ip_rel[ipv] += 1
        for r in ip["results"]:
            md = r["metadata"]
            integrated = md.get("integrated") or ""
            for prot in r.get("proteins", []):
                for loc in prot.get("entry_protein_locations") or []:
                    for k, fr in enumerate(loc.get("fragments", [])):
                        dom_rows.append({
                            "accession": a,
                            "source_database": md["source_database"],
                            "entry_accession": md["accession"],
                            "entry_name": md.get("name", ""),
                            "entry_type": md.get("type", ""),
                            "integrated_into": integrated if isinstance(integrated, str) else json.dumps(integrated),
                            "start": fr.get("start"),
                            "end": fr.get("end"),
                            "fragment_index": k,
                            "n_fragments": len(loc.get("fragments", [])),
                            "dc_status": fr.get("dc-status", ""),
                            "score": loc.get("score", ""),
                            "model": loc.get("model", ""),
                            "protein_length_interpro": prot.get("protein_length", ""),
                            "interpro_version": ipv,
                            "interpro_downloaded_utc": ip["downloaded_utc"],
                        })
    RESULTS.mkdir(parents=True, exist_ok=True)
    for name, rows in (("annotation_uniprot.tsv", up_rows), ("domain_hits.tsv", dom_rows)):
        if rows:
            with open(RESULTS / name, "w", newline="", encoding="utf8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
    return {"uniprot_status": dict(statuses), "uniprot_releases": dict(up_rel),
            "interpro_versions": dict(ip_rel),
            "interpro_status": dict(Counter(r["interpro_status"] for r in up_rows)),
            "uniprot_seq_equals_model": dict(Counter(str(r["uniprot_seq_equals_model"]) for r in up_rows)),
            "n_domain_hit_rows": len(dom_rows),
            "source_databases": dict(Counter(r["source_database"] for r in dom_rows))}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["uniprot", "interpro"])
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    for d in (UP_DIR, IP_DIR, RESULTS):
        d.mkdir(parents=True, exist_ok=True)
    seqs = read_candidates()
    accs = list(seqs)[: args.limit] if args.limit else list(seqs)
    print(f"{len(accs)} accessions; cache in {DATA}")
    failed = []
    if args.only in (None, "uniprot"):
        fetch_uniprot(accs, failed)
    if args.only in (None, "interpro"):
        fetch_interpro(accs, failed)
    summary = build_tables(accs, seqs)
    manifest = {"script": Path(__file__).name, "run_finished_utc": now(), "n_accessions": len(accs),
                "foldtrace_run": str(RUN), "uniprot_endpoint": UNIPROT_BATCH.split("?")[0],
                "interpro_endpoint": INTERPRO_ONE.split("{")[0], **summary,
                "mixed_releases_warning": len(summary["uniprot_releases"]) > 1 or len(summary["interpro_versions"]) > 1,
                "still_failed": [f"{s}:{a}" for s, a in failed]}
    save_json(RESULTS / "p2_download_manifest.json", manifest)
    print(json.dumps(manifest, indent=1))
    if failed:
        print(f"\n{len(failed)} requests still failed. Run the script again; cached ones are skipped.")
        sys.exit(1)
    print("\nDone. Send me the printout above (or results/p2_download_manifest.json).")


if __name__ == "__main__":
    main()
