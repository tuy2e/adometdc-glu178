"""Track B: does each fusion appear on a single transcript? tblastn of the p14 junction queries against NCBI TSA.

Queries: data/p14_tblastn_queries/*.txt (last 60 aa of the AdoMetDC Pfam match + linker + first 60 aa of the
PLP-decarboxylase match; segment lengths and organisms from that folder's README table).

NCBI BLAST URL API (https://blast.ncbi.nlm.nih.gov/Blast.cgi), CMD=Put / CMD=Get:
  PROGRAM=tblastn  DATABASE=tsa_nt  ENTREZ_QUERY="<organism>[ORGN]"  HITLIST_SIZE=100  FORMAT_TYPE=XML
One job at a time; status polled every POLL seconds (>= 10, NCBI's limit); no email, tool or personal
identifier is sent. The database actually searched is read back from each XML report (BlastOutput_db).
Fallback, only when a search returns no hits at all: genus[ORGN], then no organism restriction.

Per subject:
  positive        one HSP that starts inside the AdoMetDC segment and ends inside the ODC segment, i.e. a
                  single contiguous alignment through the whole linker; gaps inside the linker are counted
                  (positive_no_junction_gap only if that count is 0)
  both_sides_split  the AdoMetDC and ODC segments are both hit on this subject but by different HSPs
                  (not contiguous - not a positive)
  adometdc_only / odc_only / linker_only   partial - NOT evidence (each species has both enzymes separately)
For positives, the subject's organism and BioProject come from NCBI nuccore (GenBank DBLINK) and are compared
with the genome BioProjects for that species in results_81/q3_assemblies.tsv.

Writes data/p15_blast_raw/<query>__<tier>.xml (+ .rid), results/p15_tblastn.tsv, results/p15_tblastn.md
Each search has a hard limit of TIMEOUT (10 min); a search still WAITING then is recorded as timed out and the
run moves to the next query (no genus/no-organism fallback after a timeout - a timeout is not "no hits").
Progress (RID and the raw Status= line of every poll) is printed as it happens; do not pipe the output.
Run:  py -3.12 p15_tblastn.py            (resumes: a query/tier whose raw XML exists is not resubmitted)
      py -3.12 p15_tblastn.py --only 01  (one query)
"""
import csv
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

from p2_fetch_annotations import DATA, RESULTS, USER_AGENT

sys.stdout.reconfigure(encoding="utf8")
QDIR = DATA / "p14_tblastn_queries"
RAW = DATA / "p15_blast_raw"
RAW.mkdir(parents=True, exist_ok=True)
Q3 = RESULTS.parent / "results_81" / "q3_assemblies.tsv"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
DB, POLL, HITLIST, TIMEOUT = "tsa_nt", 15, 100, 600  # TIMEOUT: hard limit per search, seconds


class BlastTimeout(Exception):
    pass


class BlastError(Exception):
    pass
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None


def http(url, data=None, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf8", "replace")
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            print(f"    {type(e).__name__}: {e}; retry in {30 * (k + 1)} s", flush=True)
            time.sleep(30 * (k + 1))
    raise RuntimeError(f"gave up: {url}")


# ---------------------------------------------------------------- queries
def queries():
    readme = (QDIR / "README.md").read_text(encoding="utf8")
    out = []
    for m in re.finditer(r"^\| (\d+) \| `([^`]+)` \| \*([^*]+)\* \|.*\| (\d+) aa \| (\d+) / (\d+) / (\d+) \|$", readme, re.M):
        n, stem, org, L, a, l, o = m.groups()
        seq = "".join((QDIR / stem).read_text().split())
        assert len(seq) == int(L) == int(a) + int(l) + int(o), stem
        out.append({"n": n.zfill(2), "stem": stem.removesuffix(".txt"), "organism": org.replace("sp GSL018", "sp. GSL018"),
                    "seq": seq, "L": len(seq), "ado": (1, int(a)), "linker": (int(a) + 1, int(a) + int(l)),
                    "odc": (int(a) + int(l) + 1, len(seq))})
    return out


# ---------------------------------------------------------------- BLAST
def run_blast(q, tier, entrez):
    base = RAW / f"{q['stem']}__{tier}"
    xml = base.with_suffix(".xml")
    if xml.exists():
        return xml.read_text(encoding="utf8")
    form = {"CMD": "Put", "PROGRAM": "tblastn", "DATABASE": DB, "QUERY": q["seq"], "HITLIST_SIZE": HITLIST}
    if entrez:
        form["ENTREZ_QUERY"] = entrez
    body = http(BLAST, urllib.parse.urlencode(form).encode())
    rid = re.search(r"RID = (\S+)", body)
    rtoe = re.search(r"RTOE = (\d+)", body)
    if not rid:
        raise RuntimeError(f"no RID returned:\n{body[:2000]}")
    rid = rid.group(1)
    base.with_suffix(".rid").write_text(f"{rid}\t{tier}\t{entrez}\n")
    print(f"  {q['n']} {tier}: submitted, RID {rid}, NCBI estimate RTOE={rtoe.group(1) if rtoe else '?'} s, "
          f"entrez={entrez or '(none)'}", flush=True)
    t0 = time.time()
    time.sleep(POLL)
    while True:
        info = http(f"{BLAST}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}")
        raw = " ".join(re.findall(r"(?:Status|ThereAreHits)=\w+", info)) or "(no Status line in response)"
        st = re.search(r"Status=(\w+)", info)
        st = st.group(1) if st else "?"
        print(f"    RID {rid} t={int(time.time() - t0)} s: {raw}", flush=True)
        if st == "READY":
            break
        if st in ("FAILED", "UNKNOWN"):
            raise BlastError(f"{rid}, status {st}")
        if time.time() - t0 > TIMEOUT:
            raise BlastTimeout(rid)
        time.sleep(POLL)
    out = http(f"{BLAST}?CMD=Get&FORMAT_TYPE=XML&RID={rid}")
    if "<BlastOutput" not in out:
        raise RuntimeError(f"RID {rid}: no XML report:\n{out[:2000]}")
    xml.write_text(out, encoding="utf8")
    time.sleep(POLL)
    return out


def parse(xml_text):
    root = ET.fromstring(xml_text)
    db = root.findtext("BlastOutput_db")
    hits = []
    for h in root.iter("Hit"):
        hsps = []
        for p in h.iter("Hsp"):
            g = lambda k: p.findtext(f"Hsp_{k}")
            hsps.append({"qf": int(g("query-from")), "qt": int(g("query-to")), "hf": int(g("hit-from")),
                         "ht": int(g("hit-to")), "frame": g("hit-frame"), "e": float(g("evalue")),
                         "bits": float(g("bit-score")), "ident": int(g("identity")), "alen": int(g("align-len")),
                         "gaps": int(g("gaps") or 0), "qseq": g("qseq"), "hseq": g("hseq"), "mid": g("midline")})
        hits.append({"acc": h.findtext("Hit_accession"), "def": h.findtext("Hit_def"), "len": int(h.findtext("Hit_len")),
                     "hsps": hsps})
    return db, hits


def inside(x, seg):
    return seg[0] <= x <= seg[1]


def covers(h, seg):
    return h["qf"] <= seg[1] and h["qt"] >= seg[0]


def junction_view(q, h):
    """Alignment columns whose query position lies in the linker, plus 10 aa either side; gaps counted in linker."""
    lo, hi = q["linker"][0] - 10, q["linker"][1] + 10
    qa, ha, ma = [], [], []
    pos, gaps = h["qf"], 0
    for qc, hc, mc in zip(h["qseq"], h["hseq"], h["mid"]):
        cur = pos
        if lo <= cur <= hi:
            qa.append(qc)
            ha.append(hc)
            ma.append(mc)
            if inside(cur, q["linker"]) and (qc == "-" or hc == "-"):
                gaps += 1
        if qc != "-":
            pos += 1
    return gaps, f"Q {lo}-{hi}: {''.join(qa)} | M: {''.join(ma)} | S: {''.join(ha)}"


def classify(q, hit):
    full = [h for h in hit["hsps"] if inside(h["qf"], q["ado"]) and inside(h["qt"], q["odc"])]
    if full:
        return "positive", max(full, key=lambda h: h["bits"])
    a = any(covers(h, q["ado"]) for h in hit["hsps"])
    o = any(covers(h, q["odc"]) for h in hit["hsps"])
    best = max(hit["hsps"], key=lambda h: h["bits"])
    if a and o:
        return "both_sides_split", best
    if a:
        return "adometdc_only", best
    if o:
        return "odc_only", best
    return "linker_only", best


def species_of(defline):
    m = re.match(r"TSA: (\S+ \S+(?: sp\. \S+)?)", defline or "")
    return m.group(1) if m else ""


def nuccore_meta(acc):
    """Organism and BioProject from the GenBank record (E-utilities efetch, no identifiers sent)."""
    gb = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" +
              urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "gb", "retmode": "text"}))
    time.sleep(1)
    org = re.search(r"^\s+ORGANISM\s+(.+)$", gb, re.M)
    bp = re.search(r"BioProject:\s*(PRJ\w+)", gb)
    return (org.group(1).strip() if org else ""), (bp.group(1) if bp else "")


def genome_projects():
    d = defaultdict(set)
    with open(Q3, encoding="utf8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["study_accession"]:
                d[r["species"]].add(r["study_accession"])
    return d


# ---------------------------------------------------------------- main
def main():
    gp = genome_projects()
    rows = []
    for q in queries():
        if ONLY and q["n"] != ONLY.zfill(2):
            continue
        print(f"{q['n']} {q['organism']} ({q['L']} aa; AdoMetDC {q['ado']}, linker {q['linker']}, ODC {q['odc']})", flush=True)
        tiers = [("species", f'"{q["organism"]}"[ORGN]'), ("genus", f'"{q["organism"].split()[0]}"[ORGN]'), ("none", "")]
        failed = None
        for tier, entrez in tiers:
            try:
                db, hits = parse(run_blast(q, tier, entrez))
            except (BlastTimeout, BlastError) as e:
                failed = ("timed out after 10 min" if isinstance(e, BlastTimeout) else "NCBI error") + f" (RID {e})"
                print(f"  {q['n']} {tier}: {failed} - giving up on this query", flush=True)
                break
            if db != DB:
                raise RuntimeError(f"report says database {db!r}, expected {DB!r}")
            print(f"  {tier}: {len(hits)} subjects", flush=True)
            if hits:
                break
        if failed:
            rows.append({"query": q["stem"], "organism": q["organism"], "query_len": q["L"], "class": failed,
                         "search_tier_with_result": tier, "entrez_query": entrez})
            continue
        base = {"query": q["stem"], "organism": q["organism"], "query_len": q["L"],
                "segments_ado_linker_odc": f"{q['ado'][0]}-{q['ado'][1]} / {q['linker'][0]}-{q['linker'][1]} / {q['odc'][0]}-{q['odc'][1]}",
                "search_tier_with_result": tier if hits else "none (all three tiers empty)", "entrez_query": entrez,
                "database": db, "rid": (RAW / f"{q['stem']}__{tier}.rid").read_text().split("\t")[0]}
        if not hits:
            rows.append({**base, "class": "no hits"})
            continue
        for hit in hits:
            cls, h = classify(q, hit)
            r = {**base, "class": cls, "subject": hit["acc"], "subject_def": hit["def"],
                 "subject_species": species_of(hit["def"]), "subject_len": hit["len"],
                 "best_hsp_query": f"{h['qf']}-{h['qt']}", "query_coverage_pct": round(100 * (h["qt"] - h["qf"] + 1) / q["L"], 1),
                 "pct_identity": round(100 * h["ident"] / h["alen"], 1), "evalue": h["e"], "bits": h["bits"],
                 "frame": h["frame"], "n_hsps": len(hit["hsps"])}
            if cls == "positive":
                gaps, view = junction_view(q, h)
                org, bp = nuccore_meta(hit["acc"])
                same = bp and bp in gp.get(q["organism"], set())
                r.update({"subject_species": org or r["subject_species"], "junction_gaps": gaps,
                          "junction_alignment": view, "tsa_bioproject": bp,
                          "genome_bioprojects_q3": ";".join(sorted(gp.get(q["organism"], []))) or "none in q3",
                          "independent_of_genome": ("NO - same BioProject as the genome" if same else
                                                    "yes (different BioProject)" if bp else "unknown (no BioProject on record)")})
            rows.append(r)
    write(rows)


def write(rows):
    cols = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    name = "p15_tblastn" + (f"_{ONLY.zfill(2)}" if ONLY else "")
    with open(RESULTS / f"{name}.tsv", "w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    by = defaultdict(list)
    for r in rows:
        by[r["query"]].append(r)
    L = ["# p15: tblastn of the p14 junction queries against NCBI TSA", "",
         f"Generated by `scripts/p15_tblastn.py`. Database `{DB}` (confirmed from each report), tblastn, "
         f"hitlist {HITLIST}, NCBI default expect. Raw XML in `data/p15_blast_raw/`.", "",
         "**Positive** = one HSP running from the AdoMetDC segment through the whole linker into the ODC segment. "
         "Partial hits (one side only, or both sides by separate HSPs) are **not evidence**: every species has both "
         "enzymes separately. **No hits means nothing**: many of these species have little public RNA data.", "",
         "| query | organism | tier that gave a result | subjects | positives | both sides, split | AdoMetDC only | ODC only |",
         "|---|---|---|---|---|---|---|---|"]
    for qn, rs in by.items():
        c = defaultdict(int)
        for r in rs:
            c[r["class"]] += 1
        n = 0 if rs[0]["class"] == "no hits" or rs[0]["class"].startswith(("timed out", "NCBI error")) else len(rs)
        if n == 0 and rs[0]["class"] != "no hits":
            L.append(f"| {qn} | *{rs[0]['organism']}* | {rs[0]['class']} at tier {rs[0]['search_tier_with_result']} | - | - | - | - | - |")
            continue
        L.append(f"| {qn} | *{rs[0]['organism']}* | {rs[0]['search_tier_with_result']} | {n} | {c['positive']} | "
                 f"{c['both_sides_split']} | {c['adometdc_only']} | {c['odc_only']} |")
    L += ["", "## Positives", ""]
    pos = [r for r in rows if r["class"] == "positive"]
    if not pos:
        L.append("None.")
    for r in pos:
        L += [f"### {r['query']} -> {r['subject']} (*{r['subject_species']}*)", "",
              f"- identity {r['pct_identity']} %, query coverage {r['query_coverage_pct']} % ({r['best_hsp_query']} of "
              f"{r['query_len']}), e-value {r['evalue']:.2g}, gaps in linker {r['junction_gaps']}",
              f"- search tier: {r['search_tier_with_result']}",
              f"- TSA BioProject {r['tsa_bioproject'] or '?'}; genome BioProjects (q3) {r['genome_bioprojects_q3']}; "
              f"independent of the genome: {r['independent_of_genome']}",
              "", "```", *r["junction_alignment"].split(" | "), "```", ""]
    (RESULTS / f"{name}.md").write_text("\n".join(L), encoding="utf8")


if __name__ == "__main__":
    main()
