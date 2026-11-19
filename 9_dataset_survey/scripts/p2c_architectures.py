"""Phase 2c: merge the three annotation sources and take a census of domain architectures. Census only.

Annotation source per protein (one source per protein, never mixed):
    interpro            active UniProt accession: InterPro API matches (results/domain_hits.tsv)
    replacement_record  deleted accession with an active same-sequence twin: the twin's InterPro matches
                        (data/interpro_successors/), positions valid because the sequence is identical
    uniparc_features    other deleted accessions: UniParc sequence features (results/domain_hits_uniparc.tsv)
    none                nothing available
EXCLUDED from everything position-based: A0A8C5BDJ7 (UniProt sequence 227 aa, AlphaFold model 245 aa).

Architecture unit = Pfam, because it is the one database reported with the same accessions by both
InterPro and UniParc. A Pfam match counts as a confirmed domain. Regions >= 40 aa not covered by any
Pfam domain but matched by a structural-domain database (CATH-Gene3D, SUPERFAMILY, CDD, SMART,
PROSITE profiles) are listed separately as uncertain segments ([...] in arch_extended).
AdoMetDC Pfam families: PF01536 (Adenosylmethionine decarboxylase) and PF02675 (bacterial AdoMet_dc).

Every count is also given for the InterPro-only set (source = interpro), so a pattern that appears
only in the UniParc-derived or replacement data is visible.

Writes
    results/domain_hits_merged.tsv          every match from all three sources, with a source column
    results/architecture_per_protein.tsv    one row per candidate: domains in order, positions, linkers
    results/fusion_architectures.tsv        one row per distinct Pfam architecture, with counts
    results/architecture_census_summary.txt
Run (seconds; fetches names for a few Pfam accessions not seen in InterPro data, cached):
    py -3.12 p2c_architectures.py
"""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from p2_fetch_annotations import DATA, RESULTS, get, save_json

EXCLUDE_POSITIONS = {"A0A8C5BDJ7": "UniProt sequence 227 aa != AlphaFold model 245 aa"}
ADOMETDC = {"PF01536", "PF02675"}
STRUCTURAL = {"cathgene3d", "ssf", "cdd", "smart", "profile"}
MIN_SEGMENT = 40
RARE = 5
UNIPARC_DB = {"Pfam": "pfam", "Gene3D": "cathgene3d", "SUPFAM": "ssf", "PANTHER": "panther", "NCBIfam": "ncbifam",
              "PROSITE": "prosite", "PIRSF": "pirsf", "PIRSR": "pirsr", "FunFam": "funfam", "CDD": "cdd",
              "SMART": "smart", "PRINTS": "prints", "HAMAP": "hamap"}


def read(path):
    with open(path, encoding="utf8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(name, rows, fields=None):
    with open(RESULTS / name, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def interpro_rows(acc, results, source, src_acc, version, date):
    out = []
    for r in results:
        md = r["metadata"]
        for prot in r.get("proteins", []):
            for loc in prot.get("entry_protein_locations") or []:
                for fr in loc.get("fragments", []):
                    out.append({"accession": acc, "source": source, "source_accession": src_acc,
                                "database": md["source_database"], "entry_accession": md["accession"],
                                "entry_name": md.get("name", ""), "entry_type": md.get("type", ""),
                                "start": int(fr["start"]), "end": int(fr["end"]),
                                "version": f"InterPro {version}", "downloaded_utc": date})
    return out


def main():
    master = {r["accession"]: r for r in read(RESULTS / "master_candidates.tsv")}
    annot = {r["accession"]: r for r in read(RESULTS / "annotation_uniprot.tsv")}
    inact = {r["accession"]: r for r in read(RESULTS / "inactive_resolution.tsv")}
    ip_hits = defaultdict(list)
    for r in read(RESULTS / "domain_hits.tsv"):
        ip_hits[r["accession"]].append(r)
    up_hits = defaultdict(list)
    for r in read(RESULTS / "domain_hits_uniparc.tsv"):
        up_hits[r["accession"]].append(r)

    merged, source = [], {}
    for a, m in master.items():
        if annot[a]["uniprot_status"] == "active":
            source[a] = "interpro"
            for r in ip_hits[a]:
                merged.append({"accession": a, "source": "interpro", "source_accession": a,
                               "database": r["source_database"], "entry_accession": r["entry_accession"],
                               "entry_name": r["entry_name"], "entry_type": r["entry_type"],
                               "start": int(r["start"]), "end": int(r["end"]),
                               "version": f"InterPro {r['interpro_version']}", "downloaded_utc": r["interpro_downloaded_utc"]})
            continue
        rep = inact[a]
        if rep["replacement_interpro"] == "ok":
            source[a] = "replacement_record"
            c = json.loads((DATA / "interpro_successors" / f"{rep['replacement_used']}.json").read_text(encoding="utf8"))
            merged += interpro_rows(a, c["results"], "replacement_record", rep["replacement_used"],
                                    (c["versions"] or {}).get("interpro_version", ""), c["downloaded_utc"])
        elif up_hits[a]:
            source[a] = "uniparc_features"
            for r in up_hits[a]:
                merged.append({"accession": a, "source": "uniparc_features", "source_accession": r["uniparc_id"],
                               "database": UNIPARC_DB.get(r["source_database"], r["source_database"].lower()),
                               "entry_accession": r["entry_accession"], "entry_name": "", "entry_type": "",
                               "start": int(r["start"]), "end": int(r["end"]),
                               "version": f"UniProt {r['uniparc_release']}", "downloaded_utc": rep["downloaded_utc"]})
        else:
            source[a] = "none"
    for r in merged:
        r["position_based_excluded"] = EXCLUDE_POSITIONS.get(r["accession"], "")

    # Pfam names: from InterPro rows, else from the InterPro API (cached)
    names_p = DATA / "pfam_names.json"
    names = json.loads(names_p.read_text(encoding="utf8")) if names_p.exists() else {}
    for r in merged:
        if r["database"] == "pfam" and r["entry_name"]:
            names.setdefault(r["entry_accession"], r["entry_name"])
    for pf in sorted({r["entry_accession"] for r in merged if r["database"] == "pfam"} - set(names)):
        st, _, body = get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}/")
        md = json.loads(body)["metadata"] if st == 200 and body.strip() else {}
        nm = md.get("name")
        names[pf] = (nm.get("name") if isinstance(nm, dict) else nm) or "(name not found)"
    save_json(names_p, names)
    for r in merged:
        if r["database"] == "pfam" and not r["entry_name"]:
            r["entry_name"] = names[r["entry_accession"]]
    write("domain_hits_merged.tsv", merged, ["accession", "source", "source_accession", "database", "entry_accession",
                                             "entry_name", "entry_type", "start", "end", "version", "downloaded_utc",
                                             "position_based_excluded"])

    by_acc = defaultdict(list)
    for r in merged:
        by_acc[r["accession"]].append(r)

    # ---------- per protein
    prot_rows = []
    for a, m in master.items():
        L = int(m["length"])
        row = {"accession": a, "fold_ok": m["fold_ok"], "source": source[a], "organism": m["organism"],
               "length": L, "uniprot_status": annot[a]["uniprot_status"],
               "position_based_excluded": EXCLUDE_POSITIONS.get(a, "")}
        if a in EXCLUDE_POSITIONS or source[a] == "none":
            row.update(arch_pfam_ids="(excluded)" if a in EXCLUDE_POSITIONS else "(no annotation)",
                       arch_pfam_names="", n_pfam_domains="", adometdc_order="", partner_pfam="",
                       domains_with_positions="", linkers="", arch_extended="", uncertain_segments="")
            prot_rows.append(row)
            continue
        pf = sorted((r for r in by_acc[a] if r["database"] == "pfam"), key=lambda r: (r["start"], r["end"]))
        ids = [r["entry_accession"] for r in pf]
        # linkers: N-terminal, between consecutive domains (negative = overlap), C-terminal
        if pf:
            link = [f"N:{pf[0]['start'] - 1}"]
            link += [f"{ids[i]}>{ids[i + 1]}:{pf[i + 1]['start'] - pf[i]['end'] - 1}" for i in range(len(pf) - 1)]
            link.append(f"C:{L - pf[-1]['end']}")
        else:
            link = []
        # uncertain segments: >= MIN_SEGMENT aa outside Pfam, covered by a structural-domain database
        covered = [False] * (L + 2)
        for r in pf:
            for k in range(r["start"], min(r["end"], L) + 1):
                covered[k] = True
        segs = []
        for r in by_acc[a]:
            if r["database"] not in STRUCTURAL:
                continue
            free = sum(1 for k in range(r["start"], min(r["end"], L) + 1) if not covered[k])
            if free >= MIN_SEGMENT:
                segs.append((r["start"], r["end"], f"{r['database']}:{r['entry_accession']}"))
        merged_segs = []
        for s, e, lab in sorted(segs):
            if merged_segs and s <= merged_segs[-1][1]:
                merged_segs[-1][1] = max(merged_segs[-1][1], e)
                merged_segs[-1][2].add(lab)
            else:
                merged_segs.append([s, e, {lab}])
        items = [(r["start"], r["entry_accession"]) for r in pf] + [(s, f"[uncertain {s}-{e}]") for s, e, _ in merged_segs]
        ad = [i for i, x in enumerate(ids) if x in ADOMETDC]
        partners = [x for x in ids if x not in ADOMETDC]
        if not ids:
            order = "no Pfam domain"
        elif not ad:
            order = "no AdoMetDC Pfam"
        elif not partners:
            order = "AdoMetDC only" if len(ad) == 1 else "AdoMetDC repeated"
        elif ad[0] == 0 and ad[-1] == len(ids) - 1:
            order = "AdoMetDC at both ends"
        elif ad[0] == 0:
            order = "AdoMetDC first"
        elif ad[-1] == len(ids) - 1:
            order = "AdoMetDC last"
        else:
            order = "AdoMetDC in the middle"
        row.update(arch_pfam_ids=" > ".join(ids) if ids else "(no Pfam domain)",
                   arch_pfam_names=" > ".join(names.get(x, x) for x in ids),
                   n_pfam_domains=len(ids), adometdc_order=order,
                   partner_pfam="; ".join(sorted({f"{x} {names.get(x, '')}" for x in partners})),
                   domains_with_positions="; ".join(f"{r['entry_accession']} {r['start']}-{r['end']}" for r in pf),
                   linkers=" ".join(link),
                   arch_extended=" > ".join(x for _, x in sorted(items)) if items else "(none)",
                   uncertain_segments="; ".join(f"{s}-{e} ({', '.join(sorted(l))})" for s, e, l in merged_segs))
        prot_rows.append(row)

    census = [r for r in prot_rows if r["arch_pfam_ids"] not in ("(excluded)", "(no annotation)")]
    n_arch = Counter(r["arch_pfam_ids"] for r in census)
    n_ext = Counter(r["arch_extended"] for r in census)
    for r in prot_rows:
        r["n_proteins_same_arch"] = n_arch.get(r["arch_pfam_ids"], "")
        r["unique_pfam_layout"] = n_arch.get(r["arch_pfam_ids"]) == 1
        r["unique_extended_layout"] = n_ext.get(r["arch_extended"]) == 1 if r in census else ""
    write("architecture_per_protein.tsv", prot_rows,
          ["accession", "fold_ok", "source", "uniprot_status", "organism", "length", "arch_pfam_ids", "arch_pfam_names",
           "n_pfam_domains", "adometdc_order", "partner_pfam", "domains_with_positions", "linkers",
           "uncertain_segments", "arch_extended", "n_proteins_same_arch", "unique_pfam_layout",
           "unique_extended_layout", "position_based_excluded"])

    # ---------- per architecture
    groups = defaultdict(list)
    for r in census:
        groups[r["arch_pfam_ids"]].append(r)
    arch_rows = []
    for k, rs in groups.items():
        src = Counter(r["source"] for r in rs)
        lens = sorted(r["length"] for r in rs)
        n_ip = src.get("interpro", 0)
        arch_rows.append({
            "arch_pfam_ids": k, "arch_pfam_names": rs[0]["arch_pfam_names"], "adometdc_order": rs[0]["adometdc_order"],
            "partner_pfam": rs[0]["partner_pfam"],
            "n_proteins": len(rs), "n_fold_ok": sum(r["fold_ok"] == "true" for r in rs),
            "n_interpro": n_ip, "n_replacement_record": src.get("replacement_record", 0),
            "n_uniparc_features": src.get("uniparc_features", 0),
            "n_fold_ok_interpro": sum(r["fold_ok"] == "true" and r["source"] == "interpro" for r in rs),
            "seen_in_interpro_set": n_ip > 0,
            "only_in_non_interpro_sources": n_ip == 0,
            "rare_le5": len(rs) <= RARE, "unique": len(rs) == 1,
            "length_min": lens[0], "length_median": lens[len(lens) // 2], "length_max": lens[-1],
            "organisms": "; ".join(f"{r['organism']} ({r['accession']}, {r['source']}, fold_ok={r['fold_ok']})"
                                   for r in rs) if len(rs) <= RARE else
                         "; ".join(f"{o}" for o, _ in Counter(r["organism"] for r in rs).most_common(5)) + " ...",
        })
    arch_rows.sort(key=lambda r: -r["n_proteins"])
    write("fusion_architectures.tsv", arch_rows)

    # ---------- summary
    def summarise(rows, label):
        a = Counter(r["arch_pfam_ids"] for r in rows)
        o = Counter(r["adometdc_order"] for r in rows)
        rare = [k for k, n in a.items() if n <= RARE]
        uniq = [k for k, n in a.items() if n == 1]
        return [f"== {label}: {len(rows)} proteins, {len(a)} distinct Pfam architectures, "
                f"{len(rare)} rare (<= {RARE} proteins), {len(uniq)} unique (1 protein)",
                "   AdoMetDC position: " + ", ".join(f"{k} {v}" for k, v in o.most_common())]
    fo = [r for r in census if r["fold_ok"] == "true"]
    lines = [f"sources (all 1,000): {dict(Counter(source.values()))}",
             f"excluded from position-based work: {', '.join(EXCLUDE_POSITIONS)}",
             f"no annotation: {sum(1 for r in prot_rows if r['arch_pfam_ids'] == '(no annotation)')} "
             f"(fold_ok: {sum(1 for r in prot_rows if r['arch_pfam_ids'] == '(no annotation)' and r['fold_ok'] == 'true')})"]
    lines += summarise(census, "all annotated")
    lines += summarise([r for r in census if r["source"] == "interpro"], "InterPro-only set")
    lines += summarise(fo, "fold-passers")
    lines += summarise([r for r in fo if r["source"] == "interpro"], "fold-passers, InterPro-only")
    lines.append("\nArchitectures (all annotated proteins; n / fold_ok / of which InterPro-source):")
    for r in arch_rows:
        flag = "  <- only in UniParc/replacement data" if r["only_in_non_interpro_sources"] else ""
        lines.append(f"  {r['n_proteins']:>4} / {r['n_fold_ok']:>3} / {r['n_interpro']:>3}  {r['arch_pfam_ids']}"
                     f"  [{r['adometdc_order']}]{flag}")
    txt = "\n".join(lines) + "\n"
    (RESULTS / "architecture_census_summary.txt").write_text(txt, encoding="utf8")
    print(txt)


if __name__ == "__main__":
    main()
