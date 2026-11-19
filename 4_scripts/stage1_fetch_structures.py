"""Stage 1 of the pocket-chemistry survey: download the go/no-go structures.

Run this on your own PC (Python 3, no extra packages needed):

    python stage1_fetch_structures.py

It writes into  <repo>/3_data/structures_stage1/  next to this script's parent folder,
and prints a manifest table at the end. Nothing is overwritten: a file that is already
there and non-empty is skipped, so you can re-run it safely after a dropped connection.

Sources: RCSB (https://files.rcsb.org) for crystals, AlphaFold DB
(https://alphafold.ebi.ac.uk) for models. Both are public and need no login.

Why you and not Claude: the cloud workspace Claude runs in cannot reach RCSB,
NCBI, UniProt or AlphaFold DB — the proxy refuses those hosts.
"""

import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "3_data", "structures_stage1")

# --- what to fetch -------------------------------------------------------------
# The list is fixed by 1_reports/PREREG_pocket_chemistry_measurements.md. Do not add
# stramenopile / cryptophyte / oomycete structures here; those belong to Stage 2.

CRYSTALS = [
    ("1MSV", "human S68A proenzyme + putrescine (reference frame)"),
    ("1JL0", "human H243A proenzyme + putrescine"),
    ("1I7B", "human processed enzyme + putrescine"),
    ("1JEN", "human, APO (no ligand) - separates lid closure from crystal-vs-model"),
    ("1MHM", "potato, not stimulated (Bennett 2002)"),
    ("5TVO", "T. brucei AdoMetDC monomer (Volkov 2016)"),
    ("5TVF", "T. brucei AdoMetDC/prozyme heterodimer"),
]

MODELS = [
    ("P17707", "human - crystal-vs-model control"),
    ("P21182", "Saccharomyces cerevisiae"),
    ("Q1K889", "Neurospora crassa"),
    ("O02655", "Caenorhabditis elegans"),
    ("Q27883", "Onchocerca volvulus"),
    ("O76240", "Trypanosoma cruzi, full length"),
    # Q587A7 (T. brucei) is NOT fetched: it is already in the read-only FOLDTRACE
    # model folder, so it is read from there instead of downloaded again.
    ("Q04694", "potato - crystal-vs-model control"),
    ("Q9U8D4", "Plasmodium falciparum AdoMetDC-ODC fusion (large file)"),
]

AF_VERSIONS = [6, 5, 4, 3, 2]  # newest first; the first one that exists is kept
TIMEOUT = 120
RETRIES = 3


def fetch(url):
    """Return the bytes at url, or None on a 404. Retries other failures."""
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "adometdc-stage1/1.0"}
            )
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:  # network blip, DNS, timeout
            last = e
        time.sleep(2 * (attempt + 1))
    print("    ! gave up: %s" % last)
    return None


def save(name, data):
    path = os.path.join(OUT, name)
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def already(name):
    path = os.path.join(OUT, name)
    return os.path.exists(path) and os.path.getsize(path) > 1000


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = []
    missing = []

    print("Writing to: %s\n" % OUT)

    print("Crystals from RCSB")
    for pdb_id, note in CRYSTALS:
        name = "%s.pdb" % pdb_id
        if already(name):
            print("  %-6s already there, skipped" % pdb_id)
            manifest.append((name, os.path.getsize(os.path.join(OUT, name)), "cached", note))
            continue
        print("  %-6s %s" % (pdb_id, note))
        data = fetch("https://files.rcsb.org/download/%s.pdb" % pdb_id)
        if data is None:
            print("    ! not retrieved")
            missing.append(pdb_id)
            continue
        save(name, data)
        manifest.append((name, len(data), "rcsb", note))

    print("\nModels from AlphaFold DB")
    for acc, note in MODELS:
        name = "AF-%s.pdb" % acc
        if already(name):
            print("  %-8s already there, skipped" % acc)
            manifest.append((name, os.path.getsize(os.path.join(OUT, name)), "cached", note))
            continue
        print("  %-8s %s" % (acc, note))
        got = None
        for v in AF_VERSIONS:
            url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.pdb" % (acc, v)
            data = fetch(url)
            if data:
                got = (data, v)
                break
        if got is None:
            print("    ! no model found at any version (tried v%s)"
                  % ", v".join(str(v) for v in AF_VERSIONS))
            missing.append(acc)
            continue
        data, v = got
        save(name, data)
        manifest.append((name, len(data), "alphafold v%d" % v, note))

    # --- manifest ---------------------------------------------------------------
    lines = ["file\tbytes\tsha256(first16)\tsource\tnote"]
    for name, size, src, note in manifest:
        with open(os.path.join(OUT, name), "rb") as fh:
            h = hashlib.sha256(fh.read()).hexdigest()[:16]
        lines.append("%s\t%d\t%s\t%s\t%s" % (name, size, h, src, note))
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "MANIFEST.tsv"), "w") as fh:
        fh.write(text)

    print("\n" + "-" * 72)
    print(text)
    print("%d files in %s" % (len(manifest), OUT))
    if missing:
        print("\nNOT retrieved: %s" % ", ".join(missing))
        print("Tell Claude which ones failed; they can be worked around.")
        return 1
    print("\nAll structures retrieved. Tell Claude it's done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
