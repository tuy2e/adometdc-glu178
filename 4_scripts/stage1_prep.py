"""Stage 1: put every go/no-go structure into the 1MSV frame.

Same method as the 10 Nov docking round (4_scripts/docking_prep.py): map each chain to human
AdoMetDC by pairwise MAFFT, superimpose on the Calpha of human residues within 18 A of the
crystal putrescine, then refine the core twice to atoms under 2 A.

Differences from docking_prep.py, needed because Stage 1 includes crystals and multi-chain
structures (see the clarification in 1_reports/PREREG_pocket_chemistry_measurements.md):
  * every chain is tested against human; chains that map (>=20% identity, >=60 mapped residues)
    are kept in the output, so a prozyme or a processed beta chain stays in the structure
  * the superposition itself uses only the best-matching chain
  * waters, ions, ligands and the bound putrescine are stripped

Usage:  python stage1_prep.py <structure_dir> <out_dir>
Writes <out_dir>/<name>_aligned.pdb (protein only) and <out_dir>/prep_stage1.json
"""
import glob
import json
import os
import subprocess
import sys
import warnings

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer

warnings.filterwarnings("ignore")

# Human AdoMetDC, UniProt P17707, as used throughout this project.
HUMAN = ("MEAAHFFEGTEKLLEVWFSRQQPDANQGSGDLRTIPRSEWDILLKDVQCSIISVTKTDKQEAYVLSESSMFVSKRRFILKT"
         "CGTTLLLKALVPLLKLARDYSGFDSIQSFFYSRKNFMKPSHQGYPHRNFQEEIEFLNAIFPNGAAYCMGRMNSDCWYLYTL"
         "DFPESRVISQPDQTLEILMSELDPAVMDQFYMKDGVTAKDVTRESGIRDLIPGSVIDATMFNPCGYSMNGMKSDGTYWTIH"
         "ITPEPEFSYVSFETNLSQTSYDDLIRKVVEVFKPGKFVTTLFVNQSSKCRTVLASPQKIEGFKRLDCQSAMFNDYNFVFTS"
         "FAKKQQQQQS")

PUT_CENTRE = np.array([61.324, 11.189, 7.513])
POS = [13, 15, 80, 109, 111, 113, 174, 176, 178, 256, 285]

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
         # pyruvoyl group left by self-processing, and common modified residues
         "PYR": "S", "MSE": "M", "SEP": "S", "PTR": "Y"}

parser = PDBParser(QUIET=True)


class ProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.get_resname() in THREE and residue.id[0] in (" ", "H_MSE")

    def accept_atom(self, atom):
        return atom.element != "H"


def chain_residues(chain):
    return [r for r in chain if r.get_resname() in THREE and "CA" in r]


def mafft_map(target_seq, workdir):
    """human position (1-based) -> index into target residue list (1-based)."""
    fa = os.path.join(workdir, "_pair.fa")
    with open(fa, "w") as fh:
        fh.write(">H\n%s\n>T\n%s\n" % (HUMAN, target_seq))
    out = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", fa],
                         capture_output=True, text=True).stdout
    seqs, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    if "H" not in seqs or "T" not in seqs:
        return {}, 0.0
    mapping, hc, tc, ident, naln = {}, 0, 0, 0, 0
    for x, y in zip(seqs["H"], seqs["T"]):
        if x != "-":
            hc += 1
        if y != "-":
            tc += 1
        if x != "-" and y != "-":
            mapping[hc] = tc
            naln += 1
            if x == y:
                ident += 1
    return mapping, (ident / naln if naln else 0.0)


def main(indir, outdir):
    os.makedirs(outdir, exist_ok=True)
    ref = parser.get_structure("ref", os.path.join(indir, "human_1MSV_aligned.pdb"))[0]
    href = {}
    for ch in ref:
        for r in chain_residues(ch):
            href.setdefault(r.id[1], r)
    # human residues lining the site, used as the superposition target
    site_pos = [n for n, r in href.items()
                if np.min(np.linalg.norm(np.array([a.coord for a in r]) - PUT_CENTRE, axis=1)) < 18]

    results = {}
    targets = sorted(glob.glob(os.path.join(indir, "*.pdb")))
    for path in targets:
        name = os.path.basename(path)[:-4]
        if name.endswith("_aligned"):
            continue  # already in the frame, from the 10 Nov round
        struct = parser.get_structure(name, path)[0]

        # which chains look like AdoMetDC?
        per_chain = {}
        for ch in struct:
            res = chain_residues(ch)
            if len(res) < 40:
                continue
            seq = "".join(THREE[r.get_resname()] for r in res)
            mapping, ident = mafft_map(seq, outdir)
            nmapped = len(mapping)
            if ident >= 0.20 and nmapped >= 60:
                per_chain[ch.id] = dict(res=res, map=mapping, ident=ident, n=nmapped)
        if not per_chain:
            print("  %-12s NO CHAIN MAPPED TO HUMAN - skipped" % name)
            results[name] = dict(error="no chain mapped")
            continue

        best_id = max(per_chain, key=lambda c: (per_chain[c]["ident"], per_chain[c]["n"]))
        best = per_chain[best_id]
        tby = {i + 1: r for i, r in enumerate(best["res"])}
        mp = best["map"]

        pairs = [n for n in site_pos if n in mp and mp[n] in tby]
        if len(pairs) < 25:
            print("  %-12s only %d site residues mapped - skipped" % (name, len(pairs)))
            results[name] = dict(error="too few site residues mapped (%d)" % len(pairs))
            continue

        sup = Superimposer()
        sup.set_atoms([href[n]["CA"] for n in pairs], [tby[mp[n]]["CA"] for n in pairs])
        core = pairs
        for _ in range(2):
            rot, tran = sup.rotran
            d = [np.linalg.norm(href[n]["CA"].coord - (np.dot(tby[mp[n]]["CA"].coord, rot) + tran))
                 for n in core]
            tight = [n for n, x in zip(core, d) if x < 2.0]
            if len(tight) >= 20:
                core = tight
                sup = Superimposer()
                sup.set_atoms([href[n]["CA"] for n in core], [tby[mp[n]]["CA"] for n in core])
        rot, tran = sup.rotran
        for atom in struct.get_atoms():
            atom.transform(rot, tran)

        io = PDBIO()
        io.set_structure(struct)
        io.save(os.path.join(outdir, "%s_aligned.pdb" % name), ProteinOnly())

        res_at = {}
        for p in POS:
            if p in mp and mp[p] in tby:
                r = tby[mp[p]]
                res_at[p] = "%s%d" % (THREE.get(r.get_resname(), "X"), r.id[1])
            else:
                res_at[p] = "-"
        # AlphaFold models carry pLDDT in the B-factor column; crystals carry a real B-factor
        bvals = [tby[mp[p]]["CA"].bfactor for p in POS if p in mp and mp[p] in tby]
        results[name] = dict(
            chain_used=best_id,
            chains_kept=sorted(per_chain),
            identity=round(best["ident"], 3),
            core_ca=len(core),
            rmsd=round(float(sup.rms), 2),
            site_bfactor=round(float(np.mean(bvals)), 1) if bvals else None,
            residues=res_at,
        )
        print("  %-12s chain %s (id %.0f%%, %d chains kept) | core %3d CA, RMSD %.2f | site B/pLDDT %5.1f | %s"
              % (name, best_id, 100 * best["ident"], len(per_chain), len(core), sup.rms,
                 results[name]["site_bfactor"] or -1,
                 " ".join("%d:%s" % (k, v) for k, v in res_at.items())))

    with open(os.path.join(outdir, "prep_stage1.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    print("\nWrote %s" % os.path.join(outdir, "prep_stage1.json"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
