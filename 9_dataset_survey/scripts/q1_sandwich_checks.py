"""Job 2 (handoff 19 Nov): is the HMG-box > AdoMetDC > VTT protein A0A8S1GM37 (Caenorhabditis auriculariae) a real gene?

Two checks, both offline. The inputs are pasted below from the public records (retrieved 19 Nov 2026), so the
script needs no network and can be re-run to reproduce the tables:

1. Exon structure. The ENA CDS record CAD6184317.1 (locus CAUJ_LOCUS236, contig CAJGYM010000001.1, assembly
   GCA_904845305, gene model "source:AUGUSTUS") gives a 23-exon join on the plus strand. Each intron is placed
   in protein coordinates (cumulative CDS nt / 3) and labelled with the domain or linker it falls in
   (Pfam positions from results/architecture_per_protein.tsv).
2. Is a normal-length AdoMetDC gene sitting inside the model? C. elegans SMD-1 (UniProt O02655, 368 aa,
   reviewed) is aligned locally to the full 865 aa protein (BLOSUM62, gap open 11, extend 1).

The five FOLDTRACE site calls for this protein are not recomputed: they are read from master_candidates.tsv.

Writes results_81/q1_sandwich_exons.tsv and prints the SMD-1 alignment summary.
Run (seconds):  py -3.12 q1_sandwich_checks.py
"""
import csv
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "results_81"

ACC = "A0A8S1GM37"
CDS_JOIN = """join(CAJGYM010000001.1:1351782..1351906,CAJGYM010000001.1:1351950..1352085,
CAJGYM010000001.1:1352618..1352692,CAJGYM010000001.1:1353818..1353911,CAJGYM010000001.1:1353954..1354048,
CAJGYM010000001.1:1354625..1354702,CAJGYM010000001.1:1354746..1354781,CAJGYM010000001.1:1354821..1354885,
CAJGYM010000001.1:1354931..1355176,CAJGYM010000001.1:1355216..1355336,CAJGYM010000001.1:1355382..1355585,
CAJGYM010000001.1:1355630..1355725,CAJGYM010000001.1:1355834..1355932,CAJGYM010000001.1:1355972..1356090,
CAJGYM010000001.1:1356138..1356267,CAJGYM010000001.1:1356313..1356391,CAJGYM010000001.1:1357848..1357978,
CAJGYM010000001.1:1358380..1358453,CAJGYM010000001.1:1358498..1358585,CAJGYM010000001.1:1359250..1359411,
CAJGYM010000001.1:1359453..1359591,CAJGYM010000001.1:1359654..1359706,CAJGYM010000001.1:1360240..1360392)"""

SEQ = ("MEEFDIYRIPPESPDWYLIEPTFFCSQTQKTAKEKPKEEVSRGVSAYALFFRERQLAEKQRTPRASFGELSKRIARQWDSLSQQQKMNYKRRCEMSRSLDMTKA"
       "VEERARRLETPRSGESLEGKPYRTDLWVLCALKTATANLGSQHSQNPDGASAAFTRRVATSPSSIGRLRRSRQPIQSKVNAIEATRRRPTFSIPLIEDCMGPDTKTRILPSPNVKSFK"
       "NQLHAIFLQKPKKMSVGVVSDCTQTLQELENVPIEEYFFEGAEKLLEIWFTSEQGALTKGLRNIPRSEIDRMLEVACCKVLHAKSNELIDSYVLSESSLFVSSHRIIIKTCGATRLLA"
       "ALPVLINLASKYGQLTQVQSVYYSRKNFLRPMLQPNVHRNFDSEVAYLDSFFEDGHAYCLGSLKQDRWYLYTMQREMPHPVAKRPDHTLEILMSELDDDVLRKFNRDVCIDGPDCFVR"
       "SGIDKIIPAGAEVHDELFDPCGYSMNAYMNETDQYATIHVTPEKEFSFASFETNQDLVCLRQQTLKVLKCFRPGKVLMTVFANDVSLKGKDAQQHLWERELPGYRRTNVQFLRLETET"
       "LVYANFIRKSGCDTSSDEDEGERKESTYRGSTEFWKKQMRELVCLVIIFILSSSALYLIWVLGPLPEGGKRPSFPRDLAGLQQLARNLNSYKNDHYTYTMVLFCSAYIFKQTFAIPGS"
       "FFMNLLAGALCGPWLGVLLVTTLTCFGASFCFILSSICARPLVERCLASRLAVLRAKMFGQRHRLWLFLLSARVFPFTPHWLLNVMSPALDVPLATHAFTVSIGMIPYNTLCVRAGGV"
       "LSKVKSVSDVLDTWTLAELLFVSVVLLTIARFSKRISQPRGISYLPTEHASVD")
SMD1_CELEGANS = ("MSATSATNFAVQTHPVKAPDEEYFFEGAEKLLELWFCSSTQNETRSLRIIPREEIDAMLDIARCKILHSKHNESIDSYVLSESSLFISDNRVILKTCGTTRLLAALPVIMQLAGAYAGLD"
                 "QVQSVYYSRKNFLRPDLQPSLHKNFDAEVEYLDSFFVDGHAYCLGSLKQDRWYLYTFHREVEFPAHKQPDHTLEILMSDLDEEVLHKFTKDYAVDGNDCFMRAGIDKIIPAGADVHDELF"
                 "DPCGYSMNAYMNDTDQYATIHVTPEKAFSFASFETNQDLVCLYSQTRKVLQCFRPNKILMTVFANDISEKGKDAQQQLWDRELPGYRRTNVQFVRLETETLVYAHFVRKASTGQDSSSSD"
                 "EDDGERSD")

# Pfam positions (results/architecture_per_protein.tsv, InterPro 111.0)
REGIONS = [("N-terminus", 1, 42), ("HMG-box PF09011", 43, 96), ("linker HMG>AdoMetDC", 97, 258),
           ("AdoMetDC PF01536", 259, 582), ("linker AdoMetDC>VTT", 583, 690), ("VTT PF09335", 691, 809),
           ("C-terminus", 810, 865)]


def region(pos):
    for name, a, b in REGIONS:
        if a <= pos <= b:
            return name
    return "?"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert len(SEQ) == 865
    exons = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", CDS_JOIN)]
    cds_nt = sum(b - a + 1 for a, b in exons)
    assert cds_nt == 3 * (len(SEQ) + 1), cds_nt          # 865 codons + stop
    rows, cum = [], 0
    for i, (a, b) in enumerate(exons):
        start_aa = cum // 3 + 1
        cum += b - a + 1
        intron = exons[i + 1][0] - b - 1 if i + 1 < len(exons) else ""
        after = round(cum / 3, 1) if intron != "" else ""
        rows.append({"exon": i + 1, "genomic_start": a, "genomic_end": b, "exon_nt": b - a + 1,
                     "protein_start": start_aa, "protein_end": min(len(SEQ), cum // 3),
                     "intron_after_nt": intron, "intron_after_residue": after,
                     "intron_falls_in": region(int(after)) if after != "" else ""})
    with open(OUT / "q1_sandwich_exons.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(exons)} exons, CDS {cds_nt} nt; wrote {OUT / 'q1_sandwich_exons.tsv'}")
    for r in rows:
        if r["intron_after_nt"] != "" and r["intron_after_nt"] >= 1000:
            print(f"  intron of {r['intron_after_nt']} nt after residue {r['intron_after_residue']} ({r['intron_falls_in']})")
    m = SEQ.find("KKMSVG") + 3
    print(f"in-frame Met at residue {m}; exon starting there: "
          f"{[r['exon'] for r in rows if r['protein_start'] in (m - 1, m)]}")

    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError:
        print("Biopython not installed: SMD-1 alignment skipped")
        return
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    a = al.align(SEQ, SMD1_CELEGANS)[0]
    t, q = a.aligned
    n = sum(te - ts for ts, te in t)
    ident = sum(SEQ[ts + k] == SMD1_CELEGANS[qs + k] for (ts, te), (qs, qe) in zip(t, q) for k in range(te - ts))
    print(f"C. elegans SMD-1 (O02655, {len(SMD1_CELEGANS)} aa) residues {q[0][0] + 1}-{q[-1][1]} align to "
          f"{ACC} residues {t[0][0] + 1}-{t[-1][1]}: {ident}/{n} identical ({100 * ident / n:.1f}%)")


if __name__ == "__main__":
    main()
