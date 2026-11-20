# Track B validation: AdoMetDC - PLP-decarboxylase fusions (17 Nov 2026)

Script: `scripts/p7_fusion_validation.py` (tables `p7_gene_models.tsv`, `p7_same_species.tsv`, `p7_linkers.tsv`,
`p7_partner_integrity.tsv`); structural comparison of the four divergent partners to ODC/ADC/DAPDC models in
`p7_four_structural_vs_refs.txt`. 23 proteins = the 21 in the three multi-protein architectures + 2 singletons.
*Micromonas* C1E5E7 is the one already reported (Green 2011).

## 1. Gene-model validity (rewritten 19 Nov 2026)

**How to read gene structure.** All 20 genome records below are *predicted* gene models. **A predicted intron
between the two domains is not evidence that the fusion is real.** When gene-prediction software merges two
neighbouring genes, it does it by predicting an intron (or several) across the gap between them. That is exactly
what the probable merge in q1 looks like: the HMG-box > AdoMetDC > VTT protein, `results_81/q1_sandwich_report.md`.
The stronger evidence is the opposite case: both domains encoded in one **uninterrupted stretch of genome
sequence**, and above all a **single-exon gene**. Two separate genes could only read through into one exon if they
sat side by side, on the same strand, in frame, with no stop codon between them.

Two caveats still apply to the uninterrupted cases:

- An assembly error, such as a lost stop codon or an indel in an unpolished long-read assembly, can create a
  read-through.
- An uninterrupted ORF does not show that it is transcribed as one mRNA.

RNA-seq or EST reads across the junction are the direct test. They have not been checked for any of the 23.

| Accession | Organism | Evidence | Exons | Across the junction | Same taxon, active UniProt: standalone AdoMetDC / standalone PLP-DC / fusions | Partner integrity (functional residues in place / TM to ODC) |
|---|---|---|---|---|---|---|
| A0A7S1KLJ4 | *Percolomonas cosmopolitus* | **TSA only** (2 TSA records) | - | - | none active (entry deleted) | 2/13, 0.82 |
| A0A9C7PQ51 | *Galdieria partita* | genome | 9 | predicted intron in linker (53 nt) | 0 / 4 / 2 | 1/13, 0.71 |
| A0A1Y1HLG4 | *Klebsormidium nitens* | genome | 6 | predicted intron in linker (355 nt) | 0 / 3 / 1 | 13/13, 0.88 |
| A0AAV8UZF7 | *Rhodosorus marinus* | genome | 5 | predicted intron in linker (73 nt) | 0 / 2 / 1 (+ deleted standalone AdoMetDC A0A7S0BUF8 in this dataset) | 2/13, 0.75 |
| A0AAW1PYA1 | *Myrmecia bisecta* | genome | 8 | predicted intron in linker (342 nt) | 0 / 4 / 1 | 11/13, 0.83 |
| A0A9D4ZGS2 | *Adiantum capillus-veneris* | genome | **1** | **uninterrupted (single exon)** | **1** / 7 / 1 | 13/13, 0.86 |
| A0ABD3ICG1 | *Riccia sorocarpa* | genome | 2 | **uninterrupted (junction inside an exon)** | 0 / 2 / 2 | 13/13, 0.89 |
| A0A2P6V5T4 | *Micractinium conductrix* | genome | 20 | predicted introns in linker (139, 208, 174 nt) | 0 / 3 / 1 | 13/13, 0.87 |
| A0A087SFF4 | *Auxenochlorella protothecoides* | genome | 3 | predicted intron in linker (222 nt) | 0 / 2 / 1 | 12/13, 0.69 (65-residue ODC stretch unaligned) |
| A0A1D2AC22 | *Auxenochlorella protothecoides* | **TSA only** | - | - | same taxon as above | 1/13, 0.15 (partner only 65 aa) |
| A0AAW1RFQ0 | *Apatococcus lobatus* | genome | 8 | **uninterrupted (junction inside an exon)** | 0 / 2 / 1 | 13/13, 0.86 |
| A0AAE0FQH6 | *Cymbomonas tetramitiformis* | genome | 3 | **uninterrupted (junction inside an exon)** | 0 / 4 / 1 | 13/13, 0.87 |
| A0A061RZV3 | *Tetraselmis* sp. GSL018 | **TSA only** (several TSA records) | - | - | none active (entry deleted) | 13/13, 0.86 |
| C1E5E7 | *Micromonas commoda* (reported) | genome | 3 | **uninterrupted (junction inside an exon)** | 0 / 1 / 1 | 13/13, 0.88 |
| A0A9P3J9N9 | *Closterium* sp. NIES-68 | genome (2 records) | 2 | **uninterrupted (junction inside an exon)** | none active (entry deleted) | 13/13, 0.86 |
| A0A5B8MZ45 | *Chloropicon primus* | genome | **1** | **uninterrupted (single exon)** | 0 / 2 / 1 | 13/13, 0.87 |
| A0AAD9ING6 | *Prototheca wickerhamii* | genome | 4 | predicted intron in linker (68 nt) | 0 / 1 / 1 | 5/13, 0.53 (147-residue ODC stretch unaligned) |
| A0AAW1T1C3 | *Apatococcus fuscideae* | genome | 9 | predicted introns in linker (465, 549 nt) | 0 / 2 / 2 | 13/13, 0.84 |
| A0A2V3IZY9 | *Gracilariopsis chorda* | genome | **1** | **uninterrupted (single exon)** | 0 / 3 / 1 | 3/13, 0.69 |
| A0A830HJR7 | *Pycnococcus provasolii* | genome | **1** | **uninterrupted (single exon)** | 0 / 1 / 1 | 13/13, 0.87 |
| A0A8T2T6R3 | *Ceratopteris richardii* | genome | **1** | **uninterrupted (single exon)** | **4** / 9 / **8** | 13/13, 0.87 |
| A0A0G4EYD2 | *Vitrella brassicaformis* | genome | 14 | predicted introns in linker (79, 98 nt) | 0 / 2 (one named ODC, one DAPDC) / 1 | 13/13, 0.87 |
| A0AAW1STH3 | *Apatococcus fuscideae* | genome (partial gene, contig end) | 6 | predicted introns in linker (455, 531 nt) | same taxon as above | 1/13, 0.30 (partner only 143 aa) |

**Counts, corrected.** **3 are transcriptome-only** (*Percolomonas*, *Tetraselmis* GSL018, *Auxenochlorella*
A0A1D2AC22) and **20 are genome gene models** (19 unreported plus *Micromonas* C1E5E7).

- **Stronger group: 10 of the 20 are uninterrupted across the junction.** 5 are single-exon genes (*Adiantum*,
  *Chloropicon*, *Gracilariopsis*, *Pycnococcus*, *Ceratopteris*). In the other 5 the junction falls inside an
  exon of a multi-exon gene (*Riccia*, *Apatococcus lobatus*, *Cymbomonas*, *Closterium*, *Micromonas*). Leaving
  out *Micromonas* (already reported), that is **9 of the 19 unreported**.
- **Weaker group: 10 of the 20 rely on a predicted intron in the linker** (*Galdieria*, *Klebsormidium*,
  *Rhodosorus*, *Myrmecia*, *Micractinium*, *Auxenochlorella* A0A087SFF4, *Prototheca*, *Apatococcus fuscideae*
  ×2, *Vitrella*). All 10 are among the 19 unreported. Their gene structure says nothing for or against the
  fusion until the junction is shown in RNA.

(The earlier version of this section had these two groups the wrong way round: "10 without an intron, 9 with".
The table itself was right. In *A. lobatus* the intron near the junction falls at residue 324, inside the
AdoMetDC domain, which ends at 347, so the linker itself is uninterrupted.)

**Junction-intron size (19 Nov).** A merge across intergenic DNA usually leaves one oversized "intron" at the
junction. In q1 they were 1,125 and 1,456 nt, against a typical 39-62 nt in the same gene. So for the 10 genes
with a linker intron, each junction intron was compared with the other introns of the same gene (from the
ENA CDS joins already cached in `data/ena/`):

| Accession | Organism | Junction intron(s), nt | Other introns in the gene: n, median, max (nt) | Largest in the gene? |
|---|---|---|---|---|
| A0A9C7PQ51 | *Galdieria partita* | 53 | 7, 52, 59 | no (typical) |
| A0A1Y1HLG4 | *Klebsormidium nitens* | 355 | 4, 303, 490 | no |
| A0AAV8UZF7 | *Rhodosorus marinus* | 73 | 3, 79, 90 | no (smallest) |
| A0AAW1PYA1 | *Myrmecia bisecta* | 342 | 6, 390, 1,544 | no |
| A0A2P6V5T4 | *Micractinium conductrix* | 139, 208, 174 | 16, 186.5, 270 | no |
| A0A087SFF4 | *Auxenochlorella protothecoides* | 222 | 1, 275, 275 | no |
| A0AAD9ING6 | *Prototheca wickerhamii* | 68 | 2, 218, 392 | no (smallest) |
| A0AAW1T1C3 | *Apatococcus fuscideae* | 465, 549 | 6, 472.5, 2,589 | no |
| A0A0G4EYD2 | *Vitrella brassicaformis* | 79, 98 | 11, 87, 136 | no |
| A0AAW1STH3 | *Apatococcus fuscideae* (partial gene) | 455, 531 | 3, 561, 1,092 | no |

**None of the 10 has an oversized junction intron.** Several are the smallest intron in their gene. So the
simplest merge signature (the one q1 shows) is absent.

**This is a missing warning sign, not evidence that the fusions are real.** A merge can also turn intergenic or
UTR DNA into predicted exons, with normal-sized predicted introns between them. Most linkers are low-confidence
in the AlphaFold models (pLDDT 37-77), which fits either reading.

One observation for later: the two *A. fuscideae* fusions have the same linker structure, two introns of similar
size (465/549 and 455/531 nt). That fits a duplicated real gene better than two independent merges, but it has
not been tested.

**Still to add from `scripts/q3_fusion_gene_models.py` (not yet run):**

- the same junction introns against a random sample of 400 genes from the same genome annotation (percentile;
  outlier = largest intron in its gene *and* above the genome's 95th percentile);
- for every species, the other ENA assemblies and UniProt proteomes, with each extra proteome scanned for the
  fusion versus separate AdoMetDC and ODC genes.

Results will go in `results_81/q3_*.tsv`.

A standalone AdoMetDC from the same taxon was found only for *Adiantum* (1), *Ceratopteris* (4) and *Rhodosorus*
(1, a deleted entry in this dataset). For every other taxon with active entries, the fusion is the only
AdoMetDC-Pfam protein in UniProt. Standalone PLP-decarboxylase proteins exist in all taxa with active entries.
They were not classified as ODC/ADC/DAPDC here, because Pfam alone does not separate them.

**Limits.** Genome gene models are predictions. No RNA-seq or EST support for any junction has been checked.
Deleted UniProt entries cannot be searched, so "0 standalone" means "none among active entries".

## 2. The four with a replaced ODC lysine

| | *Percolomonas* | *Galdieria* | *Rhodosorus* | *Gracilariopsis* | controls |
|---|---|---|---|---|---|
| partner length, pLDDT | 435 aa, 74.6 | 354 aa, 78.5 | 395 aa, 78.5 | 357 aa, 75.8 | |
| TM to human ODC (by ODC / by partner) | 0.82 / 0.80 | 0.71 / 0.83 | 0.75 / 0.80 | 0.69 / 0.80 | *Klebsormidium* 0.88 / 0.83; *T. brucei* ODC 0.94; AZIN1 0.90 |
| TM to plant ADC / *E. coli* SpeA / LysA (by partner) | 0.76 / 0.76 / 0.74 | 0.78 / 0.78 / 0.80 | 0.78 / 0.76 / 0.77 | 0.77 / 0.77 / 0.80 | *Klebsormidium* 0.76 / 0.75 / 0.76 |
| sequence identity to human ODC | 29% | 20% | 24% | 21% | *Klebsormidium* 44%; AZIN1 54% |
| nearest reference by sequence | AZIN1, then ODC (score 371 vs 354) | AZIN1, then ODC (218 vs 199) | ODC human (243), AZIN1 (229) | ODC human (141), tobacco ODC (137) | none of the four is nearest to ADC, DAPDC or CANSDC |
| ODC functional residues identical and in place (of 13) | 2: C360, D361 | 1: D361 | 2: P275, G276 | 3: G276, C360*, D361* | *T. brucei* ODC 13; AZIN1 9; AZIN2 10 |
| K69 (PLP Schiff base) | S | R | Q (structure) / R (sequence alignment) | R | |
| C360 (catalytic) | C | **G** | C (4.4 A off) | C (3.5 A off) | |
| PLP-binding E274 / R277 | T / S | N / S | D / S | V / A | |
| substrate Y331 / D332 | V / M (displaced 6-9 A) | T / E | D / A (displaced) | M / G | |

Reading: all four partner domains are complete, confidently predicted PLP-decarboxylase-fold domains (fold coverage
and TM comparable to the intact *Klebsormidium* partner), but they keep almost none of the ODC residues that bind PLP,
form the Schiff base or bind substrate - fewer than the catalytically inactive AZIN1 keeps (9/13). Structure does not
discriminate ODC from ADC/DAPDC folds (all 0.74-0.84). By sequence their nearest references are ODC and the
antizyme inhibitor, never ADC. Whether they are inactive, or active on another substrate with another residue set, is
not decided by this comparison.

Also degraded or incomplete among the 23 (not in the four): *Prototheca* A0AAD9ING6 (147-residue ODC stretch missing,
5/13), *Auxenochlorella* A0A087SFF4 (65-residue stretch incl. K69 region missing), and the two truncated partners
(A0A1D2AC22, 65 aa; A0AAW1STH3, 143 aa).

## 3. Organism-specific literature (Europe PMC, 17 Nov 2026)

Per organism: (a) organism name + polyamine / ODC / AdoMetDC / spermidine / putrescine, any field; (b) genome or
transcriptome papers by title, with every open-access full text searched for those terms.

| Organism | Polyamine-related reports | Genome/transcriptome papers checked in full text | Fusion mentioned? |
|---|---|---|---|
| *Klebsormidium* | none specific (55 co-mention hits, none on polyamine enzymes of *Klebsormidium*) | Hori 2014 Nat Commun genome (PMC4052687): no polyamine terms | not found |
| *Galdieria* (*partita*) | none specific | Hirooka 2022 PNAS *G. partita* genome (PMC9565259, linked from the UniProt entry; found 19 Nov, earlier version of this table wrongly said none): no terms; *G. sulphuraria* organelle/splicing genome papers: no terms | not found |
| *Vitrella* | Woo 2015 eLife chromerid genomes (PMC4501334): "pathways shared between chromerids and ... Plasmodium, polyamine synthesis" (text; Fig. 2 supplement 3 not checked) | same | **uncertain**: text says chromerids share polyamine synthesis with *Plasmodium* but does not mention a fusion; the figure supplement may show more |
| *Percolomonas* | Li 2026 JBC: thermospermine synthases from *Percolomonas* TSA, nothing on AdoMetDC/ODC | no genome paper | not found |
| *Rhodosorus* | none | none | not found |
| *Gracilariopsis chorda* | none | Lee 2018 draft genome (PMID 29688518): not open access, not checked | not found (genome paper unchecked) |
| *Myrmecia* / *Apatococcus* | none (2024 lichen-alga phylogenomics hit is a buffer recipe) | none | not found |
| *Micractinium conductrix* | none | 2026 genome (PMC13543684): no terms; 2018 genome (PMID 29178410) not open access | not found (one genome paper unchecked) |
| *Auxenochlorella* / *Prototheca* | none | 2024 draft genome note, 2021/2022/2026 *Prototheca* genomes, plastid genomes: no terms | not found |
| *Cymbomonas* | none | 2023 long-read genome (PMC10675990): no terms | not found |
| *Tetraselmis* | none specific | 4 genome papers (other species): no terms | not found |
| *Closterium* | none | 2023 *C. peracerosum* genome (PMC10407611): no terms | not found |
| *Chloropicon primus* | none | 2019 genome (PMC6731263): no terms | not found |
| *Pycnococcus provasolii* | none | 2023 genome note (PMC11130579): no terms | not found |
| *Ceratopteris* | none | 2019 genome (PMC6890710), 2021 transcriptome: no terms | not found |
| *Adiantum* | none | no genome paper in Europe PMC by title | not found |
| *Riccia sorocarpa* | none | 2025 T2T genome (PMC11735767): no terms | not found |
| *Micromonas* | Green 2011 (reported) | 2009 genomes paper (PMID 19359590) not open access | reported |
| Lin & Lin 2018 microalgae polyamine review (PMC6356823) | mentions *Chlorella* ODC/ADC activities, PBCV-1 ADC, red-alga SAMDC inhibitor work; no AdoMetDC-ODC fusion | - | - |

## Linker / hinge between the two domains

`p7_linkers.tsv`: residues between the last AdoMetDC Pfam match and the first PLP-decarboxylase Pfam match.
Lengths 41-202 aa (median 133); mean pLDDT 37-77, with 0-91% of linker residues below pLDDT 50 (most are largely
low-confidence). *P. falciparum* Q8IJ77, measured the same way: 359 aa, with 21% Asn/Gln (the 23 algal linkers have
1-17% N/Q). Per-protein sequences are in the table.
