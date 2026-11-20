# Job 2: the sandwich protein A0A8S1GM37 (19 Nov 2026)

*Caenorhabditis auriculariae* (UniProt spelling; the handoff wrote *auricularia*). 865 aa, Pfam layout
`HMG-box (43-96) > AdoMetDC (259-582) > VTT (691-809)`. The only protein among the 773 fold-passers with
AdoMetDC in the middle.

## Verdict: probably a merged gene model, not a real three-domain protein. Closed.

The AdoMetDC part looks like a complete, normal *smd-1* gene with its own start and end. The HMG-box and VTT
parts are attached outside those ends. That is what you see when gene-prediction software joins neighbouring
genes into one model. It is not proven: RNA-seq reads across the two junctions would settle it (see the end).
Either way, one protein from one gene model in one genome cannot carry a project, so this lead is closed.

## 1. Is the gene model real?

**Source.** Genome, not transcriptome. ENA CDS `CAD6184317.1`, locus `CAUJ_LOCUS236`, on scaffold 1 of
assembly GCA_904845305 (strain NKZ352, submitted by T. Kikuchi, Univ. Miyazaki). It is the only assembly for
this species in ENA. The genome paper (Dayi et al. 2021, *Sci Rep*, PMC7991662, open access) predicted genes
with BRAKER2, using intron hints from one mixed-stage RNA-seq library plus protein hints from 9 nematodes.
Assembly: 109.5 Mb, 491 scaffolds, BUSCO 89.8%. The record says `source:AUGUSTUS` but does not say which
introns had RNA-seq support.

**Exon structure** (`q1_sandwich_exons.tsv`). 23 exons, plus strand, about 8.6 kb of genome. Introns fall
inside both linkers: 6 in HMG>AdoMetDC, 3 in AdoMetDC>VTT. **In an ab initio model this is not evidence that
the protein is real.** When software joins two neighbouring genes, it does so by predicting an intron between
them. (The p7 Track B check counted linker introns as support. That argument is weak for predicted models,
and the same caveat applies there.)

**What points to a merge:**

- **A normal-length AdoMetDC gene sits inside the model.** *C. elegans* SMD-1 (O02655, 368 aa, reviewed)
  aligns end to end onto residues 255-601: 262/347 identical (75.5%), and it covers the last residue of
  SMD-1. Upstream there is an in-frame **Met236**, at the very start of exon 9, 19 residues before the
  alignment starts (SMD-1 has 18 residues before its aligned part). **Met236 to Glu601 is 366 aa.** The nine
  other *Caenorhabditis* AdoMetDCs in UniProt are 366-373 aa.
- **The two biggest introns in the model sit at the two junctions:** 1,125 nt after residue 112
  (HMG>AdoMetDC linker) and 1,456 nt after residue 599, two residues before the end of the SMD-1 match. Most
  other introns are 39-62 nt. The end of the AdoMetDC gene lines up with an exon boundary.
- **The linker before the AdoMetDC is disordered:** residues 97-235 have mean pLDDT 33, and none reaches 70.
  The 42-residue N-terminus is the same (31). The HMG-box (80), AdoMetDC (Pfam part 92) and VTT (85) are each
  confidently folded on their own.

**Other AdoMetDC in *C. auriculariae*.** None. UniProt has one AdoMetDC-Pfam entry for taxon 2777116, and it
is this protein (an active entry in a reference proteome). So the species' only AdoMetDC is inside this
locus. If the model is a merge, the real *smd-1* is the Met236-Glu601 stretch.

**Other *Caenorhabditis* species** (UniProt, genus 6237, AdoMetDC Pfam, active entries):

| Species | Entry | Length | Layout |
|---|---|---|---|
| *C. elegans* | O02655 (smd-1, reviewed) | 368 | AdoMetDC only |
| *C. briggsae* | A8XU83, A0AAE9E3Z7 | 371, 371 | AdoMetDC only |
| *C. remanei* | E3LLT5 | 371 | AdoMetDC only |
| *C. remanei* | A0A261AI60 | 510 | AdoMetDC > VTT |
| *C. latens* | A0A261BMB9 | 371 | AdoMetDC only |
| *C. nigoni*, *C. japonica*, *C. tropicalis* | A0A2G5VP03, A0A8R1DRJ5, A0A1I7T6W7 | 366 each | AdoMetDC only |
| *C. brenneri*, *C. angaria* | G0MJA0, A0A9P1MUD8 | 367, 373 | AdoMetDC only |
| *C. bovis* | A0A8S1EUL9, A0A8S1F128 | 270, 98 | AdoMetDC only (gene apparently split in two) |
| *C. auriculariae* | A0A8S1GM37 | 865 | HMG-box > AdoMetDC > VTT |

- **No other *Caenorhabditis* has the HMG-box.**
- **AdoMetDC > VTT appears in one other *Caenorhabditis*, and only in one of its two annotations.**
  *C. remanei* has two entries from two genome projects (`CRE_` and `FL82_` gene names). One is a standalone
  371-aa AdoMetDC; the other is AdoMetDC > VTT. When two annotations of one species disagree like this, it
  points to annotation rather than biology.
- Outside the genus, two more nematodes in our 1,000 hits have AdoMetDC > VTT: *Parelaphostrongylus tenuis*
  (A0AAD5WKP5) and *Trichostrongylus colubriformis* (A0AAN8F2B4). This is the "VTT (nematodes)" partner in
  the Track B list. **Because it recurs, the AdoMetDC/VTT pairing may be worth one look in lead 2**: either
  real fusions or a gene that often sits next to *smd-1*. I could not check gene order. None of the three
  *C. elegans* VTT genes (`egli-1` D2013.10, `bus-19` T07F10.4, Y71A12C.2) is on the *smd-1* cosmid F47G4,
  and WormBase blocked automated access, so the neighbourhood question is still open.
- *C. auriculariae* is the basal-most *Caenorhabditis* (that is the genome paper's title claim), so it has no
  close relative to compare with. Missing HMG-boxes in the other species are therefore weak evidence on their
  own.

## 2. Does the AdoMetDC part keep the five sites?

Yes, all five, all high confidence (FOLDTRACE calls, `master_candidates.tsv`; whole-chain TM-score 0.957 vs
1MSV chain A):

| Site | Residue here | CA offset (Å) | pLDDT | Call |
|---|---|---|---|---|
| Ser68 (pyruvoyl) | S319 | 1.45 | 91.7 | retained |
| Ser229 (processing) | S482 | 0.34 | 98.0 | retained |
| His243 (processing) | H497 | 0.30 | 98.3 | retained |
| Cys82 (catalytic) | C333 | 0.53 | 95.7 | retained |
| Glu247 (substrate) | E501 | 0.45 | 94.9 | retained |

PROSITE's AdoMetDC signature (PS01336) also matches at 313-323. Nothing about the AdoMetDC part is unusual.
It is a normal, apparently active nematode AdoMetDC.

## 3. Literature: HMG-box or VTT fused to AdoMetDC?

Nothing found. Europe PMC searches (19 Nov):

| Query | Hits | Relevant |
|---|---|---|
| (adenosylmethionine decarboxylase OR AdoMetDC OR SAMDC OR AMD1) AND (HMG box OR high mobility group) AND (fusion OR fused OR chimeric OR domain) | 91 | 0. All are full-text co-mentions in expression or review papers |
| (adenosylmethionine decarboxylase OR AdoMetDC OR SAMDC) AND (VTT domain OR DedA OR SNARE-associated OR TMEM41 OR TVP38) | 9 | 0. *E. coli* gene lists and antibiotic-resistance screens |
| "adenosylmethionine decarboxylase" AND "HMG-box" AND "VTT" | 0 | — |
| "Caenorhabditis auriculariae" | 4 | Genome paper only (Dayi 2021); does not discuss AdoMetDC |

General web search was not available in this session (blocked for the organisation). So this covers Europe
PMC (abstracts plus open-access full text), not Google Scholar or preprint servers beyond Europe PMC's
coverage.

## What would settle it (optional, not needed to close)

The genome project's RNA-seq run is public: **DRR252169** (Illumina MiSeq, BioProject PRJDB10634). Spliced
reads joining exon 8 to exon 9 (genome 1354885 / 1354931) or exon 16 to exon 17 (1356391 / 1357848) would mean
the junctions are real transcripts. Reads that start at exon 9 or stop after exon 16 would mean separate
genes. This needs a read aligner (HISAT2) and one scaffold, which is more than an hour's work. I would only do
it if the VTT pairing becomes interesting in lead 2.

## What failed and what I did instead

- The cloud workspace cannot reach UniProt, ENA, NCBI, WormBase or RCSB directly (proxy refusal). I used the
  web-fetch tool for UniProt, and the in-app browser on your computer for ENA and Europe PMC.
- ENA rate-limited the web-fetch tool, so I read the ENA record in the browser instead.
- NCBI efetch refused the fetch tool (robots rule). Not needed once ENA worked.
- WormBase put up a bot check (Cloudflare). I did not try to get past it, so the *C. elegans* gene
  neighbourhood of *smd-1* is unchecked.
- I could not fetch the raw scaffold sequence to look for a stop codon right after exon 16. The browser
  download was blocked, and the scaffold is several megabases.
- General web search is turned off for this organisation, so the literature check is Europe PMC only.

Files: `q1_sandwich_exons.tsv` (exon table with each intron placed in protein coordinates),
`../scripts/q1_sandwich_checks.py` (reproduces the exon table and the SMD-1 alignment offline; inputs pasted
in from the records, so it needs no network).
