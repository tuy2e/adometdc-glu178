# p23: long-read (and other) RNA evidence for the fusion junctions

Queried 23 Nov 2026: ENA portal API `result=read_run` (searches the INSDC/SRA run archive), `tax_tree(<taxid>)`,
platform PACBIO_SMRT or OXFORD_NANOPORE; separately every TRANSCRIPTOMIC run on any platform; and `result=tsa_set`.
No identifier sent. Raw API answers are in `data/p23/`. GenBank records checked with NCBI efetch.

## Why p22 had no such files

The p22 task I received had three parts (tblastn queries, ColabFold heterodimers, census) and no long-read or
unmodelled-sequence part, so these files were never attempted; nothing failed. If the p22 prompt you sent had
more parts, they did not arrive.

## Answer per species

| species (taxid) | long-read RNA runs | long-read DNA runs | other RNA | status |
|---|---|---|---|---|
| *Galdieria partita* (83374) | **none** | 1 PacBio (genome project PRJDB12742) | **59 Illumina RNA-seq runs** (paired, 2 x 150), no TSA assembly | needs a check, see below |
| *Gracilariopsis chorda* (448386) | **none** | none | 23 RNA-seq runs | no long-read RNA |
| *Klebsormidium nitens* (105231) | **none** | none | 173 RNA-seq runs; TSA (see below) | **done** by TSA |
| *Ceratopteris richardii* (49495) | **4 PacBio RS II transcriptome runs** | 183 (PacBio and Nanopore, genome) | 370 RNA-seq runs | long-read RNA exists |

### Ceratopteris richardii: long-read RNA

| run | platform | reads (as reported) | bases | BioProject | title |
|---|---|---|---|---|---|
| SRR12605702 | PacBio RS II | 287,447,857 | 255.7 Gb | PRJNA651774 | *C. richardii* M. Whitten #5841 Transcriptome |
| SRR12605704 | PacBio RS II | 123,257,023 | 317.7 Gb | PRJNA651774 | same |
| SRR12605705 | PacBio RS II | 208,691,015 | 195.2 Gb | PRJNA651775 | *C. richardii* Hn-n (M. Whitten #5841) Transcriptome |
| SRR12605706 | PacBio RS II | 128,896,418 | 298.2 Gb | PRJNA651775 | same |

The read counts are what SRA reports for these RS II runs and are probably subreads, not full-length transcripts;
library strategy is "OTHER" (not labelled Iso-Seq). Neither BioProject is the genome BioProject in q3
(PRJNA511033, PRJNA729743), so a hit would be independent of the genome. The *Ceratopteris* fusion is a single-exon
gene (p21f), so the open question there is read-through or assembly error, not a merge across an intron.

### Galdieria partita: no long-read RNA, no TSA, but short-read RNA exists

Your note is right that *G. partita* has nothing in TSA (ENA lists no TSA set for any *Galdieria* species) and that
the two full-coverage nt hits are RefSeq models: XM_005709307.1 and XM_005709306.1 are *G. sulphuraria* 074W
records "derived from an annotated genomic sequence (NW_005178475)", not observed transcripts.

But "no RNA" is not right: SRA holds 59 *G. partita* RNA-seq runs (all Illumina, paired 2 x 150, mean read length
about 300 nt per spot):

| BioProject | runs | instrument | strain | example run | median reads/run |
|---|---|---|---|---|---|
| PRJDB19545 | 42 | NovaSeq 6000 | NBRC 102759 | DRR620710 | 12.3 M |
| PRJDB12743 | 12 | NovaSeq 6000 | NBRC 102759 | DRR337147 | 16.8 M |
| PRJNA507473 | 5 | HiSeq 2500 | not given | SRR8259111 | 13.2 M |

NBRC 102759 is the genome strain (genome BioProject PRJDB12742). PRJDB12743 is the next number and the same strain,
so it is probably the genome project's companion RNA-seq; evidence from it would not be fully independent.
PRJDB19545 and PRJNA507473 are separate projects.

The premise that short reads cannot help is only partly true here. No single read spans the whole 234-aa query,
but the fusion's gene model puts one **53-nt intron in the linker** (p7). Short reads that splice across that
intron, together with read pairs whose mates fall in the AdoMetDC and ODC exons, test the junction directly. How to
search them (not run):

- NCBI web BLAST, tblastn, database **"SRA"** (Sequence Read Archive), entering one run accession (e.g. DRR620710),
  query `data/p22_tblastn_queries/01_Galdieria_partita_A0A9C7PQ51.txt`. Hits are individual reads; look for reads
  covering the linker.
- Better: download one run and align it to the *G. partita* genome contig BQMJ01000002.1 with a splice-aware aligner
  (e.g. HISAT2 or STAR); then look for reads split at the predicted linker intron and pairs spanning it.

### Galdieria yellowstonensis: the one long-read RNA run in the genus

| run | platform | reads | bases | BioProject | title |
|---|---|---|---|---|---|
| SRR19760108 | PacBio Sequel, cDNA (Iso-Seq) | 21,609 | 37.3 Mb | PRJNA851236 | "Galdieria sulphuraria 108.79 E11 PacBio Iso-Seq" |

Strain 108.79 E11 is the strain of the *G. yellowstonensis* fusion A0AAV9I966 (its ENA CDS, p18), now classified as
*G. yellowstonensis*; the run title uses the old name. Its ODC half is also degraded (p10). This is not *G. partita*,
but it is the only long-read RNA in *Galdieria*, from the same strain as a degraded fusion. It is small (21,609
reads), so absence of a hit would mean little. How to search it (not run): tblastn against SRA run SRR19760108 on
the NCBI web form, or download the 37 Mb of reads and search locally; use a junction query built from A0AAV9I966
by the p14 recipe (the *G. partita* query should also hit at about 72 % identity, p18).

### Klebsormidium nitens: done, with a naming caveat

The three TSA records from your BLAST were checked:

| record | organism as filed | BioProject | type |
|---|---|---|---|
| JO257453.1 | *Klebsormidium flaccidum* strain UTEX 321 | PRJNA51159 | TSA mRNA, 2,307 bp |
| GBSO01022698.1 | *Klebsormidium flaccidum* | PRJNA242253 | TSA, 3,595 bp |
| GFXY01006605.1 | *Klebsormidium flaccidum* | PRJNA399177 | TSA, 3,033 bp |

Three separate projects, none the genome project (PRJDB718). All are filed as *K. flaccidum*, not *K. nitens*
NIES-2285 (the genome strain, formerly also called *K. flaccidum*); JO257453 is strain UTEX 321. Matches at 100 %
mean the same sequence, so for practical purposes the transcript is the fusion's, but write it as "transcripts from
*Klebsormidium* (filed as *K. flaccidum*)" rather than "from *K. nitens*". Identity and coverage are from your
manual BLAST; I did not rerun it.

### Gracilariopsis chorda

No long-read data of any kind. 23 short-read RNA-seq runs exist (`data/p23/rna_any_448386.tsv`); its fusion is
single-exon, so a splice test does not apply there.
