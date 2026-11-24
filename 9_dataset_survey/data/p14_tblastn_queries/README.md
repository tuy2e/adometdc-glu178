# tblastn queries: does the fusion appear on a single transcript?

Built 22 Nov 2026 from the AlphaFold models in Marc's read-only `stuy__adometdc/models/`
(sequence read from CA atoms) and the domain boundaries in `results/p7_linkers.tsv`.
Marc's folder was only read, never written to.

## What these are

Each query spans the domain junction: the **last 60 residues of the AdoMetDC Pfam
match**, the **entire linker**, and the **first 60 residues of the PLP-decarboxylase
match**. Querying the linker alone would prove nothing - a hit there does not show the
two enzymes are on one molecule. A query that covers AdoMetDC *and* linker *and* ODC
only aligns end-to-end against a transcript that genuinely carries both.

## How to run it

1. blast.ncbi.nlm.nih.gov -> **tblastn**
2. Paste the contents of one `.txt` file into the query box (no header needed)
3. Database: **Transcriptome Shotgun Assembly (TSA)**
4. Organism box: type the organism name below; NCBI autocompletes it. If nothing comes
   back, retry with the genus alone, then with the database set to **Nucleotide
   collection (nr/nt)**.
5. BLAST

## What counts as a positive

**One TSA contig aligning across the whole query, contiguously, with no gap at the
junction.** That means a single RNA molecule carried both halves, and the gene is real.

Partial hits covering only the AdoMetDC side or only the ODC side prove nothing - those
are just the two separate enzymes, which every one of these species has anyway.

## Caveats to record with any result

- If a species' only TSA assembly comes from the same sequencing project as its genome,
  a hit is **not independent evidence** - the same reads built both. Check the TSA
  BioProject against the genome BioProject in `results/q3_assemblies.tsv` before
  claiming independence.
- Absence of a hit means nothing on its own. Many of these species have little or no
  public RNA data.

## The queries

| # | file | organism | why it is a priority | query length | AdoMetDC / linker / ODC |
|---|---|---|---|---|---|
| 1 | `01_Adiantum_capillus-veneris_A0A9D4ZGS2.txt` | *Adiantum capillus-veneris* | single-exon gene - strongest candidate, cannot be a gene merge | 219 aa | 60 / 99 / 60 |
| 2 | `02_Chloropicon_primus_A0A5B8MZ45.txt` | *Chloropicon primus* | single-exon gene - strongest candidate, cannot be a gene merge | 264 aa | 60 / 144 / 60 |
| 3 | `03_Gracilariopsis_chorda_A0A2V3IZY9.txt` | *Gracilariopsis chorda* | single-exon gene - strongest candidate, cannot be a gene merge | 200 aa | 60 / 80 / 60 |
| 4 | `04_Pycnococcus_provasolii_A0A830HJR7.txt` | *Pycnococcus provasolii* | single-exon gene - strongest candidate, cannot be a gene merge | 309 aa | 60 / 189 / 60 |
| 5 | `05_Ceratopteris_richardii_A0A8T2T6R3.txt` | *Ceratopteris richardii* | single-exon gene - strongest candidate, cannot be a gene merge | 230 aa | 60 / 110 / 60 |
| 6 | `06_Percolomonas_cosmopolitus_A0A7S1KLJ4.txt` | *Percolomonas cosmopolitus* | transcriptome-only - a transcript is its ONLY possible evidence | 191 aa | 60 / 71 / 60 |
| 7 | `07_Tetraselmis_sp_GSL018_A0A061RZV3.txt` | *Tetraselmis sp GSL018* | transcriptome-only - a transcript is its ONLY possible evidence | 278 aa | 60 / 158 / 60 |
| 8 | `08_Auxenochlorella_protothecoides_A0A1D2AC22.txt` | *Auxenochlorella protothecoides* | transcriptome-only - a transcript is its ONLY possible evidence | 265 aa | 60 / 165 / 40 |

## One query is weaker than the others

`08_Auxenochlorella_protothecoides_A0A1D2AC22` runs to the very end of its model, so it
carries only ~40 residues of the PLP domain rather than 60. That protein's partner
region is just 65 aa in total - it is the truncated one flagged in p9. A hit still tells
you the transcript exists, but it is thin evidence that both *domains* are on it.

## Suggested order

Run 1-5 first (the single-exon genes). A hit on any of those settles that protein
outright. Then 6-8, where a transcript is the only evidence available at all.
