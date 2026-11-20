# AdoMetDC survey of the 1,000-hit run — package for Marc

Tui, 19 November 2026. Everything here comes from the `stuy__adometdc` FOLDTRACE run.
`aim1_settings_used.json` and the `models/` folder were never written to. Nothing in this
package has been sent to anyone else.

## What this is

A broad survey of all 1,000 hits, rather than only scoring the five prespecified sites.
**773 pass the fold gate** (TM ≥ 0.5 against 1MSV chain A); everything below uses those.

Before anything else: **326 of the UniProt accessions had been deleted** since the models
were built — 319 of them for "not part of a reference proteome", which is housekeeping
rather than a judgement on the protein. 211 of those are fold-passers. All 211 were
recovered through UniParc, where the sequences match the models exactly, so nothing was
dropped. The deleted set skews toward under-studied organisms, which matters for any
argument of the form "there is no other copy in UniProt".

## Architecture census

Of the fold-passers: **590 standalone · 81 with two PF01536 matches in one chain ·
54 partner-then-AdoMetDC · 46 AdoMetDC-then-partner · 1 sandwiched**
(`A0A8S1GM37`, HMG-box > AdoMetDC > VTT, *Caenorhabditis auricularia*).

61 distinct Pfam architectures in all; 43 held by a single protein. Every architecture
with three or more proteins also appears in the InterPro-only subset, so the rescued
annotation is not generating the pattern.

## One negative result, stated first

I tested whether fused proteins lose the prespecified catalytic and processing sites more
often than standalone ones.

**They do not.** Standalone 25/590 (4.2%) vs fused 3/101 (3.0%), Fisher p = 0.785;
within-phylum permutation p = 0.717; MH odds ratio 0.96. InterPro-only subset: 3.6% vs
3.6%, p = 1.

An earlier version of this using a length-based proxy for "fused" gave p = 0.004. That was
wrong — the proxy counted long N-terminal extensions with no partner domain as fusions. I
am recording it because it nearly became a claim. Fused models also have lower pLDDT than
standalone ones (75.1 vs 83.5), so the original comparison was confounded in the direction
that would manufacture the result.

Site losses concentrate in Euglenozoa (7 of 11), which is the kinetoplastid story we
already had.

→ `p5_loss_vs_fusion_report.txt`

## The finding

**23 proteins are AdoMetDC fused to an ODC-family PLP decarboxylase.**

That fusion is described in *Plasmodium* (Müller 2000, Krause 2000, Wrenger 2001,
Birkholtz 2004, Hart 2016) and reported once in *Micromonas* (Green 2011, sequence only —
and that exact protein, `C1E5E7`, is in this set).

**The other 22 did not turn up in the literature.** They span green algae (*Myrmecia*,
*Micractinium*, *Apatococcus*, *Cymbomonas*, *Tetraselmis*, *Chloropicon*, *Prototheca*,
*Pycnococcus*, *Auxenochlorella*), streptophytes including two ferns and a liverwort
(*Klebsormidium*, *Closterium*, *Adiantum*, *Ceratopteris*, *Riccia*), red algae
(*Galdieria*, *Rhodosorus*, *Gracilariopsis*), plus *Percolomonas* (Heterolobosea) and
*Vitrella* (chromerid).

16 of the 23 partners are confidently ODC by alignment: 44–54% identity, with Lys69 and
Cys360 both present.

I checked Michael's own papers specifically, since his lab is the obvious pre-emption risk
— Green 2011, Li 2023 *JBC*, Li 2026 *JBC*, Li 2023 *PNAS*, Lin & Lin 2018, Fuell 2010 —
and the genome paper for nearly every organism on the list.

→ `trackB_fusion_partners.md`, `p6_decarboxylase_partner.tsv`

## Are they real genes?

**19 of the 22 are genome gene models.** Three are transcriptome-only (*Percolomonas*,
*Tetraselmis* GSL018, *Auxenochlorella* `A0A1D2AC22`) and I have set those aside.

Three things argue against assembly artifacts:

- **In every taxon except *Adiantum*, *Ceratopteris* and *Rhodosorus*, the fusion is the
  only AdoMetDC-Pfam protein among active UniProt entries** — the same situation as
  *Plasmodium*. Where a merged gene model is to blame, the separate component genes are
  usually still annotated somewhere in the genome.
- **Multiple copies in one genome.** *Ceratopteris richardii* has eight; *Galdieria*,
  *Riccia* and *Apatococcus fuscideae* have two each. Eight independent mis-merges in one
  genome is not credible.
- **Ten of the 19 have no intron at the domain junction at all**, five of those being
  single-exon genes (*Adiantum*, *Chloropicon*, *Gracilariopsis*, *Pycnococcus*,
  *Ceratopteris*). Two domains encoded contiguously cannot be produced by merging two
  genes separated by intergenic DNA.

**Unresolved: the 9 with a predicted intron in the linker.** I had this backwards in an
earlier draft. Gene-prediction software merges neighbouring genes by calling the
intergenic DNA an intron, so a *predicted* intron at the junction is the signature of the
artefact rather than evidence against it — this is exactly how the sandwich protein
(`A0A8S1GM37`, a separate lead) was exposed, its extra domains hanging off across the two
largest introns in the model. Those 9 need their junction intron compared against the
intron-length distribution of the same gene and genome.

Two further caveats: "no standalone copy" can only mean "none among active entries", per
the deletion caveat above; and no RNA-seq or EST support across any junction has been
checked. Nor has the obvious test of whether any of these species has a second genome
annotation lacking the fusion — that test is what exposed a nematode AdoMetDC–VTT fusion
in this dataset as an artefact.

→ `trackB_fusion_validation.md`, `p7_*.tsv`

## A second observation, possibly separate

**Four of the 23 have a structurally complete but catalytically gutted partner domain.**

| | ODC functional residues kept (of 13) | Schiff-base Lys69 | Catalytic Cys360 |
|---|---|---|---|
| *Percolomonas* | 2 | Ser | present |
| *Galdieria* | 1 | Arg | **Gly** |
| *Rhodosorus* | 2 | Arg | 4.4 Å displaced |
| *Gracilariopsis* | 3 | Arg | 3.5 Å displaced |
| AZIN1 (known dead ODC paralogue) | 9 | | |
| *T. brucei* ODC | 13 | | |

Fold TM-scores and partner pLDDT (75–79) are comparable to the intact partners, so this is
not a modelling failure. A dead enzyme domain fused to a live one.

Caveat: two of the four score nearest to AZIN1 by sequence, but that most likely just means
"divergent" — the antizyme system is not known outside animals, so I am not reading them as
antizyme inhibitors.

## Does the fused AdoMetDC half keep the putrescine site? Mostly a negative

Eight positions in 1MSV numbering mapped onto 23 ODC fusions, 78 other fused proteins and
590 standalone, same TM-align method as the residue grid, pLDDT ≥ 70 or the call is
unresolved.

- **Glu178 is kept in all 23** — but standalone keeps it 514/579 (89%) anyway, so there is
  no room for a difference. Fisher p = 0.159, permutation p = 0.517.
- **Leu13→Arg (20 of 23) and Phe111→Lys (15 of 23) are not a fusion signal.** Standalone
  proteins do the same: 133 Arg and 58 Lys at position 13 among 590. Fisher p = 5.5e-07
  collapses to permutation p = 0.667. Same for Phe111, Phe285, Asp174, Thr176. The basic
  filled pocket is a plant and algal trait, not a consequence of fusion.
- **One survivor: Glu15.** Standalone 424/459 kept (92%); the 23 fusions 13/22 (59%), five
  of the nine losses to Ser. Fisher p = 4.05e-05, **within-phylum permutation p = 0.0055**.
  The other 78 fused proteins show nothing (p = 0.647), so it is specific to this fusion
  rather than to fusion generally. It holds in the InterPro-only subset (p = 0.043).

**I am treating Glu15 as a lead, not a result.** Only five phyla contain both groups, and
eight positions were tested, so it barely clears correction for eight and not for
everything run that day.

→ `p8_putrescine_site_report.txt`, `p8_putrescine_site_tests.tsv`

## What is still open

1. **The thorough literature search you asked for** has not run yet. The search so far was
   organism by organism; the next one comes at the fusion as a concept, your full
   publication list for Michael and the structural people, the 23 accessions as literal
   strings, and preprints.
2. **Woo et al. 2015 (*eLife*, PMC4501334), Figure 2 supplement 3, is unchecked.** The text
   says chromerids share polyamine synthesis with *Plasmodium*. *Vitrella* is the organism
   most at risk of already being described.
3. Three genome papers are paywalled: Lee 2018 (*Gracilariopsis*), the 2018 *Micractinium*
   genome, the 2009 *Micromonas* genomes paper.
4. **Whether the AdoMetDC half still self-processes.** Ser68, Ser229, His243, Cys82,
   Glu247 across all 23 — not yet checked, and it changes the reading if any are dead.
5. **One fusion event or many.** Planned both ways: phylogeny of each half separately, and
   structurally by superposing all 23 on their AdoMetDC domain to see whether the partner
   lands in the same orientation.
6. **Base rate.** The same census on ODC, spermidine synthase and SAM synthetase, so
   "13% of fold-passers are fused" has something to compare against.
7. **A new FOLDTRACE run using a fusion as the query**, as you suggested. This would be its
   own run with its own settings file, not filed under the Aim 1 protocol.

## Questions for you

1. Is the two-tree plus structural-orientation approach the right way to settle one fusion
   versus many, or is there something better?
2. Does this stand on its own, or is it a case study for the FOLDTRACE paper?
3. On pre-emption — would you ask Michael directly, or would you rather I did not flag it
   yet?
4. Is the dead-partner observation worth separating from the main fusion story?

The Glu178 regulatory work is unaffected. One thing from this run that bears on it: of the
65 standalone proteins that lost Glu178, **30 went to Asn** — the diatom pattern, sitting
in this survey's control group.

## Files

| File | What it is |
|---|---|
| `RESEARCH_LOG.md` | Full running record; this survey is Thread 6 |
| `survey_README.md` | What every script does |
| `architecture_census_summary.txt` | The 61 architectures |
| `fusion_architectures.tsv` | Per-protein fusion layouts |
| `p5_loss_vs_fusion_report.txt` | Track A, the refuted one |
| `trackB_fusion_partners.md` | The 23 fusions and the literature check |
| `p6_decarboxylase_partner.tsv` | Partner identity vs ODC/ADC/DAPDC references |
| `trackB_fusion_validation.md` | Gene models, the four gutted partners, per-organism literature |
| `p7_gene_models.tsv` | Genome vs transcriptome, exons, introns in the linker |
| `p7_same_species.tsv` | Standalone copies in the same taxon |
| `p7_partner_integrity.tsv` | The 13 ODC functional residues, per protein |
| `p7_linkers.tsv` | Linker sequence, length and pLDDT |
| `p7_four_structural_vs_refs.txt` | The four divergent partners vs ODC/ADC/DAPDC |
| `p8_putrescine_site_report.txt` | Putrescine-site counts and tests |
| `p8_putrescine_site_tests.tsv` | The test table |
| `p8_putrescine_site.tsv` | Per-protein, per-position calls for all 691 |

Code and data are in the repo; this package is the readable subset.
