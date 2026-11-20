# Track B step 3: the red-algal pattern in the ODC half

Script: `scripts/p10_red_algal_pattern.py`.

**Status (20 Nov 2026): preliminary.** This write-up comes from the check run (`--smoke`), which covered:

- the full search;
- the 7 red-algal fusions that were new;
- the four species' enzyme inventories;
- parts 4 and 5.

Two things wait for the full run:

- the structure checks on the other new fusions (dinoflagellates, a glaucophyte, a brown alga, CRuMs, *Mesodinium*, *Plasmodium*);
- the tree, which needs `pyfamsa`.

Generated tables:

| file | contents |
|---|---|
| `p10_fusion_search.tsv` | every sequence with both domains |
| `p10_all_fusions.tsv` | per-fusion status |
| `p10_degraded_substitutions.tsv` | the 13 ODC sites in each degraded protein |
| `p10_degraded_pairs.tsv` | pairwise comparison of the substitutions |
| `p10_species_enzymes.tsv` | enzyme inventory of the four species |
| `p10_correlates_23.tsv` | part 5 variables for the 23 |
| `p10_report_smoke.md` | the check run's generated report; the full run writes `p10_report.md` |

No p-values are computed anywhere.

## The pattern, stated precisely

As first stated, the pattern ("every green-lineage fusion has an intact ODC half") does not hold.

- **Green fusions whose ODC half is not intact:**
  - Three green-lineage fusions have an incomplete ODC half: *Prototheca*, *Auxenochlorella* AC22 and *Apatococcus* STH3. They are truncated or have a 147-residue gap.
  - Four have an undetermined half: *Myrmecia*, *Auxenochlorella* SFF4, *Closterium* and *Vitrella*. *Vitrella* is not green.
- **What does hold in every protein with a model so far:** every confident in-place substitution at the 13 human-ODC functional residues is in a red alga or *Percolomonas*. A confident substitution is a residue that sits in place at pLDDT ≥ 70 but is not the human residue.
  - In our 23 that means 4 of 4 proteins. With the new red algae it is 8 of 8: 7 Rhodophyta plus *Percolomonas*.
  - None of the 18 green-lineage fusions has one.

"Degraded" below means exactly this: at least one confident substitution. "Incomplete" and "undetermined" are separate categories (definitions in the script docstring).

## 1. Grow the sample

The Foldseek cap did hide fusions, and many of them.

**Searches:**
- **UniParc**, Pfam features with no structure filter: `(feature_id:PF01536 OR feature_id:PF02675) AND (feature_id:PF02784 OR feature_id:PF00278)`. It returns 311 sequences: 23 are ours and 288 are new.
- **UniProtKB**, the same query on Pfam cross-references, as a check: 86 active entries, all of them among the UniParc hits.

| group | new sequences | AdoMetDC before the PLP domain |
|---|---|---|
| Alveolata | 114 | 114: about 70 *Plasmodium* strain/species sequences; about 40 dinoflagellates (*Alexandrium*, *Scrippsiella*, *Heterocapsa*, *Karlodinium*, *Lingulaulax* and others); *Chromera velia*; *Cardiosporidium*; *Haemoproteus*; *Mesodinium pulex* (ciliate) |
| Viridiplantae | 92 | 92 |
| Bacteria | 66 | 2 (64 are the actinobacterial ODC→AdoMetDC order, a different architecture) |
| **Rhodophyta** | **10** | 10 |
| Glaucophyta | 1 | 1: *Gloeochaete wittrockiana* |
| Stramenopiles | 1 | 1: *Heribaudiella fluviatilis* (brown alga) |
| other | 3 | 3: *Collodictyon triciliatum* ×2, *Diphylleia rotans* (CRuMs) |
| Metazoa | 1 | 1: *Piliocolobus* (monkey genome; most likely a *Hepatocystis* parasite contaminant, unchecked) |

**The 10 new red-algal sequences:**

| class | sequences |
|---|---|
| Cyanidiophyceae | *Galdieria partita* (a second fusion, A0A9C7PSW2); *G. sulphuraria* ×2; *G. yellowstonensis* |
| Bangiophyceae | *Porphyra umbilicalis*; *Pyropia yezoensis*; *Pyropia vietnamensis* |
| Florideophyceae | *Gracilaria gracilis*; *Agarophyton chilense*; *Delisea pulchra* |

- The red algae with a fusion now span three classes, four with *Rhodosorus* (Stylonematophyceae).
- No Cryptophyta was found.
- Heterolobosea is still only *Percolomonas*.

**The 7 new red-algal fusions with an AlphaFold model, analysed with the p9 gates:**

| accession | organism | AdoMetDC half | ODC half | ODC kept / lost / unresolved |
|---|---|---|---|---|
| A0A9C7PSW2 | *Galdieria partita* (2nd copy) | intact | degraded | 1/7/5 |
| M2YAF1 | *Galdieria sulphuraria* | intact | degraded | 1/6/6 |
| A0AAV9I966 | *Galdieria yellowstonensis* | intact | degraded | 1/6/6 |
| A0A1X6PGY4 | *Porphyra umbilicalis* | undetermined | degraded | 0/6/7 |

- Every red-algal fusion with a model is degraded: 7 of 7, counting our three.
- The other six red-algal sequences have no AlphaFold model; five of the six have no UniProt entry at all.
  - By nearest reference (p6 method) their ODC halves are closest to ODC, AZIN1 or *T. brucei* ODC, and in one case (*Gracilaria*) to the plant LDC.
  - This tells us nothing about their sites. The full run writes their sequence-read sites into `p10_all_fusions.tsv`.

**Too small to conclude:** red-algal fusions with models are n = 7, from 4 genera; Heterolobosea is n = 1.

## 2. One event or many

### Substitutions, residue by residue

Notation: upper case = confident substitution; lower case = unresolved; `.` = human residue kept; `-` = nothing aligned.

| site | *Percolomonas* | *G. partita* PQ51 | *G. partita* PSW2 | *G. sulphuraria* | *G. yellowstonensis* | *Rhodosorus* | *Gracilariopsis* | *Porphyra* |
|---|---|---|---|---|---|---|---|---|
| K69 | S | R | R | R | R | q | R | R |
| H197 | V | S | S | s | c | Q | E | E |
| S200 | f | - | - | - | - | e | - | a |
| G237 | s | s | s | s | n | s | k | g |
| E274 | T | N | N | N | N | D | V | G |
| P275 | S | V | V | A | A | . | V | A |
| G276 | A | D | D | D | D | . | . | S |
| R277 | S | S | S | S | S | S | A | P |
| Y331 | v | t | t | t | t | d | m | l |
| D332 | m | e | e | e | e | a | g | g |
| C360 | . | G | G | G | g | c | . | c |
| D361 | . | . | . | . | . | E | . | d |
| Y389 | v | f | f | f | F | T | q | - |

**What is shared and what isn't:**
- **K69→R** (the PLP Schiff-base lysine) is shared by all three red-algal classes that can be read: *Galdieria* (Cyanidiophyceae), *Porphyra* (Bangiophyceae) and *Gracilariopsis* (Florideophyceae).
  - *Rhodosorus* has Q there, but at pLDDT 69.1, so it is unresolved.
  - *Percolomonas* has S.
- **R277→S** is shared by *Galdieria*, *Rhodosorus* and *Percolomonas*. *Gracilariopsis* has A and *Porphyra* has P.
- **Elsewhere the red-algal classes mostly differ** at the same positions:
  - E274: N / D / V / G
  - H197: S / Q / E
  - P275: V or A
- **Within *Galdieria* the substitutions are nearly identical.** The two *G. partita* copies agree at all 7 confident sites, and the three species share 4–5 of 5–6 sites.
- ***Percolomonas* shares only R277S** with any red alga.

**Reading these patterns:**
- The shared derived residues (K69R across three red classes; R277S) fit one early event in red algae. But K→R and R→S are common substitutions, so convergence can't be excluded from 13 sites.
- Different residues at the same positions do not by themselves show repeated independent decay. Once a domain has lost function, every lineage keeps drifting at those positions, so differing residues are expected whether there was one event or several.
- The tree is the test that can separate the two. It has not been run yet.

### Tree

This is pending the full run (`pip install pyfamsa`). The method:

- Align the ODC-half region (first to last PLP Pfam match) of every fusion with a region ≥ 250 aa, plus standalone references.
  - The references include the four species' own standalone ODCs and DAPDCs.
  - Distant references go in the same alignment.
- Build an NJ tree on Kimura distances with 1,000 bootstraps.
- Report support for the exact splits: red-algal fusions; red + *Percolomonas*; degraded; green; apicomplexan.

It is a distance tree, not maximum likelihood.

## 3. The biological consequence: decarboxylases in the four species

I searched UniParc for each taxon and every sequence carrying any of these Pfam families:

- PLP decarboxylase (PF02784/PF00278)
- AdoMetDC (PF01536/PF02675)
- pyruvoyl ADC (PF01862)
- agmatine iminohydrolase, AIH (PF04371)
- arginase/agmatinase (PF00491)
- nitrilase/CN-hydrolase (PF00795, the NCPAH family)

PLP decarboxylases were classified by the p6 method, with LDC and plant DAPDC references added. The 13 ODC sites were read by sequence alignment, and by structure with the pLDDT gate wherever an AlphaFold model exists.

| | *Galdieria partita* | *Rhodosorus marinus* | *Gracilariopsis chorda* | *Percolomonas cosmopolitus* |
|---|---|---|---|---|
| proteome (UniProt) | UP001061958, 7,492 proteins, BUSCO 210 complete / 13 fragmented / 32 missing of 255 | UP001157974, 8,528 proteins, BUSCO 215 / 10 / 30 | UP000247409, 10,593 proteins, BUSCO 206 / 13 / 36 | none; transcriptome (TSA) only |
| **standalone ODC-family, catalytic residues present** | **yes, 2**: A0A9C7PSP6 and A0A9C7Q3X8 (nearest human ODC, 45–47% identity; structure 11 kept / 1 lost / 1 unresolved; K69 and C360 kept; the loss is Y331→S/T) | **yes**: A0AAV8UQI0 (nearest human ODC, 48%; structure 13/13 kept, TM 0.90), plus the same protein in TSA (A0A7S2ZKI5, 13/13) | **no**: A0A2V3IYX3 is ODC-family (nearest *T. brucei* ODC, 43%) but **C360→R** (catalytic Cys) is a confident loss; structure 10/1/2 | **probably**: two TSA ODC-family sequences, 12/13 by sequence (Y331→H or F; K69 and C360 present); no model |
| DAPDC | 2 | 1 genome + TSA fragments | 2 | 1 |
| ADC (SpeA / plant-ADC fold, or pyruvoyl) | none found | none found | none found | none found |
| agmatine iminohydrolase (AIH) | 4 | none found | 1 | none found |
| arginase / agmatinase family | none found | none found | none found | 1 |
| NCPAH candidate (nitrilase family) | 1 at 58% identity to *Arabidopsis* NLP1, i.e. a probable NCPAH ortholog, plus several at 24–29% | none above 35% | 1 at 57% | none above 32% |

**Answers to the question:**
- ***Galdieria partita* and *Rhodosorus*:** each has a standalone ODC-family protein that keeps the PLP lysine and the catalytic cysteine.
  - *Rhodosorus*' standalone keeps all 13 sites, confidently modelled.
  - *Galdieria*'s standalones keep 11 of 13; the confident loss is Y331.
  - For these two species the fusion's degraded ODC half is not their only putative ODC.
- ***Gracilariopsis chorda*:** its only standalone ODC-family protein has lost the catalytic Cys360 (to Arg), and no ADC of any recognised family was found.
  - By these searches it has no protein with a complete ODC catalytic set.
  - It does have AIH and a probable NCPAH, the downstream half of the arginine route, but no ADC to feed them.
  - This is the one result that points to a real gap. It depends on a proteome that is 81% BUSCO-complete, so "not found" here is not "absent from the genome".
- ***Percolomonas*:** two transcriptome ODC-family sequences with K69 and C360. A transcriptome can show presence; it cannot show absence.

**Limits of the classification:**
- Nearest-reference similarity cannot tell an ODC from an ODC-fold arginine decarboxylase. The PBCV-1 ADC is an ODC-fold enzyme, and the plant LDCs derive from ODC.
- So an "ODC-family" enzyme here is not proven to use ornithine.
- No SpeA- or plant-ADC-fold enzyme was found in any of the four. If these algae decarboxylate arginine, it would have to be by an ODC-family enzyme or by a gene missing from these proteomes.

### Literature: measured polyamine pathways in red algae

Europe PMC and PubMed, 20 Nov 2026; the four queries and their hit counts (27, 198, 5 and 24) are appended to `literature_query_log.tsv`. Abstracts were read for all of the following; none was read in full text.

- **Guo et al. 2026, *Plant Physiol Biochem* (PMID 42679707), *Pyropia yezoensis*.**
  - ODC, ADC and N-carbamoylputrescine amidohydrolase transcripts rise during conchospore development.
  - The ODC inhibitor DFMO and siRNA silencing of *PyODC* block development; putrescine partly rescues it.
  - So the ODC route is functional and needed in this red alga.
  - *P. yezoensis* is also one of the new fusion species (A0ACC3CDY0). Whether *PyODC* is the fusion or a standalone gene can't be told from the abstract, and the paper is not open access.
- **Schweikert et al. 2014, *J Phycol* (PMID 26988184), *Pyropia cinnamomea*.** Both ADC and ODC activities were measured. ADC activity rose under UV-B stress; ODC activity correlated with putrescine under low stress. Activity was assigned with irreversible inhibitors.
- **García-Jiménez et al. 2009, *J Plant Physiol* (PMID 19524320), *Grateloupia imbricata*.** The ODC gene was cloned (GenBank FJ223132), and the authors state that "ODC is predominant in the synthesis of polyamines". I did not check whether FJ223132 is standalone or fused.
- **Hamana & Matsuzaki 1982, 1985 (PMIDs 7096289, 4030724).** Putrescine and spermidine were found in all red algae examined; spermine and sym-homospermidine only in Rhodophyta.
- **Hamana & Niitsu 2006 (PMID 17116972).** Cellular polyamines of Glaucophyta, Rhodophyta, Cryptophyta, Haptophyta and **Percolozoa** (= Heterolobosea). It has no abstract in either database, and I could not read it.
- **Heterolobosea/*Naegleria* + polyamine enzymes:** 5 hits, none measuring a pathway.

**Bottom line from the literature:** red algae do make putrescine. In *Pyropia* both ADC and ODC activities are measured and ODC is required. None of these papers mentions a fused AdoMetDC–ODC.

## 4. Retained or rotting

| ODC status | n | partner length (aa) | TM to human ODC | longest unaligned stretch | mean pLDDT | % partner < 70 |
|---|---|---|---|---|---|---|
| degraded | 8 | 354–435 (median 372) | 0.65–0.82 (0.71) | 14–35 (22.5) | 65–79 (75) | 19–42 (30) |
| intact | 12 | 388–484 (400) | 0.84–0.89 (0.87) | 13–29 (17) | 81–95 (90) | 2–22 (8) |
| undetermined | 4 | 327–489 (380) | 0.69–0.87 (0.85) | 17–65 (20) | 81–89 (85) | 7–24 (15) |
| incomplete | 3 | 65–254 (143) | 0.15–0.53 (0.30) | 19–147 (43) | 65–82 (80) | 9–61 (20) |

**Reading the table:**
- **The degraded domains are retained, not truncated.** They are close to full length: 346–379 of 422 human-ODC positions align, and the longest gap (≤ 35) is no larger than in intact domains.
- **But they are less ODC-like and less confidently modelled.** TM is 0.65–0.82 against 0.84–0.89, and a median 30% of residues sit below pLDDT 70, against 8%.
- **This is the opposite of the green incomplete ones**, which are short or gapped.
- **Lower pLDDT has more than one possible cause.** It can reflect real disorder, or simply fewer close homologs in AlphaFold's alignment for a diverged sequence. These data cannot separate the two.

## 5. Correlates across the 23 (descriptive only)

| variable | degraded (n = 4) | the other 19 |
|---|---|---|
| AdoMetDC half (p9) | 4 intact | 15 intact, 1 not intact (*Riccia*), 3 undetermined |
| linker length (aa) | 71, 80, 114, 124 | 41–202 |
| protein length (aa) | 730–809 | 411–1039 |
| gene-model source | 3 genome (1, 5 and 9 exons), 1 TSA | 17 genome, 2 TSA |
| putrescine-site residues 13 / 111 / 285 | N, R or K / R, Q, S or I / F, H, V or Y | almost all R / K / S |

**What the table shows:**
- **Most variables give no separation.** AdoMetDC integrity, linker length, protein length and gene-model source overlap completely between the degraded four and the rest.
- **The putrescine-site residues do differ.** The degraded four depart from the green consensus (R13, K111, S285) at 111 and 285.
- **This is not an independent correlate.** These four are also the only non-green lineages except *Vitrella*, so any residue difference between them and the green fusions is confounded with phylogeny. n = 4 can't separate the two.

## What is too small to support a conclusion

- **Degradation in red algae:** the claim that degradation is red-algal-specific rests on 7 red-algal fusions with models from 4 genera, and 1 heterolobosean.
  - The dinoflagellate, glaucophyte, CRuMs and brown-algal fusions are needed to test whether degradation is red-specific or simply non-green. The full run checks those with models.
  - Many of them, mostly TSA-derived, have no AlphaFold model and will get sequence reads only.
- **One event or many:** this can't be decided from substitutions alone (see section 2). The tree is the test.
- **Gracilariopsis lacking a working ODC:** this is one species and rests on an incomplete proteome.

## Run

```
py -3.12 -m pip install pyfamsa
cd scripts
py -3.12 p10_red_algal_pattern.py
```

- **What it does:** analyses the remaining new fusions, downloading their AlphaFold models to `data/p10_models/` (never Marc's folder), then builds the tree. Caches are reused.
- **How long:** expect tens of minutes. Most of the time goes on downloads and the 1,000 bootstraps.
- **What it writes:** `p10_report.md`, `p10_odc_tree.nwk`, `p10_tree_support.tsv` and `p10_odc_alignment.fasta`, and it overwrites the other `p10_*.tsv` with the complete set.
