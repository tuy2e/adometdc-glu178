# Round 4 (10 Nov 2026): the BLAST sampling added to both trees

**Short version**

- **Spermidine synthase (245 sequences, up from 143).** With 12 photosynthetic cryptophytes instead of 2, **kinetoplastids and cryptophytes now form a clade (SH-aLRT 93.5 / UFBoot 90)**. That is the same pairing as AdoMetDC group (b), so the partner gene now echoes it clearly.
  - Euglenids don't join that pairing: they group with diplonemids (96.7/100).
  - The colourless cryptophyte *Goniomonas* sits outside the other cryptophytes.
  - Group (a) is still rejected (AU p = 0.001), and so is the 4-group lineage (p = 0.00006).
  - Forcing euglenids into the cryptophyte–kinetoplastid group is rejected too (p = 0.026).
- **AdoMetDC (200 enzymes, up from 188).**
  - The 4 new euglenid sequences (*Euglena longa*, *Rhabdomonas costata*, a second *Eutreptiella* copy, plus *E. gracilis* from round 3) all sit with kinetoplastids inside group (b): euglenids 95.1/100, + kinetoplastids 99.6/100, + cryptophytes 99/100.
  - **All 8 diplonemids keep Glu178** and group with heteroloboseans, away from kinetoplastids. The Euglenozoa split in AdoMetDC now rests on 8 diplonemids instead of 1. Forcing all Euglenozoa together is rejected (AU p = 0.017).
  - Both centrohelids keep Glu178.
- **Glu178 is still lost once and never regained.**
  - One wrinkle: in the best tree, the 4 fast-evolving *Blastocystis* sequences jump out of the lineage to sit with heteroloboseans (UFBoot 48), a typical long-branch artifact. That splits the lineage and makes the ancestral reconstruction count a second, *Blastocystis*-only loss.
  - Putting *Blastocystis* back is not rejected (p = 0.27), and on that tree there is a single Glu→Asn at the lineage base (P 0.994).
  - Support for the lineage without *Blastocystis* drops to UFBoot 65 (was 90).
- **New residue facts.**
  - *E. longa* matches *E. gracilis* exactly at the pocket (R13, **Q15**, Y109, D174, Q176, **N178**). The Q15 is shared by the two *Euglena* species, not unique to *E. gracilis*.
  - ***Rhabdomonas* has Lys at 178**: a positive charge where human has Glu. That is reminiscent of plant AdoMetDCs, which have built-in Arg/Lys residues that stand in for putrescine. This is a hypothesis, not a result.

---

## What went in

All hit sequences came from the 13 BLAST reports in `3_data/blast_results/` (Download All → Text). Each hit's aligned stretch was read from the `Sbjct` lines. When a hit had several HSPs that didn't overlap on the query, they were joined in query order, and stop codons became X. No sequences were fetched separately.

- **Check:** the *E. gracilis* AdoMetDC (search G) and SpdSyn (search H) HSPs read this way are exact substrings of the full ORFs used in round 3.
- **Limit:** for transcriptome hits this gives the aligned region (the domain), not the full ORF.

### Spermidine synthase

| Step | Result |
|---|---|
| Candidates from A, A2, B, C, D2, D3 (≥ 150 aa; *E. gracilis* taken from round 3's full ORF) | 598 |
| Collapse near-identical copies (CD-HIT 97%) | 314; 1 already in the reference set |
| Classify on the full aminopropyltransferase family tree (the 375-sequence `c1` set, FastTree). A hit counts as SpdSyn if its nearest reference and at least 2 of its 3 nearest references are in the 9 Nov SpdSyn clade | 163 SpdSyn; the rest are mostly thermospermine synthase-like (105) and other PABS-domain proteins |
| Contamination screen (same rule as before: drop if the closest relative is in another supergroup, ≥ 0.50 identity and closer than anything in its own group), plus the dinotoms *Durinskia* and *Kryptoperidinium* | 21 dropped, 142 kept. Examples: *Teleaulax* ↔ *Phaeocystis* pairs at 90% (cross-contamination between a cryptophyte and a haptophyte transcriptome); *Bolivina plicata* 100% identical to a stramenopile; *Sorites* 90% to a dinoflagellate (its algal symbiont) |
| Redundancy reduction (90% within group; the 142 original sequences always kept) | **102 new**: 12 cryptophytes (11 species), 1 *Roombia*, 3 euglenids, 6 diplonemids, 18 dinoflagellates/*Perkinsus*, 42 Rhizaria, 21 haptophytes |
| Final set | **245** = 142 original + *E. gracilis* + 102 new |

- **Alignment:** MAFFT L-INS-i on all 245 from scratch (a new dataset, so a fresh alignment), then ClipKIT gappy 0.5, giving **286 columns**.
- **Tree:** IQ-TREE 2.0.7, ModelFinder (LG/WAG/JTT/Q.pfam × G4/R4–R7), which again picked **LG+R6**. UFBoot 1000 + SH-aLRT 1000.

### AdoMetDC

From searches F and E: one sequence per species (the longest, highest-scoring HSP), trimmed to the human-aligned domain.

- ***Diplonema papillatum*** (GJNJ01033060) is 100% identical to the one already in the tree, so it was not added.
- ***Eutreptiella*** HBJA01101101 is only **69% identical** to the *Eutreptiella* sequence already in the tree (UPI0016DBF433). That is too different for strain variation, so it is probably a second gene copy. It was added as "copy 2".
- The contamination check was passed by all 12: each new sequence's closest relative is in its own group, and the centrohelids' closest foreign match is 0.43.
- **Alignment:** added to the round-3 alignment with `mafft --add --keeplength`, keeping the original 340-column mask, so all 188 earlier rows are unchanged. Result: **200 × 340**.
- **Tree:** IQ-TREE Q.pfam+R7, UFBoot 1000 + SH-aLRT 1000, with ancestral reconstruction.

## Spermidine synthase tree

Figure: `2_figures/Fig_spdsyn_tree_round4.png`

| Clade | 9 Nov (142) | Round 3 (+ *E. gracilis*, 143) | Round 4 (245) |
|---|---|---|---|
| Kinetoplastids | 100 / 98 | 100 / 50 | 100 / 100 |
| Cryptophytes (all, incl. *Goniomonas*) | 100 / 100 (2 seqs) | 100 / 100 | not monophyletic: *Goniomonas* (2) sits apart |
| **Kinetoplastids + 12 photosynthetic cryptophytes** | not in tree (AU p 0.51) | not in tree (p 0.48) | **93.5 / 90** |
| Euglenids | – | 1 seq | 100 / 100 (3) |
| Euglenids + diplonemids | – | 99.2 / 100 (2 seqs) | 96.7 / 100 (3 + 7) |
| All Euglenozoa | – | 92.2 / 41 | not in tree. Euglenids + diplonemids (+ *Tritrichomonas*, 19.7/53) are sister to kinetoplastids + cryptophytes (64.7/89) |
| Stramenopiles + dinoflagellates (group a) | not in tree | not in tree | not in tree |
| 4-group lineage | not in tree | not in tree | not in tree |

- *Roombia*, a katablepharid (a single fragment), sits near apicomplexans and *Paulinella* with weak support, not with cryptophytes.
- Some foram (Rhizaria) hits cluster with a *Helicobacter pylori* sequence near the base of the tree. These may come from bacteria living in or on the forams. They don't touch the groups tested here.

**AU tests** (10,000 RELL, LG+R6):

| Tree (constraint) | ΔlnL | p-AU | 9 Nov (142 seqs) |
|---|---|---|---|
| ML (kinetoplastids + 12 photosynthetic cryptophytes, 90) | 0 | 0.637 | – |
| Cryptophytes (all 14, incl. *Goniomonas*) + kinetoplastids | 9.5 | **0.457, not rejected** | 0.51 |
| Cryptophytes + all Euglenozoa (kinetoplastids + euglenids + diplonemids) | 8.8 | 0.485, not rejected | – |
| Cryptophytes + kinetoplastids + euglenids (AdoMetDC group b exactly) | 63.2 | **0.026, rejected** | – |
| Cryptophytes + red algae | 75.1 | 0.044, rejected (borderline; shorter search, see note) | 0.016 |
| Stramenopiles + dinoflagellates/*Perkinsus* (group a) | 138.1 | **0.001, rejected** | 0.002 |
| 4-group lineage (stramenopiles + dinos + kinetoplastids + euglenids + cryptophytes) | 242.2 | **0.00006, rejected** | 0.0006 |

*Note:* the last three constrained trees were searched with `-nstop 30` instead of the default 100 to save time. A shorter search can only make a constrained tree look worse, so those three rejections are slightly less conservative. The first three used the full search.

## AdoMetDC tree

Figure: `2_figures/Fig_tree_Glu178_round4.png`

| Clade | Round 3 (188) SH/UFB | Round 4 (200) SH/UFB |
|---|---|---|
| Euglenids (5 sequences, 4 species) | 100 / 100 (2) | 95.1 / 100 |
| *E. gracilis* + *E. longa* | – | 100 / 100 |
| Kinetoplastids + euglenids | 99.9 / 100 | 99.6 / 100 |
| Group (b): + cryptophytes | 99.5 / 100 | 99 / 100 |
| Cryptophytes | 81.1 / 98 | 90.1 / 99 |
| Dino/*Perkinsus* Asp copy | 100 / 100 | 99.9 / 100 |
| Stramenopiles | 65.6 / 74 | 94.1 / 94 without *Blastocystis* (with it: not a clade) |
| 4-group lineage | 67.3 / 90 | 95.7 / 65 without *Blastocystis* (with it: not a clade) |
| Diplonemids (8) | 1 seq | 100 / 100; with 4 heteroloboseans + 1 amoebozoan (82 / 44) |
| Centrohelids (2) | – | 95.3 / 89; near red algae, unsupported |

**AU test** (10,000 RELL, Q.pfam+R7):

| Tree | ΔlnL | p-AU |
|---|---|---|
| ML tree (*Blastocystis* outside the lineage) | 0 | 0.772 |
| Whole lineage forced together, *Blastocystis* included | 17.5 | **0.273 (not rejected)** |
| All Euglenozoa forced together (kinetoplastids + 5 euglenids + 8 diplonemids), rest of the lineage together | 52.9 | **0.017 (rejected)** |

**Residues** (multiple alignment, human numbering):

| | 13 | 15 | 109 | 111 | 174 | 176 | 178 | 256 |
|---|---|---|---|---|---|---|---|---|
| human | L | E | S | F | D | T | **E** | E |
| *Euglena gracilis* | R | Q | Y | Q | D | Q | **N** | E |
| *Euglena longa* | R | Q | Y | Q | D | Q | **N** | E |
| *Rhabdomonas costata* | R | E | F | S | D | Q | **K** | E |
| *Eutreptiella* (both copies) | R | E | Y | Q | D | Q | **N** | E |
| Diplonemids (7 of 8) | K | E | Y | S | D | T | **E** | E |
| *Lacrimia lanifica* | R | E | F | S | D | T | **E** | E |
| Centrohelids (2) | T | V | F | H/N | Q/E | T/V | **E** | E |

**Ancestral states** (best tree, posterior probabilities):

- Euglenid ancestor: R13, E15, Y109 (0.61), D174, Q176, N178 (0.99).
- *E. gracilis* + *E. longa* ancestor: Q15 (1.00), Y109, N178.
- Kinetoplastid + euglenid ancestor: R13, W109 (0.94), D174, Q176, N178 (0.99).
- Group (b) ancestor: R13, W109, D174, N178 (1.00).
- Diplonemid ancestor: K13, Y109, D174, T176, **E178 (1.00)**.
- On the tree with *Blastocystis* back in place: **one** E(0.999)→N(0.994) at the lineage base. Then N→D in the dino copy, N→H in oomycetes, N→S in kinetoplastids and N→K in *Rhabdomonas*. **No return to Glu anywhere.**

## What it means

1. **Group (b) is now echoed by the partner gene.** In SpdSyn, kinetoplastids sit with the photosynthetic cryptophytes (UFBoot 90). Before, that pairing was only "not rejected" and rested on 2 cryptophytes.
   - So two genes of the same pathway now pair cryptophytes with kinetoplastids.
   - Li 2015's neighbour-joining spermidine synthase tree showed the same pairing.
   - That points to a shared history for this part of the pathway (gene transfer between these lineages, or an ancient shared origin), or to an artifact common to both genes. It no longer looks like an AdoMetDC-only quirk.
2. **Euglenids differ between the two genes.** In AdoMetDC they sit inside group (b) with kinetoplastids. In SpdSyn they sit with diplonemids. So whatever moved the cryptophyte/kinetoplastid genes did not involve the euglenid SpdSyn.
3. **Diplonemids are firmly outside the Glu-less group** (8 of 8 keep Glu178). Kinetoplastids + euglenids lost Glu178; their euglenozoan relatives the diplonemids didn't.
4. **Group (a) and the 4-group lineage still have no echo in SpdSyn.** For stramenopiles + the dino copy, the AdoMetDC pattern is gene-specific.
5. **The prediction stays "stramenopiles only".** *E. longa* matches *E. gracilis*, so euglenids look alike. *Rhabdomonas*'s Lys178 is a new oddity worth one sentence in the paper.

## Limits

- The new sequences are BLAST HSP regions from transcriptomes, not full ORFs. Ends are missing (2–20% gaps in the trimmed alignments), and stops were coded as X.
- The family-tree classification uses FastTree and nearest references: good enough to sort sequences into SpdSyn and the rest, not to name the others.
- Transcriptome contamination is screened with a fixed rule but not ruled out. The foram/*Helicobacter* cluster and one *Lingulaulax* sequence among diatoms show how much noise this data carries.
- AdoMetDC support for the lineage dropped (UFBoot 90 → 65, without *Blastocystis*), and *Blastocystis* moves. Both are typical of adding fast-evolving taxa. The AU test shows the intact lineage is still a fine explanation.
- The cryptophyte–kinetoplastid clade in spermidine synthase (UFBoot 90) is a single gene's signal, like AdoMetDC's. Two genes agreeing is stronger than one, but both could share the same artifact. The next test would be a third pathway gene, or a check for compositional bias in these two.

## Files

- **Trees:** `3_data/trees/main_tree_round4_200.treefile`, `spermidine_synthase_tree_round4_245.treefile`
- **Alignments:** `3_data/alignments/adometdc_round4_200_trimmed.aln`, `spermidine_synthase_round4_245_trimmed.aln`
- **Sequences:** `3_data/sequences/round4_new_adometdc.fa` (12), `round4_new_spermidine_synthase.fa` (102), `round4_blast_hits_all.json` (all 839 parsed hits)
- **Reports:** `3_data/test_results/main_tree_round4_report.iqtree`, `AU_main_tree_round4.iqtree`, `spermidine_synthase_round4_report.iqtree`, `AU_spermidine_synthase_round4.iqtree`, `ancestral_main_tree_round4_site178.json`
- **Figures:** `2_figures/Fig_tree_Glu178_round4.png`, `2_figures/Fig_spdsyn_tree_round4.png`
- **Scripts:** `4_scripts/round4_*.py`
- **Everything else:** `8_full_workspace_backup/workspace_part4_round4.zip`
