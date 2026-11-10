# *Euglena gracilis* added to the AdoMetDC and spermidine synthase trees (10 Nov 2026)

**Short version**

- **AdoMetDC:** the *E. gracilis* enzyme sits with *Eutreptiella* (100/100), then with kinetoplastids (99.9/100), then with cryptophytes (group b, 99.5/100). The 4-group Glu-less lineage stays together (UFBoot 90). Glu178 is never regained.
- **Spermidine synthase:** the *E. gracilis* SpdSyn sits with *Diplonema* (99.2/100), then with kinetoplastids. This is the normal Euglenozoa grouping, as in the species tree. Forcing *Euglena* into the exact AdoMetDC group (b), cryptophytes + kinetoplastids + *Euglena*, is rejected: ΔlnL 34.9, AU p 0.016 (0.049 in the 5-tree set that matches 9 Nov). Cryptophytes + kinetoplastids with *Euglena* left free is still not rejected (p 0.48), and neither is cryptophytes + all Euglenozoa (p 0.47).
- **So:** the two genes now disagree about euglenids. AdoMetDC splits Euglenozoa (euglenids + kinetoplastids go with cryptophytes, and *Diplonema* goes with heteroloboseans). SpdSyn keeps Euglenozoa together. This is more evidence that the AdoMetDC pattern belongs to that gene alone.
- **Residues:** the *Euglena* AdoMetDC has R13, **Q15**, Y109, D174, **Q176**, **N178** and E256. It is putrescine-activated (Aleksijevic 1979), yet it lacks three of the residues usually listed for the putrescine site (Glu15, Thr176, Glu178). It also has Tyr, not Trp, at 109. The ancestral reconstruction puts the switch to Tyr109 at the euglenid ancestor.
- **One caution:** the order of branches inside the lineage changed. The new best tree joins stramenopiles to group (b), with the dino/*Perkinsus* Asp copy branching first. The 9 Nov tree had stramenopiles + Asp copy together. The two arrangements are equally good: ΔlnL 0.7, AU p 0.58 vs 0.52. On the new arrangement the first residue after Glu is uncertain (N 0.41 / E 0.30 / D 0.29). On the old arrangement it is still a clean E→N (0.996).

---

## Sequences

All 5 proteins come from `3_data/sequences/euglena_gracilis_from_tblastn.fa`. They match the backup copy `audit/euglena/euglena_proteins.fa` exactly.

| ID | Search | Used in | Notes |
|---|---|---|---|
| Egr_AdoMetDC_GDJR01034331 | G (F5EJKSE8014) | AdoMetDC tree | Longest ORF; the N-terminus is uncertain. The human-anchored domain trim kept the whole 401 aa. |
| Egr_AdoMetDC_GEFR01006935 | G | not used | Same gene from the second assembly, with a shorter N-terminus. |
| Egr_SpdSyn_GDJR01022181 | H (F5FJXSNR014) | SpdSyn tree | Full 285-aa ORF. |
| Egr_APT2_GDJR01022779 | H | classified only | Falls in the **thermospermine synthase (ACL5-like)** group, beside dinoflagellate and haptophyte thermospermine synthases (SH-like 0.995). |
| Egr_APT3_GEFR01015322 | H | classified only | Falls in a small, divergent PABS-domain group (13 tips: haptophytes, dinoflagellates, *Pelagomonas*; SH-like 0.98), outside the SpdSyn clade. Identity unknown. |

How the APT hits were classified: I added them to the full 375-sequence aminopropyltransferase set (`audit/spdsyn/c1`: all of PF01564, including thermospermine and spermine synthases), then ran MAFFT `--add`, ClipKIT gappy 0.5 (289 columns) and FastTree `-lg`. Only the true SpdSyn falls inside the SpdSyn clade, which was the 9 Nov rule for including a sequence, so APT2 and APT3 stay out of the SpdSyn tree.

## Method, and one change from the original pipeline

- The tools match the earlier work: IQ-TREE 2.0.7, the same models and the same commands as `expA` and `spdA`/`C_h*`/`au_spd`, with seed 1.
- **I did not realign from scratch.** MAFFT v7.505 here doesn't reproduce the earlier alignments: realigning the original 142 SpdSyn sequences gives 902 columns, not 866. A fresh realignment would therefore mix two changes, the new sequence and the new aligner version.
- Instead I added each *Euglena* sequence to the **original untrimmed alignment** with `mafft --localpair --maxiterate 1000 --add --keeplength`. Then I applied the **original ClipKIT column mask**. I checked that ClipKIT gappy 0.5 on the original alignments gives the original trimmed files exactly.
- `--keeplength` removed 32 *Euglena* AdoMetDC residues (positions 4–29 and 37–42 of the uncertain N-terminal extension). No other sequence has residues in those columns, so ClipKIT would have trimmed them anyway. No SpdSyn residues were removed.
- Result: the trimmed alignments are the 9 Nov ones (byte-identical rows) plus one *Euglena* row. That gives 188 × 340 for AdoMetDC and 143 × 292 for SpdSyn. Any change in the trees comes from *Euglena*.

## AdoMetDC tree (188 enzymes, Q.pfam+R7, UFBoot 1000 / SH-aLRT 1000)

The figure is `2_figures/Fig_tree_Glu178_expanded_with_Euglena.png`.

| Clade | 9 Nov (expA, 187) SH/UFB | Now (expE, 188) SH/UFB |
|---|---|---|
| *Euglena* + *Eutreptiella* | (n/a) | 100 / 100 |
| Kinetoplastids + euglenids | 99.3 / 100 | 99.9 / 100 |
| Group (b): kinetoplastids + euglenids + cryptophytes | 99.4 / 100 | 99.5 / 100 |
| Cryptophytes | 78 / 97 | 81.1 / 98 |
| Stramenopiles | 56.8 / 68 | 65.6 / 74 |
| Dino/*Perkinsus* Asp copy | 99.9 / 100 | 100 / 100 |
| Stramenopiles + Asp copy (group a) | 7.2 / 56 | not in the ML tree |
| Stramenopiles + group (b) | not in the ML tree | 15.8 / 49 |
| 4-group Glu-less lineage | 93.5 / 88 | 67.3 / 90 |

*Diplonema* still keeps Glu178 and still sits with heteroloboseans (61.2/98).

**AU test** (10,000 RELL, Q.pfam+R7, 3 trees):

| Tree | ΔlnL | p-AU |
|---|---|---|
| New ML tree (stramenopiles + group b) | 0 | 0.578 |
| 9 Nov arrangement, *Euglena* grafted next to *Eutreptiella* | 0.71 | 0.518 (equally good) |
| H_kin + *Euglena*: kinetoplastids + both euglenids + *Diplonema* forced together (Euglenozoa), rest of lineage together | 48.3 | **0.019 (rejected)** |

*(The 9 Nov H_kin had ΔlnL 54.6 and p 0.001 / 0.0005 in the 5- and 7-tree sets. The p-values aren't directly comparable because the sets of trees differ.)*

**Residues** (human numbering, multiple alignment):

| | 13 | 15 | 109 | 111 | 174 | 176 | 178 | 256 |
|---|---|---|---|---|---|---|---|---|
| human | L | E | S | F | D | T | **E** | E |
| *E. gracilis* | R | **Q** | Y | Q | D | Q | **N** | E |
| *Eutreptiella* | R | E | Y | Q | D | Q | N | E |
| kinetoplastids (6/6) | R | E | W | S | D | Q | S | E |
| cryptophytes (7) | K (R in *Goniomonas*) | E | W/Y | Q | D/R/E | M/S/T | N | E |

**Ancestral reconstruction**, new ML tree, posterior probabilities:

- Euglenid ancestor: R13 1.00, E15 0.97, **Y109 0.98**, D174 1.00, Q176 1.00, N178 1.00. *Euglena*'s Q15 is a change on its own tip branch.
- Kinetoplastid + euglenid ancestor: R13 1.00, W109 0.99, D174 1.00, Q176 0.99, N178 0.99. This is the same inherited R13/D174/N178 pocket reported on 10 Nov.
- Group (b) ancestor: R13 0.99, W109 1.00, D174 1.00, N178 1.00.
- Kinetoplastid ancestor: S178 0.98. **Ser178 is shared by all 6 kinetoplastids**, so it isn't specific to *T. cruzi*. The log's line "*T. cruzi* later changed Asn→Ser" is better read as "kinetoplastids changed Asn→Ser".
- Base of the 4-group lineage: Glu lost (parent E 0.997), but the first replacement is uncertain (N 0.41 / E 0.30 / D 0.29). After that, N is 1.00 at the stramenopile + group (b) node, and the Asp copy goes to D (0.97).
- On the 9 Nov arrangement with *Euglena* added: E(0.999)→N(0.996) at the lineage base, then N→D in the Asp copy, N→S in kinetoplastids and N→H in oomycetes.
- **Glu178 is never regained on either tree.**

## Spermidine synthase tree (143 sequences, LG+R6, UFBoot 1000 / SH-aLRT 1000)

- *Euglena* SpdSyn + *Diplonema*: **99.2 / 100**
- That pair + kinetoplastids (Euglenozoa): **92.2 / 41**
- Kinetoplastids: 100 / 50 (was 100/98; *Euglena* and *Diplonema* attach at different points across the bootstrap trees)
- Cryptophytes: 100 / 100. They are not next to Euglenozoa; their sister is a large group (63.9/51).
- On 9 Nov, kinetoplastids were sister to *Tritrichomonas* + *Diplonema* (62.9/50). Now *Diplonema* moves to *Euglena*.

**AU test** (10,000 RELL, LG+R6; rows 1–5 are the same 5 hypotheses as on 9 Nov, with *Euglena* placed by its AdoMetDC group):

| # | Tree (constraint) | ΔlnL | p-AU, 7-tree set | p-AU, 5-tree set (rows 1–5 only) | 9 Nov (no *Euglena*) |
|---|---|---|---|---|---|
| 1 | ML | 0 | 0.604 | 0.970 | 0.538 |
| 2 | Cryptophytes + kinetoplastids + *Euglena* (AdoMetDC group b exactly) | 34.9 | **0.016 (rejected)** | 0.049 | CR+KE: Δ0.25, p 0.51 |
| 3 | Cryptophytes + red algae | 37.8 | 0.020 (rejected) | 0.058 | 0.016 |
| 4 | Stramenopiles + dino/*Perkinsus* (group a) | 107.5 | 0.0003 (rejected) | 0.0015 | 0.002 |
| 5 | 4-group lineage (incl. *Euglena*) | 189.0 | 3e-8 (rejected) | <1e-10 | 0.0006 |
| 6 | Cryptophytes + kinetoplastids, *Euglena* left free *(new)* | 4.6 | 0.484 (not rejected) | – | – |
| 7 | Cryptophytes + all Euglenozoa (kinetoplastids + *Euglena* + *Diplonema*) *(new)* | 3.4 | 0.466 (not rejected) | – | – |

How to read this:

- SpdSyn still can't rule out cryptophytes sitting next to kinetoplastids, or next to Euglenozoa as a whole (rows 6–7, ΔlnL 3–5). So the weak echo of the cryptophyte + Discoba pairing remains.
- But the *Euglena* SpdSyn won't join the exact AdoMetDC group (b) (row 2). It stays with *Diplonema*.
- The 9 Nov rejections of group (a) and of the 4-group lineage still hold, and are stronger now.

## What it means

1. *Euglena* doesn't change the AdoMetDC picture. Group (b) and the cross-supergroup lineage hold, with the same support. Moving kinetoplastids + euglenids back into Euglenozoa is still rejected.
2. The SpdSyn gap is now filled, and it points the other way for euglenids: SpdSyn follows the species tree. So the AdoMetDC grouping of euglenids with cryptophytes is not shared by the partner gene.
3. The one-vs-two-losses question is unchanged. The order of branches inside the lineage was already weak (SH 7.2 on 9 Nov), and the two arrangements are tied. Report it as unresolved, as before.
4. Biochemistry: *Euglena* is activated by putrescine while lacking Glu15, Thr176 and Glu178, and it keeps Asp174. This supports keeping the clean prediction to **stramenopiles only**. It also fits a euglenozoan way of being activated that differs from the human one, but that is a hypothesis, not a result.

## Limits

- The *Euglena* sequences come from transcriptome assemblies. The AdoMetDC N-terminus is uncertain, but the trimmed alignment covers the domain (4% gaps).
- APT2 and APT3 were classified with FastTree SH-like supports on a broad family tree. That is good enough to keep them out of the SpdSyn tree, not for naming APT3.
- Of the AdoMetDC AU hypotheses, only H_kin was rerun with *Euglena*, because it is the one *Euglena* bears on. H_crypto, H_dino, H_all, H_li/H_sar and the C60 runs were not rerun.

## Files

- **Alignments:** `3_data/alignments/adometdc_main_trimmed_plus_Euglena.aln`, `spermidine_synthase_trimmed_plus_Euglena.aln`
- **Trees:** `3_data/trees/main_tree_Qpfam_R7_plus_Euglena.treefile` (+ `.contree`), `spermidine_synthase_tree_plus_Euglena.treefile`, `aminopropyltransferase_family_FastTree_plus_Euglena.nwk`
- **Reports:** `3_data/test_results/main_tree_plus_Euglena_report.iqtree`, `AU_main_tree_plus_Euglena.iqtree`, `spermidine_synthase_plus_Euglena_report.iqtree`, `AU_spermidine_synthase_plus_Euglena.iqtree`, `ancestral_main_tree_plus_Euglena_site178.json`
- **Figure:** `2_figures/Fig_tree_Glu178_expanded_with_Euglena.png`
- **Scripts:** `4_scripts/euglena_*.py`
- **Everything else** (constraint trees, constrained-search outputs, logs, `.state` files): `8_full_workspace_backup/workspace_part3_euglena_trees.zip`
