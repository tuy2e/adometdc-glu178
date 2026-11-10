# Spermidine synthase (SpdSyn) tree, 9 Nov 2026
Data: UniProt PF01564 for Stramenopiles, Alveolata, Cryptophyceae, Discoba, Haptista, Rhizaria, Rhodophyta, Amoebozoa, Metamonada, Apusozoa (all entries); Viridiplantae, Opisthokonta and prokaryotes (reviewed only).
- NCBI name search added Baffinella frigidus (cryptophyte), Diplonema papillatum, Acrasis kona, Naegleria, Plasmodiophora, Emiliania.
- NCBI BLAST (5 jobs, to catch unannotated transcriptome proteins) stayed "WAITING" for more than 2.5 h and never returned. So cryptophytes = 2, euglenids = 0, dinoflagellates = 3, Perkinsus = 2.
Filtering:
- single-domain (PF17284/PF01564 only), 250–450 aa, no methyltransferases;
- SpdSyn clade defined by reference enzymes (human SRM P19623, yeast SPE3 Q12074, Arabidopsis SPDS1/2, Dictyostelium, L. major, T. brucei, P. falciparum SpdSyn), excluding thermospermine synthase (ACL5 Q9S7X6) and animal spermine synthase (P52788);
- contamination rule as for AdoMetDC (4 dropped: Vitrella A0A0G4F525, Porphyridium A0A5J4Z939, Cardiosporidium A0ABQ7JBP6, Acanthamoeba L8GVA0);
- 90% redundancy reduction.
=> 142 sequences (130 eukaryotic + 12 prokaryotic), MAFFT L-INS-i, ClipKIT gappy 0.5 -> 292 columns. IQ-TREE 2: LG+R6 (BIC), UFBoot 1000 + SH-aLRT 1000. logL -39948.583.
ML tree (rooted on prokaryotes):
- Kinetoplastids monophyletic (100/98), grouped with Tritrichomonas + Diplonema (62.9/50).
- Cryptophytes (Guillardia + Baffinella, 100/100) sit next to 11 prokaryotes (83.6/53). By identity they are typical eukaryotic SpdSyn (55% to red algae, human, diatoms).
- Stramenopiles NOT monophyletic: (i) a red-alga + diatom/ochrophyte clade (97.8/81); (ii) a clade of oomycetes, Bigyra, other ochrophytes + haptophytes + Plasmodiophora (78/34).
- Dinoflagellates (3) and Perkinsus (2) branch separately; apicomplexans at the base.
AU tests (10,000 RELL, LG+R6):
| Tree | dlnL | p-AU |
| ML | 0 | 0.538 |
| Cryptophytes + kinetoplastids clade (AdoMetDC group b) | 0.25 | 0.51 (not rejected; equally good) |
| Cryptophytes + red algae clade | 35.5 | 0.016 (rejected) |
| Stramenopiles + dino/Perkinsus clade (AdoMetDC group a) | 100.2 | 0.002 (rejected) |
| Stramenopiles + dino/Perkinsus + kinetoplastids + cryptophytes clade (AdoMetDC 4-group lineage) | 117.6 | 0.0006 (rejected) |
Interpretation:
- SpdSyn does not share AdoMetDC's group (a) or its 4-group lineage, so the AdoMetDC pattern is gene-specific, not a whole-pathway or whole-genome artifact.
- SpdSyn can't separate cryptophytes+kinetoplastids from the ML placement; Li 2015's NJ also paired Guillardia with kinetoplastids. So that pairing is echoed weakly in a second gene.
- One vs two AdoMetDC lineages: not settled by SpdSyn.
