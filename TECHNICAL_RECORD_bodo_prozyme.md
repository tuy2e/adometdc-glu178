# HARDENED FINDING: Bodo saltans has an active+dead (prozyme-like) AdoMetDC tandem pair
Date: 3 Nov 2026. Sources: UniProt (Bodo saltans genome, Jackson 2016), MAFFT alignment vs human/Tb/Tc.

## Gene-level (UniProt, verifiable accessions)
- A0A0S4JHX8  "S-adenosylmethionine decarboxylase proenzyme" (EC 4.1.1.50)  gene BSAL_30285  381 aa  = ACTIVE
- A0A0S4JRG8  "S-adenosylmethionine decarboxylase, putative"               gene BSAL_30280  331 aa  = DEAD/prozyme-like
- Consecutive gene IDs (BSAL_30280 / _30285) => TANDEM pair, same arrangement as trypanosomatid AdoMetDC+prozyme (Tb927.6.4460/4470).
- Length pattern matches trypanosomatids: active longer (~370-381), prozyme shorter (~325-331).

## Catalytic residues (MAFFT, human P17707 numbering) — columns read programmatically, not by eye
pos:            Glu8 Glu11 Ser68(PYR) Cys82(donor) Ser229 His243(base) Glu247
human:           E    E     S          C            S      H            E
Tbrucei active:  E    E     S          C            S      H            E
Tcruzi active:   E    E     S          C            S      H            E
Bodo ACTIVE JHX8:E    E     S          C            S      H            E   <- all intact (working enzyme)
Tb PROZYME:      S    E     S          S            E      H            -   <- dead (lost E8, Cys82, Ser229, Glu247)
Bodo PUTATIVE JRG8: S E     C          S            R      -            -   <- DEAD (lost E8, Ser68, Cys82, Ser229, His243, Glu247)

## Key points
- Bodo dead copy shares Glu8->Ser AND Cys82->Ser with the Tb prozyme = shared derived losses -> supports common/ancient origin, not independent convergence.
- Bodo dead copy is even more degraded (also lost pyruvoyl Ser68 -> cannot form cofactor at all).
- => prozyme-like pseudoenzyme pair EXISTS in a free-living kinetoplastid (Bodo), which branches basal to trypanosomatids.
- This DIRECTLY extends the prozyme beyond trypanosomatids (Volkov 2016 "only in trypanosomatids").

## Still needed to make it bulletproof (honest)
1. Phylogeny with these exact seqs: is Bodo dead copy ORTHOLOGOUS to the trypanosomatid prozyme (one ancient duplication) or an independent duplication? (Earlier exploratory tree put Bodo at base of trypanosomatid prozymes - rerun cleanly.)
2. Repeat in Neobodo + other bodonids/Paratrypanosoma to map how far back the pair goes.
3. Structural (FOLDTRACE) confirmation of the dead engine (His243 was a gap in alignment - confirm structurally).
4. Functional caveat: cannot claim the dead copy ACTIVATES the live one (true prozyme function) without experiments; call it "prozyme-like pseudoenzyme" structurally.

## KINETOPLASTID SWEEP (UniProt, 3 Nov) — active + short dead-paralog tandem pairs
Format: ACTIVE (proenzyme, ~360-391 aa) + DEAD (shorter, ~314-331 aa, often "proenzyme-like")
- Bodo saltans (FREE-LIVING): BSAL_30285 (381) + BSAL_30280 (331)  [A0A0S4JHX8 / A0A0S4JRG8]  *** key novel case ***
- Novymonas esmeraldas (tryp): NESM_000475000 (374) + NESM_000475100 (319, "proenzyme-like")
- Leptomonas pyrrhocoris (tryp): ABB37_03063 (391) + 03064 (316, "proenzyme-like")
- Leptomonas seymouri (tryp): ABL78_3891 (390) + 3890 (316, "proenzyme-like")
- Porcisia hertigi (tryp): JKF63_03102 (363) + 03103 (317)
- Angomonas deanei: S9VCY8 (314, "putative") [short copy present]
- Strigomonas culicis: S9WBP4 (382) + S9UUR4 (242)
=> ALL trypanosomatids have the pair (known). Among NON-trypanosomatid/free-living kinetoplastids in UniProt, ONLY Bodo saltans is present — and it HAS the pair. UniProt notably labels the paralogs "proenzyme-LIKE" (= prozyme).
- Paratrypanosoma confusum, Neobodo, Trypanoplasma, Perkinsela, Cryptobia, Azumiobodo, Dimastigella, Rhynchomonas: NO AdoMetDC in UniProt OR NCBI protein (not full proteomes there). Paratrypanosoma (199 NCBI proteins only) => need TriTrypDB / genome assembly to check. THESE ARE THE KEY NEXT TARGETS (Paratrypanosoma = basal trypanosomatid bridge; Neobodo = free-living).

## OTHER THREADS TO PURSUE (per Tuy's "keep finding unique things")
1. Structures: NO experimental structure of any free-living-kinetoplastid AdoMetDC/prozyme exists (RCSB EC 4.1.1.50 = human, potato, T. maritima, Aquifex, + T. brucei AdoMetDC/prozyme complex 5TVM/5TVO/5TVF). => Bodo pair is structurally unexplored. Model both (AlphaFold DB has A0A0S4JHX8/JRG8) + FOLDTRACE-compare the dead copy's fold to the Tb prozyme complex.
2. Palpitomonas (non-kinetoplastid protist): earlier exploratory tree showed it lost the engine INDEPENDENTLY (Ser at His243 position, not the Arg the prozymes share) => candidate CONVERGENT pseudoenzymization. Verify.
3. Candidate labs: Phillips (UT Southwestern) = prozyme discovery + wet-lab assays (THE validation lab); Lukeš (Czech Acad Sci, Inst Parasitology) = Bodo/kinetoplastid genomics; Ealick (Cornell) = AdoMetDC structures. 

## PHYLOGENY (3 Nov) — MAFFT + FastTree(LG), 17 seqs (active + dead across kinetoplastids + human outgroup)
- All 8 DEAD copies form EXACTLY ONE clade (monophyletic; verified: MRCA-of-dead contains only the 8 dead, zero active taxa pulled in). => single ancient gene duplication, inherited; NOT independent per-lineage losses.
- Within the dead clade, BODO_DEAD is the most basal branch (sister to all trypanosomatid prozymes incl. Tb prozyme; node support 0.999). => duplication predates the trypanosomatids (present in the free-living ancestor).
- Internal arrangement of trypanosomatid prozymes: (Lseymouri,Lpyrrhocoris)0.999 +Novymonas, +Porcisia (0.995); (Angomonas,Strigomonas)0.991 — minor nodes weak (0.45-0.64), fine.
- CAVEAT: human outgroup sits on a very long branch and is pulled next to the fast-evolving dead clade = long-branch attraction (LBA). Affects ROOT placement only, not the two core conclusions. For the paper: add more outgroups/bodonids, use IQ-TREE2 with model selection + UFBoot, consider SR4/LBA-aware methods.
- Files (chat-only, NOT in project folder): /mnt/user-data/outputs/my_suggestions_NOT_in_project_folder/prozyme_phylogeny.png + .nwk

## PARATRYPANOSOMA CONFUSUM — BRIDGE TAXON CONFIRMED (3 Nov, Tuy's TriTrypDB BLAST)
BLASTP of Bodo active AdoMetDC vs P. confusum CUL13 proteome (TriTrypDB release-CURRENT) -> TWO hits:
- PCON_0025430 (ACTIVE): "s-adenosylmethionine decarboxylase", 390 aa, 69% id to Bodo active, E=0.0.
  location Pcon_scaffold00006_3:129639-130811(-), is_pseudo=false
- PCON_0025440 (DEAD):  "Adenosylmethionine decarboxylase domain-containing protein", 339 aa, 27% id, E=2e-34.
  location Pcon_scaffold00006_3:132535-133554(-), is_pseudo=false
=> SAME scaffold, SAME strand, ~1.7 kb apart, consecutive gene IDs => DIRECT SYNTENY (tandem pair, physically confirmed - stronger than Bodo gene-ID inference).
=> BOTH is_pseudo=false = maintained protein-coding genes; the dead one is conserved, NOT a decaying pseudogene -> consistent with a functional (prozyme) role under selection.

### Catalytic residues (MAFFT, human numbering), rigorous:
ACTIVE (Human, Tbrucei, Bodo, Pcon_0025430): ALL have E8,E11,S68,C82,S229,H243,E247 (full engine).
Pcon_0025440 DEAD: E8->S, C82->A (lost catalytic Cys), H243->R (lost general base). 
Tbrucei prozyme DEAD: E8->S, C82->S, H243->R.
Bodo DEAD: E8->S, S68->C, C82->S, H243->R.
=> SHARED derived losses across all three independent species' dead copies: Glu8->Ser AND His243->Arg. Same synapomorphy = single ancient origin (matches the monophyletic dead clade in the tree).

## LADDER COMPLETE
free-living Bodo saltans (pair) -> basal Paratrypanosoma confusum (pair, syntenic) -> parasitic trypanosomatids (pair, known).
=> The AdoMetDC prozyme is ANCESTRAL to the whole kinetoplastid lineage and PREDATES parasitism. Literature (Volkov 2016): "only in trypanosomatids" -> now extended to free-living + basal taxa.
Residual honest limits: (1) still computational (sequence+synteny), functional activation untested; (2) tree root needs LBA-robust methods; (3) could add Neobodo/Trypanoplasma for even denser sampling.

## DEFINITIVE RE-RUN (3 Nov) — resolves audit items #2,#5,#7
- Comprehensive MAFFT align, 25 seqs (16 active incl. 6 eukaryotic outgroups + 9 dead), IQ-TREE2 LG+G4, 1000 UFBoot + 1000 SH-aLRT.
- RESIDUE TABLE (now stable): ALL 16 active copies = E8,E11,S68,C82,S229,H243,E247 intact. ALL 9 dead copies lost catalytic Cys82 (->A/D/S, never C) AND His243 (8/9 ->Arg, Angomonas->Lys). Glu8 lost in all dead. (Pyruvoyl S68 lost in some dead, retained in others, but engine dead via C82+H243 loss regardless.)
- TREE: dead clade MONOPHYLETIC at 100/100 (max support). Branching order Bodo_DEAD (basal) -> Pcon_DEAD -> parasite prozymes = matches species phylogeny (free-living -> basal -> parasite). LBA fixed (no active/outgroup inside dead clade).
=> #2 RESOLVED: dead signature robust & universal (C82 + H243->R). #5 RESOLVED: single origin, 100/100, Bodo basal, Pcon next. #7 RESOLVED: verified active+dead pairs in Bodo, Paratrypanosoma, T.brucei, Novymonas, L.pyrrhocoris, L.seymouri, Porcisia; dead copies also in Angomonas, Strigomonas.
- #6: use "catalytically dead paralog" (verified); "prozyme/activator" = function untested (wet-lab).
STILL OPEN: #1 Bodo genomic coordinates (tandem) - need TriTrypDB/NCBI; #3 full Pcon sequences - TriTrypDB; #4 novelty vs genome papers (Jackson2016 Bodo, Flegontov2013 Paratryp); #8 DHS ortholog tree.

## #8 RESOLVED — DHS (second pseudoenzyme system), IQ-TREE LG-type, 1000 UFBoot
- UniProt annotates the trypanosomatid dead DHS directly: Q4GZD1 "Inactive deoxyhypusine synthase" (DHSp, T.brucei), A4ZZ93 "Inactive...DHSL20" (Leishmania). Active = Q38BX0 (DHSc), A4I9K3 (DHS34).
- Bodo: active=A0A0S4KFM4 (BSAL_88720, catalytic Lys329 intact); dead=A0A0S4JCE9 (BSAL_12440, Lys329->Tyr).
- TREE: dead-DHS clade (Bodo_12440 + Tb_DHSp + Leish_DHSL20) MONOPHYLETIC 100/100, Bodo basal (92.9/92 for the Tb+Leish subclade). Active-DHS clade 100/100 separate.
=> Bodo dead DHS is ORTHOLOGOUS to the trypanosomatid inactive DHS = SAME ancient origin, predates trypanosomatids. The DHS dead-paralog system is ALSO ancient & present in free-living Bodo.
=> TWO independent polyamine pseudoenzyme systems (AdoMetDC prozyme + DHS dead paralog) both predate parasitism and are present in free-living Bodo = an ancestral REGULATORY STRATEGY, not a parasite innovation. (Function still untested for both = wet-lab.)

## GAP 2 RESULT (conserved surface of the dead paralogs) — 4 Nov
- Across the 9 dead paralogs, 66 non-catalytic positions are conserved >=85%; at many, the dead copies are MORE conserved than the active enzymes (e.g., human-pos 93 P 100% dead vs 25% active; 134 Y 100% vs 44%; 135 A 100% vs 25%; 138 K 100% vs 38%; also tight blocks ~109-178 and ~255-327).
- Interpretation: the dead paralogs are under strong purifying selection at a specific residue set (not the lost active site) => maintained for a function (candidate: the partner-binding/activation surface). Honest limit: this is conservation evidence, not a structural interface map; needs the complex structure + mutagenesis to call them interface residues.
