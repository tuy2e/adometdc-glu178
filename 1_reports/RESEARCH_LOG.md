# 📋 Research log — AdoMetDC / FOLDTRACE work
**Living master index. Last updated: 2026-11-10.** One place to see every thread, what's done, what's left, and what's a dead end so we don't loop.
Legend: `[x]` done · `[ ]` to do · ⛔ dead end / don't re-chase · ⏸ waiting on you · 💡 my suggestion (not Marc's protocol).

---

## ⚡ Needs YOUR decision (these are blocking things)
- ⏸ **Send Marc the Aim 1 re-run?** The null-result zip + draft email are ready; nothing sent. (Thread 0)
- ⏸ **"Apply B"** → fix the research-note citation errors in `4_research_notes/`. Documented, not applied. (Thread 1)
- ⏸ **"Commit reports"** (or you drop them in) → move corrected reports into `1_reports/`. Not committed. (Thread 1)
- ⏸ **Pick the lead project:** diatom "stuck-on" (#1 bet) vs. keep driving Bodo prozyme. (Threads 4 & 5)

## 🔒 Standing ground rules (so we never slip)
- Read-only on the GitHub repo / `aim1_settings.json` / `sites/`. No pull/push/edit unless you give a Marc command.
- Nothing goes to Marc without you. Suggestions live in `outputs/`, never Marc's project folder.
- Verify every fact against NCBI/UniProt/PMC/RCSB — never from memory. Don't oversell novelty.
- ⛔ Cloud workspace can't reach RCSB/NCBI/UniProt/TriTrypDB directly → use the desktop browser. (Stop trying cloud curl.)
- ⛔ No GPU here → AlphaFold runs on the AlphaFold Server, not in this workspace.
- ⛔ Full Paratrypanosoma sequences are NOT in UniProt/NCBI → only TriTrypDB (needs your login). Stop searching the other two.

---

## Thread 0 — Marc's Aim 1 re-run (1MSV / 5 sites) 📋 Marc's protocol
**Goal:** re-run Aim 1 pipeline on new reference 1MSV with 5 sites (Ser68, Ser229, His243, Cys82, Glu247).
- [x] Pulled Marc's update, renamed run folder, re-ran all six steps
- [x] Result: **NULL** — assembly loss 2.3% vs chemistry 1.9%, overlapping CIs; hypothesis not supported; Glu247 "altered" nuance noted
- [x] Produced `aim1_stuy_adometdc_20261028.zip`; preserved old 1I7B zip + saved old run as `runs/stuy__adometdc_1I7B`
- [x] Drafted email to Marc
- ⏸ [ ] **Send to Marc** — awaiting your go
**Status:** complete except sending.

## Thread 1 — Fact-check Marc's AdoMetDC_project files 📋
**Goal:** check every project file against primary sources; fix errors; never from memory.
- [x] Fact-checked the 4 reports; applied verified fixes (e.g. "Table 1"→"Table 4"; restored K80A numbers ~10% activity, Ka ~400 vs ~20 µM)
- [x] Corrected reports saved to `outputs/corrected_reports/` (6 files)
- [x] Documented research-note citation errors (the "B list": PMIDs/PMCs/DOIs for Graham 2002, Giles & Graham 2008, Ekstrom 2001, M. jannaschii ArgDC, Barker 2009 strain, prozyme/DHS fold numbers, PMC2846215)
- [x] Fixed my OWN α/β pyruvoyl-chain error (it was only in my scratch notes — confirmed NOT in the reports and NOT on GitHub)
- ⏸ [ ] **Apply B** (fix research-note citations) — awaiting your go
- ⏸ [ ] **Commit corrected reports** into `1_reports/` — awaiting your go
- [ ] Refresh `7_run_prep` files (`2_PREDICTION_worksheet.md`, `5_LAB_LOG.md`, `4_note_to_Marc.md`) — still describe the OLD 1I7B/2-site plan; flagged, not edited
**Status:** fixes staged; waiting on your approvals.

## Thread 2 — "Everything about AdoMetDC" deep sweep ✅ 💡
**Goal:** full picture — every structure, what's exceptional, residues, metals, animal use, what's new.
- [x] All **40 PDB structures** cataloged (organism/method/res/ligands) incl. newest 2025 entries
- [x] Metals answered: only *A. aeolicus* 2III has metal; eukaryotic enzymes metal-free; no E. coli structure exists
- [x] "What's new": 2025 multi-temperature crystallography (dynamic loops); 2026 AMD1 ribosome-arrest-peptide paper; AMD1 as prostate-cancer vulnerability
- [x] Animal/human use & regulation (AMD1, ESC self-renewal, translational control)
- **File:** `AdoMetDC_everything_we_know.md`
**Status:** complete (reference doc).

## Thread 3 — Dead-paralog literature review (function + why conserved) ✅ 💡
**Goal:** what the dead copy DOES and why it's been conserved.
- [x] 6 open-access prozyme-function papers synthesized (Willert 2007 → Patel 2018)
- [x] Parallel systems: DHS dead paralog (Nguyen 2013), PRMT1 prozyme (Kafková 2017 / Hashimoto 2020)
- [x] General pseudoenzyme conservation principles (Ribeiro 2019, Murphy 2017, Mace 2021)
- [x] ILL list produced (paywalled papers)
- **File:** `dead_paralog_literature_review.md`
**Status:** complete.

## Thread 4 — Bodo / Paratrypanosoma prozyme finding + interface mapping 🔬 💡 (ACTIVE — co-authorship track)
**Goal:** is the prozyme activation mechanism present/ancestral in free-living Bodo?
- [x] Established active+dead AdoMetDC pairs in *Bodo* & *Paratrypanosoma*; dead clades monophyletic 100/100; DHS parallel dead copy too
- [x] **Interface mapping (AdoMetDC):** Bodo/Pcon active enzymes keep the activation architecture (P31 switch, h1 L/M residues, D169/H172); dead copies keep the validated interface residues (M148'/Y152' → Bodo M150/Y154)
- [x] **Interface mapping (DHS):** Bodo active keeps catalytic K329; Bodo dead lost it (→Y) but keeps binding residues — same signature as validated T. brucei DHSp
- [x] AlphaFold-Multimer input package prepared (FASTAs + analysis checklist with exact Bodo residues)
- **Files:** `activation_interface_mapping_Bodo.md`, `AlphaFold_modeling_package.md`, `alphafold_inputs/*.fasta`
- [ ] **Run AlphaFold-Multimer** on AF Server (Bodo AdoMetDC pair; then DHS pair) → check M150/Y154 at interface
- [ ] Get full **Paratrypanosoma** sequences from TriTrypDB (PCON_0025430 / PCON_0025440) → upgrade "partial" calls
- [x] **Jackson 2016 (Bodo genome, Curr Biol) — CORRECTION: it is OPEN ACCESS (PMC4728078), not paywalled.** Read the full main text (5 Nov): no mention of prozyme, AdoMetDC, polyamines, DHS, or catalytically dead paralogs → does not pre-empt us. Residual: its supplementary gene-cluster spreadsheet (Data S1) not yet checked.
- [x] **Flegontov 2013 (Paratrypanosoma, Curr Biol)**: you uploaded the PDF; read the main text on 5 Nov. It covers morphology plus a phylogenomic placement (114/129/52-protein datasets). There is **no mention of AdoMetDC, prozyme, DHS or polyamines** in the main text, so it doesn't pre-empt us. Its supplementary Table S1 (list of proteins used) is not yet checked.
- [ ] Lab test (for Marc to place): express Bodo enzyme ± dead copy, measure activation
- **Novelty status (honest):** general "prozyme predates parasitism in Bodo" is PARTLY pre-empted (PRMT prozyme already shown in Bodo). The AdoMetDC+DHS+dated-phylogeny package appears unreported. Both genome papers' main texts (Jackson 2016, Flegontov 2013) are now checked and clean; only their supplements remain unchecked.
**Status:** strong and active; best as co-authorship on Marc's FOLDTRACE paper. NOT a dead end.

## Thread 5 — "WOW findings" publication assessment + the diatom lead 🔬 💡 (NEW lead)
**Goal:** find an independent first-author-type project among the FOLDTRACE WOW findings.
- [x] Ranked all WOW findings; did first-pass novelty checks
- [x] **#1 pick = diatom + brown-algae "stuck-on" AdoMetDC** (Glu178→Asn convergence across Ochrophyta; a 2nd independent route to putrescine-independence). Crisp test (human E178Q precedent), tractable model (*Phaeodactylum* B7FPB8), silica-polyamine story
- [x] Confirmed WOW diatom proteins are the **canonical standalone SAMDC** (not Michael 2011's fusion machines)
- **File:** `WOW_findings_publication_assessment.md`
- [x] **Step 1 (5 Nov) — re-analysed your raw FOLDTRACE run data (read-only copy).** Result is BIGGER than the WOW summary: **Glu178 (putrescine sensor) is replaced in 31/31 stramenopile AdoMetDCs** (21 Ochrophyta across 11 classes → Asn; 4 Bigyra → Asn; 6 Oomycota → His/Gln) vs kept in 509/527 other eukaryotes. Gly111 in 21/21 Ochrophyta. Independent losses also in Euglenozoa, cryptophytes, some cestodes, 5 bacteria.
- [x] **Independent check:** *Phaeodactylum* B7FPB8 sequence from UniProt matches all 5 FOLDTRACE-called residues; sequence-only MAFFT alignment to human gives the SAME mapping at all 7 pocket positions (E178↔N198).
- [x] **Novelty check (first pass):** no paper reports a stramenopile-wide loss of the putrescine sensor. Lin 2024 (*Biomolecules*, OA) overexpressed *Phaeodactylum* PtSAMDC but never studied its regulation/pocket (searched full text). Precedents to use as validation: *Tetrahymena* putrescine-insensitive (Pösö 1975); *Acanthamoeba* putrescine-activated (Gupta 1987) and keeps Glu178 in our data ✓.
- **File:** `stramenopile_Glu178_step1_results.md`
- [x] **Step 2 (5 Nov):** pulled all PF01536 AdoMetDCs for stramenopiles + SAR + outgroups from UniProt, plus 38 from the FOLDTRACE run → **248 proteins, 91 stramenopile** (was 31)
- [x] **Step 3:** pairwise MAFFT vs human at 16 positions. **87/87 active stramenopile enzymes lack Glu178** (Ochrophyta N 42, Bigyra N 4, oomycetes H 32 / Q 9). Agrees with FOLDTRACE 125/129 at 178 (31/31 stramenopiles), 830/874 over all pocket positions. Alveolate dinoflagellates/*Perkinsus* have **Asp** (charge kept). *Vitrella* Ala & *Raphidocelis* Arg = likely alignment errors (multiple alignment says Glu).
- [x] **Step 4:** IQ-TREE (LG+R6, 171 seqs, UFBoot + SH-aLRT, ancestral reconstruction). **Within stramenopiles: one early Glu→Asn, then Asn→His in oomycetes, His→Gln in a *Phytophthora*/downy-mildew subgroup; Glu never regained.** The gene tree groups stramenopiles with kinetoplastids + *Eutreptiella* + *Guillardia* + dinoflagellates/*Perkinsus* (UFB 95), but support drops to 63 with the 35 pocket columns removed. AU test: each species-tree grouping alone is not rejected (p 0.58 / 0.18 / 0.066); all three together are rejected (p 0.004). → **Shared vs independent loss across supergroups is UNRESOLVED**; don't claim either.
- [x] **Step 5:** validation vs published putrescine responses. 6/6 sequenced activated enzymes outside trypanosomes keep Glu (human, yeast, *Neurospora*, *C. elegans*, *Onchocerca*, *A. castellanii*). E178Q (human; Bale 2008 Table 4 = 4× basal, no stimulation) and E194Q (*C. elegans*; 100× basal) are the precedents. ⚠️ **Exception: *T. cruzi*** (Ser178, still activated via Asp174), so the prediction is "likely", not certain. *Tetrahymena* has no class-2 enzyme → NOT a validation point (step-1 note corrected).
- [x] **Step 6:** structure figure (human 1MSV vs diatom/oomycete AF models, core RMSD 0.76–0.95 Å). Diatom pocket intact + acidic, Asn at 178, empty site (E178Q-like). Oomycete pocket partly self-filled (His/Asn/Tyr), so oomycetes may differ.
- [x] Re-verified the "E178Q ~4× more active without putrescine" number directly in Bale 2008 Table 4 ✓ (no correction needed)
- **Files:** `stramenopile_Glu178_steps2-6_results.md`, `Fig_tree_Glu178.png`, `Fig_pocket_human_vs_stramenopile.png`
- [x] **Novelty lock: Michael 2011, read in full (5 Nov; you found it open-access, so my "paywalled" label was wrong).** It is about the bacterial-type (class-1b) AdoMetDC–aminopropyltransferase fusions. On the canonical enzyme it only says orthologues are present, with uORFs suggesting translational feedback. Nothing on the pocket, Glu178 or the enzyme's evolution → **not pre-empted.**
- [x] **Wider sampling (5 Nov evening):** NCBI + UniParc, which recovers the TSA/MMETSP entries UniProt dropped → **339 active enzymes**. Contamination screen set in advance: 10 removed (e.g. "*Symbiodinium*" = 100% *Cafeteria*; the dinotoms *Durinskia*/*Kryptoperidinium* carrying diatom endosymbiont copies; *Tiarina* with diatom prey). Results:
  - **Stramenopiles 97/97 lack Glu.**
  - **Cryptophytes 9/9 lack Glu** (7 species incl. *Goniomonas*).
  - Kinetoplastids 7/7 and the euglenid lack Glu; **the diplonemid *Diplonema* KEEPS Glu**.
  - Haptophytes 32/32 keep Glu, as do Rhizaria and heteroloboseans.
  - Dinoflagellates have **two copies**: an Asp copy and a Gln copy (the Gln one is often fused to ODC).
- [x] **Expanded tree (187 seqs, Q.pfam+R7, UFBoot/SH-aLRT, ancestral reconstruction):**
  - **The Glu-less AdoMetDCs of stramenopiles + cryptophytes + kinetoplastids/euglenid + the dinoflagellate/*Perkinsus* Asp copy form ONE gene lineage** (UFB 88), with a **single Glu→Asn change at its base** (P 0.99). Glu is never regained.
  - The robustness checks agree: without pocket columns, the same lineage appears (87/64); without kinetoplastids, it appears again (UFB 85).
  - **AU tests reject moving each group to its species-tree position:** cryptophytes p 0.025, kinetoplastids/euglenid p 0.001, dino Asp copy p 0.003, all three p 0.0001. Still rejected without pocket columns (p 0.034 / 0.012 / 0.002 / 0.0006).
  - → **Changes the step-4 answer: "one loss in a cross-supergroup gene lineage"** (horizontal transfer or an ancient duplicate with losses; one gene can't tell which). Stated with limits in the doc.
- [x] ⚠️ **Robustness round 2 (7 Nov): one conclusion revised.** Ran (1) a Li 2019-style alternative tree and (2) a C60 mixture model (PMSF). Record: `/root/audit/strameno/expand2/RESULTS_round2.md`.
  - **Holds under both models:** two Glu-less groups each break the species tree. Group (a) is stramenopiles + the dino/*Perkinsus* Asp copy; group (b) is kinetoplastids/euglenid + cryptophytes.
    - Standard-model AU: cryptophytes p 0.033, kinetoplastids p 0.0005, Asp copy p 0.004, all p 0.00001.
    - C60 AU: 0.008 / 0.0002 / 0.013 / 0.0005.
    - Without pocket columns the cryptophyte test is borderline (p 0.068).
  - **Not settled:** whether (a) and (b) are ONE lineage. The best tree joins them (UFBoot 88 standard, **90 under C60**), but keeping them apart is not rejected: Li-style p 0.18 standard / 0.37 C60 / 0.48 without pocket; next to the alveolate copy p 0.18 / 0.39.
  - → **"One loss" is downgraded to "1 or 2 losses".** Group (a) could be an ancient SAR duplicate, which would need no gene transfer. Group (b), cryptophyte + Discoba, is the clear cross-supergroup case, and Li 2019's tree already shows that pairing.
  - Trp109 sits in both groups (ST 36/45, Asp copy 16/16, KE 6/7, CR 5/7) and almost nowhere else. One origin, or the same pair of pocket changes made twice (convergence). Either reading is interesting.
  - Report updated (header, short version, step 4, limits, paper section). Team doc updated.
- [x] **Round 3 (9 Nov): ancestral reconstruction + spermidine synthase (SpdSyn) tree.** Records: `/root/audit/strameno/expand2/RESULTS_round2.md` and `/root/audit/spdsyn/RESULTS_spdsyn.md`.
  - **Ancestral reconstruction:** the two-group trees need **2** Glu178 losses (one per group); the best tree needs 1. Glu is never regained in any tree. Trp109 is 1 gain on the best tree but 1–3 (unstable) on the others.
  - **SpdSyn tree:** 142 seqs, LG+R6.
    - Stramenopiles + dino/*Perkinsus* forced together: rejected (p 0.002). All four groups together: rejected (p 0.0006).
    - Cryptophytes + kinetoplastids: **not** rejected (p 0.51, Δ0.25).
    - → The AdoMetDC pattern is gene-specific; the partner gene did not travel with it. The cryptophyte–kinetoplastid pairing is weakly echoed.
    - **It does not settle one-vs-two.**
  - ⚠️ NCBI BLAST never ran (queued for more than 2.5 h), so the SpdSyn tree has only 2 cryptophytes, 0 euglenids and 3 dinos. Retry later if needed. *(10 Nov: euglenid gap filled from search H; see "Euglena added to both trees" below.)*
- [x] **New shared pocket change: Ser109→Trp** in 84% of that lineage vs 1% elsewhere. In the diatom/brown-alga models the Trp sits about 4 Å from the putrescine site, roughly where human Phe111 sits (Phe111 is Gly in Ochrophyta).
- **Files:** `stramenopile_Glu178_steps2-6_results.md` (updated), `Fig_tree_Glu178_expanded.png` (new main tree)
- **Plan to December (agreed 5 Nov): keep checking until we're sure it's unexplored. Work in chunks; save after each.**
  - [x] **1A: gene-transfer and family-tree literature search.** 327 papers screened. No report of Glu178 loss, no class-2 AdoMetDC tree inside eukaryotes, no class-2 gene transfer among eukaryotes, no Ser109→Trp.
    - Closest neighbour: **Li et al. 2026 *JBC*** (Michael lab), which covers the partner enzymes across the same lineages. Must cite.
    - ✅ **Fan et al. 2020 *Sci Adv* overlap resolved: none.** Their diatom "AdoMetDC" gene-transfer hits are ~80-aa bacterial-type PF02675 fragments (*Thalassiosira* EED92531 = B8C3B1), not our class-2 enzyme.
    - ⚠️ WOW #3 (class-2 AdoMetDC in bacteria) is **pre-empted** by Li et al. 2019 *Biochem J*.
    - File: `novelty_reading_list.md` (52 papers in tiers 1–3 + 117 citing papers)
  - [x] **1B: closest papers read in full (or their full text keyword-searched).** All clean: Li 2026, Li 2023 *JBC*, Li 2023 *PNAS*, Michael 2016 *JBC* and 8 genome/gene-transfer papers (zero AdoMetDC mentions).
    - Couldn't open at first: **Li 2019 *Biochem J*** (most important; has a class-2 tree), Michael 2017 and 2016 *Biochem J*, Pegg 2009, Li 2015, Fuell 2010. *(6 Nov: you got 4 of them, now read; Li 2019 and Pegg 2009 are on ILL.)*
    - Indirect check: bacterial class-2 AdoMetDCs mostly keep Glu178 (42/54), so they're not from our lineage.
  - [x] **1C: residue and biochemistry search.**
    - Nothing on Ser109→Trp or Asn178 anywhere. Glu178 has only been studied in human mutants; Singh 2013 notes *Leishmania* lacks 1 of 3 putrescine residues.
    - Standard classification (Bale & Ealick 2010: 2a / 2b-I / 2b-II) has **no class for stramenopiles or cryptophytes**.
    - Lys80 (the 2b-I marker) is kept in stramenopiles; cryptophytes have Leu80 and kinetoplastids/euglenid Ile80.
    - *Plasmodium*'s lack of response is explained by built-in Arg/Lys (Wells 2006).
    - **No stramenopile, oomycete or cryptophyte AdoMetDC has ever been measured ± putrescine.**
    - **Library needed:** Aleksijevic 1979 *BBA* (*Euglena*; doi 10.1016/0005-2787(79)90096-0) and Pösö 1976 *BBRC* (doi 10.1016/0006-291x(76)90033-4).
  - 📥 **Paper access status (6 Nov):**
    - **You have the PDFs** (open access): Michael 2017 *Biochem J*, Michael 2016 *Biochem J*, Li 2015 *Mol Microbiol*, Fuell 2010 *Plant Physiol Biochem*. ✅ All 4 read in full 6 Nov.
    - **ILL requested:** ~~Li 2019 *Biochem J*~~ (✅ you got it 7 Nov, read), Pegg 2009 *Essays Biochem* (DOI also sent to Marc), Aleksijevic 1979 *BBA*, Pösö 1976 *BBRC*.
    - Rule from now on: **read every paper in full**, never just the abstract; ask you for a PDF if I can't open it.
  - [x] **1B-PDF: your 4 uploaded PDFs read cover to cover (6 Nov). All 4 are clean (none pre-empts us).** Notes are in `novelty_reading_list.md` under "Your uploaded PDFs".
    - **Michael 2017 *Biochem J*:** class 2 is a pre-LECA fusion "found in all eukaryotic lineages". The *Shewanella* class-2 copies are "likely" from a eukaryote, with no tree shown. Ku & Martin 2016 / Artamonova 2015 warn that bacterial-looking genes are often contamination; worth citing for our screen. No class-2 tree, no Glu178, none of our groups.
    - **Michael 2016 *Biochem J*:** same textbook view (class 2 "in all eukaryotic supergroups"). Kinetoplastids lack mRNA (upstream-ORF) control of AdoMetDC and use prozyme instead. Putrescine activation is never mentioned. Its *Shewanella* source is **Illingworth 2003**.
    - **Li 2015 *Mol Microbiol*:** ciliates lost class 2 (fits our survey). Its only tree is for **spermidine synthase** (neighbour-joining). In it, the cryptophyte *Guillardia* sits next to kinetoplastids and stramenopiles sit with red algae; no support values are shown for those nodes.
    - **Fuell 2010 *PPB*:** green algae and plants only. Carnation/potato AdoMetDC is not putrescine-stimulated (Bennett 2002, already in our step 5). No residues, no tree.
    - **New to-read (moved to Tier 1):** Illingworth 2003 *FEBS Lett* (🟢 free; may hold a class-2 tree with eukaryotes) and Green 2011 *Mol Microbiol* (🟢 free; *Micromonas* AdoMetDC–ODC fusion).
    - **Green 2011 *Mol Microbiol*: main text read in full (6 Nov; PMC author manuscript). Clean.** It covers bacterial class-1b fusions. The only class-2 content is the *Micromonas* AdoMetDC–ODC fusion (XP_002502359 = C1E5E7), independent of *Plasmodium*'s. My pairwise check: *Micromonas* keeps **Glu178** and has plant-like Arg (at L13) and Lys (at F111). Our dinoflagellate ODC–AdoMetDC fusion is not mentioned in the main text. ✅ **Supplement S1–S7 read (your PDF):** all bacterial/archaeal, and S7 is a generic pathway sketch that names no organisms. **So Green 2011 doesn't report our dinoflagellate fusion.** Still to do: a short literature search before calling it new, since it may be inherited from the ancestor shared with *Plasmodium*.
    - ✅ **Illingworth 2003 *FEBS Lett* (549:26–30): read in full (your PDF). Clean, and it has no tree.** The *Shewanella* "got it from a eukaryote" claim rests on one 3-sequence alignment (*Arabidopsis*, human, *S. oneidensis* NP_717479). The only known paper with a class-2 tree including eukaryotes is still **Li 2019** (interlibrary loan).
    - ✅ **Li 2019 *Biochem J* 476(18):2579–2594: read in full (your PDF, 7 Nov). It doesn't pre-empt us, but we must cite it and answer its tree.**
      - Its neighbour-joining class-2 tree (Fig. 2) contains our groups. Kinetoplastids + *Guillardia* are drawn as "earlier diverging" and split from all other eukaryotes at 100%; that agrees with our kinetoplastid + cryptophyte group (100).
      - Stramenopiles and *Perkinsus* + *Symbiodinium* (the *Symbiodinium* sequence is our Asp copy) sit in the core of their tree with no support shown. Ours puts them in the same lineage (UFBoot 88).
      - No Glu178, no putrescine site, and no cross-supergroup claim. Their Discussion even says eukaryote-to-eukaryote transfer is rare (Keeling & Palmer 2008, PMID 18591983), so class 2 was present in the last common ancestor of eukaryotes.
      - Their Fig. S2: *Shewanella* and *Bacteriovorax* class-2 **keep Glu178**, so the bacterial copies are not from our lineage (confirms the 42/54 indirect check).
    - 💡 **Possible chunk 2D (my idea, not started):** build a spermidine-synthase tree from the same organisms. It would tell us whether the partner gene travelled with our AdoMetDC lineage.
    - ✅ **ILL batch read in full (9 Nov, your 5 PDFs).** Notes are in `novelty_reading_list.md`. **None pre-empts us.**
      - **Aleksijevic 1979 *BBA* 565:199–207 (PMID 116684): ⚠️ *Euglena gracilis* AdoMetDC IS putrescine-activated, 3–10×** (dialysed extract; K<sub>a</sub> 0.6 µM; specific; no Mg²⁺ effect). Euglenids sit in our Glu-less Discoba group. *Eutreptiella* (our only euglenid) has Asn178, but *E. gracilis*'s own sequence is in no protein database.
        - **Likely a second Discoba exception after *T. cruzi*.** So we now state the prediction for stramenopiles (and cryptophytes), not for the whole Glu-less group. The report's Step 5 is updated.
        - **Needs:** tblastn of AdoMetDC vs *E. gracilis* TSA (taxid 3039), added to your BLAST list.
      - **Pösö 1976 *BBRC* 68:227–234 (PMID 813639):** 7 sources; the only protist is *Tetrahymena* (1.0×, already explained by enzyme class).
        - Added 2 validation rows: rat (P17708, 1.8×) and pike-perch (A0A8C9YFM5, 2.0×), both keep Glu178. That makes **8/8 sequenced**.
        - Their old idea: putrescine activation goes with having a spermine synthase. 💡 *Idea, not run:* overlay spermine-synthase presence (from Li 2026 *JBC*) on our Glu178 tree.
      - **Pegg 2009 *Essays Biochem* 46:25–45 (PMID 20095968):** review; classes IIA/IIB/IIC (≠ Bale & Ealick naming, so pick one). Good citation for the putrescine-site residues (Phe111/Phe285 stacking; Asp174/Glu15/Thr176; Glu178/Glu256/Ser113 via water). No protist data.
      - **Opperdoes 2016 *J Eukaryot Microbiol* 63:657–678 (PMID 27009761), re-read:** no AdoMetDC or prozyme in the text.
        - ⚠️ **Correction to my 3 Nov note:** Fig. 3 *does* label AdoMetDC (*speD*) as present in all 9 kinetoplastid groups (I had read only the text layer). No copy number and no prozyme, so the Bodo novelty conclusion is unchanged.
      - **Li 2019 final version vs accepted manuscript:** same text (sentence-by-sentence check) and same Fig. 2. Cite the final version.
      - Citation checks on PubMed: Mitchell & Rusch 1973 (*Physarum* not activated) = PMID 4735856. Wells *P. falciparum* model = online 2005, print **2006**, 24(4):307–318.
    - ✅ **Search G done (10 Nov, you ran it; RID F5EJKSE8014):** *E. gracilis* AdoMetDC = one gene, 10 transcripts from 2 independent assemblies, 65% identical to *Eutreptiella*. It has **Asn178**, R13, D174, E256 (also Q15). **So *Euglena* is a confirmed second Discoba exception** (no Glu178, still activated). No second AdoMetDC-like transcript, so no prozyme candidate. The prediction stays stramenopile (+ cryptophyte) specific.
    - ✅ **Search H done (10 Nov; RID F5FJXSNR014; tblastn *T. brucei* SpdSyn vs *E. gracilis* TSA):** one clear spermidine synthase (GEFR01020952.1 = GDJR01022181.1, full 285-aa ORF, 46–47% identical to *T. brucei*). Two weaker aminopropyltransferase-family genes are also present (GDJR01022779 ~41%; GEFR01015322 ~33%), unclassified. That fits *Euglena* making norspermidine/norspermine (Aleksijevic 1979), but the tree will decide what they are. The proteins are saved in the session workspace (`/root/audit/euglena/euglena_proteins.fa`) for tonight's SpdSyn and AdoMetDC trees. **This fills the "no euglenids" gap in the SpdSyn tree.**
    - 🔎 **Why *Euglena* matches *T. cruzi* (checked 10 Nov on the existing main tree's ancestral reconstruction, expA, Q.pfam+R7):**
      - The ancestor of the kinetoplastid + euglenid AdoMetDCs (Node86) already had **Arg13, Asp174, Asn178** (posterior 0.999 / 0.997 / 0.983). So the pattern is **inherited from a shared ancestor, not convergent**. *T. cruzi* later changed Asn→Ser at 178.
      - NCBI Taxonomy confirms both are Euglenozoa (Discoba; Euglenida vs Kinetoplastea). *Diplonema* is also Euglenozoa but its AdoMetDC keeps Glu and groups with heteroloboseans.
      - ⚠️ **Same reconstruction, group (b) ancestor (kinetoplastids + euglenid + cryptophytes, Node79):** Arg13, Asp174, Asn178 (all ≥0.99). The cryptophyte ancestor (Node80) had Arg13 too; modern cryptophytes have Lys13. **So cryptophytes look like the exceptions at position 13, not like stramenopiles.** The stramenopile ancestor (Node35) had Ile13 (0.67), Asp174, Asn178.
      - **Consequence:** the cleanest prediction is **stramenopiles only**; cryptophytes move to "uncertain". Whether position 13 matters at all is unknown until we read Clyne 2002 / Beswick 2006 in full. *Plasmodium*'s Arg at a nearby spot is proposed to be a built-in activator (Wells 2006), which is the opposite outcome.
    - ✅ **Euglena added to both trees (10 Nov, new session).** Record: `3_data/test_results/RESULTS_euglena.md`. Figure: `2_figures/Fig_tree_Glu178_expanded_with_Euglena.png`. Working files: `8_full_workspace_backup/workspace_part3_euglena_trees.zip`.
      - **Method note:** MAFFT here (v7.505) doesn't reproduce the 9 Nov alignments, so I didn't realign. Each *Euglena* protein was added to the original untrimmed alignment (`mafft --add --keeplength`, L-INS-i settings), then the original ClipKIT mask was applied. The original rows are byte-identical. Same IQ-TREE 2.0.7 commands and models as before.
      - **Search H's other two hits:** GDJR01022779 is a **thermospermine synthase (ACL5-like)**, grouping with dinoflagellates and haptophytes. GEFR01015322 is in a small divergent PABS group outside the SpdSyn clade. Neither goes in the SpdSyn tree.
      - **AdoMetDC (188 seqs):** *Euglena* + *Eutreptiella* 100/100. With kinetoplastids 99.9/100. Group (b) 99.5/100. 4-group lineage UFB 90. *Diplonema* still with heteroloboseans. AU: H_kin with *Euglena* is still rejected (Δ48.3, p 0.019).
        - ⚠️ The ML tree now joins stramenopiles to group (b), with the Asp copy branching first. The 9 Nov arrangement is tied: Δ0.7, AU 0.52 vs 0.58. On the new arrangement the first state after Glu is uncertain (N 0.41 / E 0.30 / D 0.29). On the old one it is still E→N 0.996. **Glu is never regained on either tree.** One vs two losses stays unresolved.
        - Residues: *E. gracilis* R13, **Q15**, **Y109**, D174, **Q176**, **N178**, E256. It is activated by putrescine yet lacks Glu15, Thr176 and Glu178. Euglenids have Tyr109; the reconstruction puts the W→Y change at the euglenid ancestor.
        - Small correction: **Ser178 is shared by all 6 kinetoplastids** (kinetoplastid ancestor S 0.98), so it isn't specific to *T. cruzi*.
      - **SpdSyn (143 seqs):** *Euglena* SpdSyn + *Diplonema* 99.2/100, then + kinetoplastids (92.2/41). That is the normal Euglenozoa grouping, as in the species tree.
        - AU (7 trees): cryptophytes + kinetoplastids + *Euglena* (AdoMetDC group b exactly) **rejected**, p 0.016 (0.049 in the 9-Nov-style 5-tree set). Cryptophytes + kinetoplastids with *Euglena* free: p 0.48. Cryptophytes + all Euglenozoa: p 0.47. Group (a) p 0.0003; 4-group lineage p 3e-8.
        - → **The two genes disagree about euglenids.** This is more support for the AdoMetDC pattern belonging to that gene alone.
      - Not rerun with *Euglena*: H_crypto, H_dino, H_all, H_li/H_sar, the C60 runs.
    - ✅ **BLAST searches A–H all checked (10–11 Nov).** Index of every RID, its settings and its verdict: `3_data/blast_results/BLAST_INDEX.md`. Reports are saved next to it (Download All → Text).
      - D (defaults) was cut off at 100 hits, so it was replaced by D2 (Rhizaria) and D3 (haptophytes), both at 500 / 1e-20.
      - Every other search was complete.
      - **Katablepharid AdoMetDC not detected** (E, E2, E3; human query at E 0.05 and *Guillardia* query at 1e-5). Only a small *Roombia* transcriptome exists, so say "not detected".
      - Centrohelids have AdoMetDC (E; Glu178 not yet checked).
      - ⚠️ **Gap:** the plan for search E (Project Map) also listed **Telonemia and Picozoa**, but the run covered only katablepharids and centrohelids, so those two are still unsearched. The optional *Physarum* search (not putrescine-activated: does it keep Glu178?) wasn't run either.
      - **Next:** add the new sampling to the trees:
        - SpdSyn: A2 cryptophytes (+ *Roombia*, labelled separately), B euglenids/diplonemids, C dinos (dinotom screen), D2 Rhizaria, D3 haptophytes (filter out other aminopropyltransferases first).
        - AdoMetDC: F *E. longa*, *Rhabdomonas* and diplonemids.
    - ✅ **Round 4 (10 Nov): all the BLAST sampling added to both trees.** Record: `3_data/test_results/RESULTS_round4.md`. Figures: `2_figures/Fig_spdsyn_tree_round4.png` and `Fig_tree_Glu178_round4.png`. Working files: `8_full_workspace_backup/workspace_part4_round4.zip`.
      - **Lost and redone:** the previous chat's sequence pull and its post-2:37 am log/report edits never reached the folder. All sequences were re-read from the saved BLAST reports instead (HSP regions, not full ORFs).
      - **Spermidine synthase (245 seqs; LG+R6 again):** 598 candidates → 314 non-redundant → 163 classified as SpdSyn on the family tree → 21 dropped by the contamination rule (e.g. *Teleaulax* ↔ *Phaeocystis* 90%) → 102 added after 90% reduction.
        - **Kinetoplastids + 12 photosynthetic cryptophytes now form a clade (93.5/90).** The group (b) pairing is now echoed in the partner gene. *Goniomonas* sits apart.
        - Euglenids group with diplonemids (96.7/100), not with kinetoplastids.
        - Group (a) and the 4-group lineage still absent.
        - AU (7 trees):
          - CR+KE (all 14 cryptophytes) p 0.457, and CR + all Euglenozoa p 0.485: both not rejected.
          - CR+KE+euglenids (group b exactly) p 0.026, rejected.
          - Group (a) p 0.001 and the 4-group lineage p 0.00006, both rejected.
          - CR + red algae p 0.044 (borderline; the last three constrained runs used -nstop 30).
      - **AdoMetDC (200):** + *E. longa*, *Rhabdomonas*, a 2nd *Eutreptiella* copy (69% to the first), 7 diplonemids, 2 centrohelids.
        - Euglenids all sit in group (b) (99/100).
        - **All 8 diplonemids keep Glu178**, grouping with heteroloboseans. Euglenozoa-together is rejected (p = 0.017).
        - Centrohelids keep Glu178.
        - ⚠️ *Blastocystis* (long branch) jumps out of the lineage (UFB 48), which drops lineage support to UFB 65. Putting it back is not rejected (p = 0.27), and then there is one E→N (0.994). Glu is never regained.
        - **New:** *E. longa* = *E. gracilis* at the pocket (Q15 shared). **_Rhabdomonas_ has Lys178** (a positive charge, plant-like; a hypothesis only).
    - ✅ **Activation papers read in full (10 Nov, in a separate chat on your PDFs):** Clyne 2002, Beswick 2006, Bennett 2002, Wells 2006, Williams 2011. Notes and data check: `5_literature/paper_notes_activation_mechanisms.md`. (Ndjonka 2003 ×2 not obtained.)
      - *T. cruzi* needs **Asp174 only** to bind putrescine. Rebuilding the human Lys80/Glu178 pair does nothing. **Arg13 is a processing determinant**, not an activation residue.
      - *Plasmodium* keeps Glu178 and the whole relay but has Lys at 174 plus built-in Arg/Lys, and is not activated. So **Glu178 is neither necessary nor sufficient**.
      - **Our alignment:** Ochrophyta keep Asp174 (20/22 acidic), Glu15, Glu256 and Lys80, and lack *T. cruzi*'s Arg13/Gln176/Met113/Ile80. Their pocket is human-like except at 178 (Asn), 111 (Gly) and 109 (Trp).
      - Oomycetes: Arg111 in 5/12 + His178 look plant-like. Cryptophytes have His285 (5/7).
      - → **Prediction downgraded from "likely" to "plausible"**: the closest precedent is human E178Q, but Asp174 binding may persist. Oomycetes get a separate (plant-like) hypothesis. The cryptophyte "uncertain" reason changes (position 13 = processing). Docking (3B/3C) moves up the list.
      - Reports and team doc still say "likely" in places; reword them when writing up.
    - ✅ **Putrescine docking (10 Nov; to-do 3B/3C partly done).** Record: `3_data/test_results/RESULTS_docking.md`.
      - Vina 1.2.5 on human 1MSV (+ E178Q in silico) and 13 FOLDTRACE AlphaFold models (read from the Downloads copy, unchanged).
      - **Controls pass:** human redock 0.7 Å; *Acanthamoeba* and lamprey dock in the site; the plant doesn't (pocket filled, 5 Å³); *T. cruzi* docks 6 Å away (shifted site).
      - **Test enzymes:** diatom, brown alga, cryptophyte and *Perkinsus* copy have 0 in-site poses (best 4–6 Å away, *T. cruzi*-like). 2 of 3 oomycetes have filled pockets (plant-like). *Cafeteria* and *Aphanomyces* dock in the site occasionally (2/9).
      - → The pockets are remodelled, but docking can't tell "not activated" from "activated via a shifted site". **The prediction stays "plausible"**, and the oomycete plant-like hypothesis is supported.
      - The haptophyte (keeps Glu) also has a filled model pocket, a caution about apo models.
      - AD4 charge-aware scoring failed the control, so it isn't used. E178Q docks like wild type (signalling ≠ binding).
      - **Bug caught:** the first model alignment applied the transform twice. It was fixed before any result was used.
    - 📝 **Publication framing (10 Nov, discussion; Marc decides).** A computational-only paper is publishable if the prediction isn't the headline.
      - Lead with the observation (Glu178 lost in all stramenopiles plus the other groups) and the evolution result (the Glu-less groups cluster against the species tree).
      - The activation prediction goes in the discussion as a testable hypothesis ("plausible").
      - A lab measurement (e.g. Lin lab PtSAMDC) would strengthen it but isn't required.
      - A molecular evolution or protist journal fits better than JCIM. State one-vs-two losses as unresolved.
    - 📄 **Model paper: Begg, Mason & Zimmer 2026, *J. Chem. Inf. Model.* 66:10266** (Luke and Maddie; open access). Computational only, built with Claude Code, untested hypotheses stated as such. Things to copy:
      - An "AI-assisted computation" methods paragraph plus a supplementary note on an AI failure. Ours: the docking double transform, the AD4 failed control, the MAFFT version mismatch.
      - Code on GitHub plus a per-entry data table.
      - Explicit handling of over-represented groups (our 90% redundancy reduction).
      - A blunt Limitations section.
      - An AlphaFold-vs-crystal check (→ to-do 3D).
    - ➡️ **Next direction (10 Nov, Tui's choice): pocket-chemistry survey**, computational only. Can structure-based pocket chemistry (electrostatics, pKa, H-bonds, filling, cavity, docking) separate the known activation types, and where do Glu-less lineages fall?
      - Go/no-go on crystals plus enzymes with known responses first, then the 1,000 FOLDTRACE models.
      - Plan and ground rules: `1_reports/HANDOFF_pocket_chemistry_survey.md`.
      - The evolution paper stays the fallback (core claims: findings 1–3).
    - ⚠️ **Open item from this round:** the report said stramenopiles have "no prozyme". **That was never searched.** I've reworded it. Four unexamined low-identity "dead" class-2-like hits (3 *Cafeteria*, 1 *Schizochytrium*) need a quick look before we claim it.
      - ✅ **Looked at, 11 Nov.** Not prozymes. The 3 *Cafeteria* copies are 85–89% identical to each other (recent triplication or assembly redundancy) and read Gly or another non-functional residue at **every** diagnostic column (S68→G, C82→V, D174→L, E178→A, S229→G, H243→G/gap, E256→G) at 20–22% identity — the signature of a bad alignment, not of a dead paralog that keeps the fold. *Schizochytrium* CAE0742880.1 (660 aa, 25%) looks unrelated or fused. **Wording stays "not searched"**: a profile/HMM search built from trypanosomatid prozymes is still the proper test. A Europe PMC search finds no prozyme reported outside kinetoplastids. Details: `1_reports/BACKGROUND_factors_pocket_chemistry_survey.md` §6.
  - ✅ **Pocket-chemistry survey, Stage 1 started (10–11 Nov).** Pre-registration (measurements + gate, fixed before any structure was measured): `1_reports/PREREG_pocket_chemistry_measurements.md`. Test-set design and the four-quadrant rationale: artifact "AdoMetDC Stage 1 Test Set". Scripts: `4_scripts/stage1_fetch_structures.py`, `stage1_prep.py`, `stage1_measure.py`, `stage1_charge_control.py`.
    - ✅ **Literature check (step 1) done**: nobody has computed or compared AdoMetDC pocket chemistry across species. Nearest prior art is **Bale & Ealick 2009** (*Amino Acids* 38:451, open access) — descriptive, no computation, only the species with crystals. Notes: `5_literature/litcheck_pocket_chemistry_survey.md`. **Not yet read in full; should be before writing up.**
    - ✅ **Trypanosome crystals exist after all**: **5TVO** (*T. brucei* monomer, 1.48 Å) and **5TVF** (+ prozyme heterodimer, 2.42 Å), Volkov 2016. Plus extra human crystals with putrescine (1JL0, 1I7B) and an **apo** human crystal (**1JEN**).
    - ✅ **Pipeline validated against the 10 Nov round**: cavity volumes reproduce exactly (human 33.9 Å³, *Selaginella* 4.6 Å³) and docking agrees to within one pose.
    - 🔶 **Early signal (6 structures, not the gate):** potential at the putrescine position separates activated (human −29, *Acanthamoeba* −37, lamprey −39, *T. brucei* full-length −25 kT/e) from the plant (**+13 kT/e**, i.e. the opposite sign). **Electrostatics sees what docking could not**: human E178Q shifts −29 → −18 although it docks identically to wild type.
    - ⚠️ **The truncated *T. cruzi* model was distorting the 10 Nov result.** The full-length *T. brucei* model (Q587A7, site pLDDT 98) gives −25.1 kT/e, cavity 29.5 Å³, 2/9 in site, versus the truncated *T. cruzi* model's −9.7, 70.1 Å³, 0/9. The missing ~20 N-terminal residues are exactly the **autoinhibitory** segment (Volkov 2016, Velez 2013).
    - ✅ **Bulk-charge confound ruled out** (`stage1_charge_control.py`): net charge varies by 18 e across the set but solvent-shell potentials stay within ±1.4 kT/e, and *Selaginella* is the **second most negative protein** in the set while having a **positive** site. The separation is local pocket chemistry, not amino-acid composition.
    - 🆕 **New result — the signalling relay is intact everywhere; the lid is not.** Bale 2008's relay (**Glu11, Lys80, Ser229, His243**) is invariant across all 19 lineages, **stramenopiles included** (22/22 at K80/S229/H243 in Ochrophyta). So stramenopiles kept the *wire* and lost the *sensor*. Two further points: **Lys80 is lost only in cryptophytes (L 7/7) and kinetoplastids (I 6/6)** — both members of group (b), a shared derived character the project had not noticed, worth testing against the ancestral reconstruction; and **ochrophytes have replaced the whole aromatic lid** (human F285/F315/Y318/F320 → M/V/C/C), a second independent change to the sensing apparatus. Lid calls are consensus-level on a trimmed alignment and need per-sequence checking.
    - 📄 **Amendment 1 to the pre-registration (11 Nov, before any measurement):** 5TVO/5TVF removed from the G2 noise estimate (they are autoinhibited vs activated **states**, not replicates); **1JEN added** so that lid closure is separated from crystal-vs-model. The six measurements and the gate are unchanged; the relay/lid positions are deliberately **not** added to the gate.
    - 💡 **Functional hypothesis for the loss (for the discussion, labelled as a hypothesis):** diatoms run a high, constitutive demand for decarboxylated AdoMet to build long-chain polyamines for biosilica (Kröger 1999/2000/2001; Michael 2011, already read; Knott 2007 on thermospermine). An enzyme serving constitutive bulk biosynthesis has little use for feedback activation, so losing the sensor may be adaptive rather than merely tolerated — the same endpoint plants reach by filling the pocket. **Also raises a question to settle from Michael 2011: is the eukaryotic class-2 copy the main flux enzyme in diatoms at all?**
    - 🧪 **Wet-lab partner lead, now concrete: Lin HY et al. 2024, *Biomolecules* 14:372** (open access). The Lin lab cloned **PtSAMDC** with PtSDS1/PtSDS2 from *Phaeodactylum*, localised them, and showed PtSAMDC overexpression alone does nothing while co-overexpression with PtSDS1 raises spermidine. They have the enzyme and the system and **did not test putrescine activation** — that is the experiment that would fill the empty quadrant.
    - ⚠️ **The empty quadrant.** No enzyme has ever been measured that both **lacks Glu178 and fails to respond to putrescine** — which is exactly where the stramenopile prediction sits. Stage 1 can only show the measurements behave sensibly on the three populated quadrants; that is necessary, not sufficient.
  - [x] **1D: citers of Bale 2008 screened (11 Nov, 36 papers).** Nothing pre-empts the survey. Four that matter are in `BACKGROUND_factors_pocket_chemistry_survey.md`: Volkov 2016 (autoinhibition/prozyme), Velez 2013 (N-terminus essential), Patel 2025 (mobile loops 20–28 and 164–174 by multi-temperature crystallography), Lin 2024 (PtSAMDC).
  - [ ] 1E: screen citers of Michael 2011 (+ Bennett 2002, Stanley 1994, Ndjonka 2003)
  - [ ] 1F: genome vs transcriptome origin of the key non-Glu sequences
  - [ ] **1G (added 6 Nov, agreed with you): background check on the dinoflagellate AdoMetDC–ODC fusion.** Is it already reported, or just inherited from the ancestor shared with *Plasmodium*?
    - (a) PubMed / Europe PMC search: dinoflagellate (and *Perkinsus*, chromerid) AdoMetDC–ODC fusion, including dinoflagellate genome papers.
    - (b) In our own data: which dinoflagellate copies are fused, the domain order (AdoMetDC first, like *Plasmodium*?), and whether they come from genomes or transcriptomes. A transcriptome fusion could be an assembly artifact, so this pairs with 1F.
    - (c) Do the fused copies group with *Plasmodium* in our tree? Already partly answered: the Gln copy sits with apicomplexans.
    - Outcome goes in the paper as background or as a short note, never as a headline claim.
    - Note from the tree figure (6 Nov), to check in 1G: "ODC-fusion copy" labels appear in three places. Most are in the Gln clade with apicomplexans, but one *Durusdinium* fusion copy sits in the Asp clade and two *Alexandrium tamarense* fusion copies sit with green algae (Glu). These could be assembly or contamination problems, so they must be checked before saying anything about the fusion.
  - [x] 2A mixture model (C60 PMSF, 7 Nov) · [ ] 2B alternative alignment/trimming · [x] 2C Li 2019-style AU (7 Nov, not rejected) + ancestral reconstruction on the two-group trees (9 Nov) · [x] 2D spermidine-synthase tree (9 Nov: gene-specific pattern; doesn't settle one-vs-two; *Euglena* added 10 Nov: SpdSyn follows the species tree for euglenids)
    - 💡 **Add to 2C (from Li 2019):** an AU test of a "Li 2019-like" tree: kinetoplastids + cryptophytes kept together, but stramenopiles and the dinoflagellate/*Perkinsus* Asp copy moved out of their lineage. A reviewer who knows Li 2019 will ask for this test.
  - [ ] 3A Trp109 figure · [x] 3B docking validation on 1MSV (10 Nov) · [x] 3C docking into diatom, *Cafeteria* and oomycete pockets (10 Nov; plus cryptophyte, *Perkinsus*, haptophyte and controls) · [ ] 3D (optional) dock the AlphaFold model of human AdoMetDC (P17707) and compare with the 1MSV crystal, to test how far AlphaFold models can be trusted (as Begg et al. did)
- [ ] Optional: composition-robust model (mixture model) as an extra artifact check
- [ ] Redraw the structure figure with Trp109; paper "Figure 1"; 1-page summary for Marc
- [ ] Cestode (tapeworm) cases: only *Rodentolepis* Gln holds up at sequence level; low priority
- (later, optional) wet-lab partner: Lin lab (Natl Taiwan Ocean Univ) already has PtSAMDC constructs
**Status:** ACTIVE. Steps 1–6 done, sampling expanded, cross-supergroup question now **supported (one lineage)** with the stated limits. Next: HGT literature check → figures → write-up; ask Luke/Maddie about the journal last.

---

## Thread 6 — Broad survey of the 1,000-hit run 🔬 (ACTIVE — current lead)

**Scope change (15 Nov, Tui's call):** stop building on the single Glu178 hypothesis and survey all 1,000 hits with no prior. Work lives in `9_dataset_survey/` (scripts `p1`–`p8`, results in `results/`, script table in `9_dataset_survey/README.md`). Marc's `aim1_settings_used.json` untouched throughout; `models/` read only.

- **773 of 1,000 pass the fold gate** (TM ≥ 0.5 vs 1MSV chain A). Everything below uses those.
- ✅ **Annotation rescue (p2, p2b).** 326 UniProt accessions had been deleted since the models were built — 319 "not part of a reference proteome", 4 redundant sequence, 3 deleted from ENSEMBL. **211 of them are fold-passers**, and all 211 recovered domain annotation through UniParc (sequences identical to the models). **Nothing was dropped.** ⚠️ The deleted set skews toward under-studied organisms, so "no standalone copy in UniProt" arguments can only ever mean "none among *active* entries".
- ✅ **Architecture census (p2c).** 61 distinct Pfam architectures across 950 annotated proteins; 43 layouts held by a single protein. Of the fold-passers: **590 standalone · 81 with two PF01536 matches · 54 partner-then-AdoMetDC · 46 AdoMetDC-then-partner · 1 sandwiched** (`A0A8S1GM37`, HMG-box > AdoMetDC > VTT, *Caenorhabditis auricularia*). Every architecture with ≥3 proteins also appears in the InterPro-only subset, so the rescued annotation is not inventing patterns.

### Track A — fusion vs site loss ⛔ REFUTED (17 Nov)
Hypothesis: fused proteins lose the five prespecified catalytic/processing sites more often than standalone ones.

- An earlier pass using a **length-based proxy** for "fused" gave 31/523 vs 1/151, p = 0.0041. **That was an artefact.** The proxy counted long N-terminal extensions with no partner domain as fusions.
- With the proper Pfam definition: standalone **25/590 (4.2%)** vs fused **3/101 (3.0%)**, Fisher **p = 0.785**; within-phylum permutation **p = 0.717**; MH odds ratio 0.96. InterPro-only: 3.6% vs 3.6%, p = 1.
- Confounded anyway: fused models have **lower** pLDDT than standalone (75.1 vs 83.5, p = 4e-21).
- Losses concentrate in Euglenozoa (7 of 11) — the known kinetoplastid story, not a new pattern.
- ⛔ **Side observation now withdrawn (19 Nov).** The 81 two-PF01536 proteins appeared to lose sites at 16%, the highest of any class. **That was AlphaFold noise.** Once calls below pLDDT 70 are excluded the rate is **5/81 (6%) vs 3.4% standalone, p = 0.21** — 8 of the 13 original "losses" were low-confidence calls. Do not cite the 16%.
- ⛔ **Do not re-chase.** Full numbers: `9_dataset_survey/results/p5_loss_vs_fusion_report.txt`.

### Track B — the AdoMetDC–ODC fusion 🆕 (the current finding)
- **23 proteins are AdoMetDC fused to an ODC-family PLP decarboxylase.** The fusion is described in *Plasmodium* (Muller 2000, Krause 2000, Wrenger 2001, Birkholtz 2004, Hart 2016) and reported once in *Micromonas* (Green 2011, sequence only — and **that exact protein, `C1E5E7`, is in this set**).
- **The other 22 were not found in the literature.** Spread: green algae (*Myrmecia*, *Micractinium*, *Apatococcus* ×2 spp., *Cymbomonas*, *Tetraselmis*, *Chloropicon*, *Prototheca*, *Pycnococcus*, *Auxenochlorella*), streptophytes (*Klebsormidium*, *Closterium*, *Adiantum*, *Ceratopteris*, *Riccia*), red algae (*Galdieria*, *Rhodosorus*, *Gracilariopsis*), plus ***Percolomonas*** (Heterolobosea) and ***Vitrella*** (chromerid).
- **16 of 23 partners are confidently ODC** by alignment: 44–54% identity, Lys69 and Cys360 both present.
- Michael-lab papers checked specifically (the main pre-emption risk): Green 2011, Li 2023 *JBC*, Li 2026 *JBC*, Li 2023 *PNAS*, Lin & Lin 2018, Fuell 2010. Also four non-polyamine partners, all unreported, 2–3 proteins each: Rab-GAP TBC (plants), VTT (nematodes), NUBPL (Sordariomycetes), epoxide hydrolase N-terminus (Dothideomycetes).
- 🔗 Partly answers standing to-do **1G** (the dinoflagellate AdoMetDC–ODC fusion question, raised 6 Nov) from a different direction.
- Files: `results/trackB_fusion_partners.md`, `results/p6_decarboxylase_partner.tsv`.

### Track B validation (p7, 17 Nov) — survives
- **19 of 22 are genome gene models**; 3 transcriptome-only (*Percolomonas*, *Tetraselmis* GSL018, *Auxenochlorella* `A0A1D2AC22`) and set aside.
- **9 of the 19 run through the junction without a break** — 5 single-exon genes, 4 where the junction sits inside one exon. The other **10 rely on a predicted intron in the linker** (the earlier "9 / 10" in this log had the two groups the wrong way round; the 10th is the partial *Apatococcus fuscideae* `A0AAW1STH3`). ⚠️ **Corrected 19 Nov — an unbroken junction is the strong evidence, not a predicted intron.** Gene-prediction software merges neighbouring genes by calling the intergenic DNA an intron, so a predicted intron at the junction is the signature of the artefact. The sandwich protein (lead 1) was exposed exactly this way.
- **Multiple copies in one genome:** *Ceratopteris richardii* **8**; *Galdieria*, *Riccia*, *Apatococcus fuscideae* 2 each. An assembler does not repeat the same mis-join eight times.
- In every taxon except *Adiantum*, *Ceratopteris* and *Rhodosorus*, **the fusion is the only AdoMetDC-Pfam protein among active UniProt entries** — the same situation as *Plasmodium*.
- ✅ **Intron-size check done (19 Nov, lead-3 chat, gene level): none of the junction introns is oversized.** *Galdieria*'s is 53 nt against a median of 52 for that gene's other introns; none is the largest intron in its gene and several are the smallest. Compare the sandwich protein, whose two junction introns were **1,125 and 1,456 nt against a typical 39–62**. So the clearest signature of a merge is **absent** from the 23.
  - ⚠️ **Not proof.** A merge can also convert intergenic DNA into predicted *exons*, leaving normal-sized introns between them. The genome-level version of this check (junction intron vs 400 random genes from the same genome) is in `scripts/q3_fusion_gene_models.py` and has not been run yet.
  - ⛔ **"The fusion is the only AdoMetDC in that genome" is NOT evidence for a real fusion.** Retracted 19 Nov. A merge swallows the AdoMetDC gene, so no separate copy remains to be found — the sandwich protein holds *C. auriculariae*'s only AdoMetDC and is still probably a merge. This argument was asserted twice in earlier versions of this log and of `trackB_fusion_validation.md`; both are corrected.
  - [ ] **Still unchecked: second genome annotations.** Whether any of these species has another assembly or annotation lacking the fusion — the test that exposed the nematode AdoMetDC–VTT case. Also in `q3_fusion_gene_models.py`.
  - [ ] No RNA-seq or EST support across any junction has been checked, for any of the 23.
  - 💡 **The strongest surviving argument is one that does not depend on gene models at all:** all 23 carry the *same* partner, ODC — the enzyme one step away in the same pathway. Random merges produce random neighbours (the sandwich got an HMG-box and a VTT domain). Recovering the same functionally related partner across green algae, red algae, ferns, a liverwort, Heterolobosea and a chromerid is not what mis-annotation looks like. ⚠️ The one way it fails: if AdoMetDC and ODC are genomically adjacent in these lineages, repeated merges become plausible — **check gene order in genomes where the two are separate.**
- 🆕 **Four partners are structurally complete but catalytically gutted:** *Percolomonas* 2/13 ODC functional residues, *Galdieria* 1/13 (Cys360→Gly), *Rhodosorus* 2/13, *Gracilariopsis* 3/13 — against **AZIN1 9/13** (a known dead ODC paralogue) and *T. brucei* ODC 13/13. Fold TM and pLDDT (75–79) comparable to the intact partners. A dead enzyme domain fused to a live one; possibly a separate finding. ⚠️ Two of the four score nearest to AZIN1 by sequence, but that probably just means "divergent" — the antizyme system is not known outside animals, so do not read them as antizyme inhibitors.
- Linkers 41–202 aa (median 133), mostly low pLDDT. *P. falciparum* Q8IJ77 measured the same way: 359 aa.
- ⚠️ **Open pre-emption risk:** Woo et al. 2015 (*eLife*, PMC4501334, chromerid genomes) says chromerids share polyamine synthesis with *Plasmodium*; **Figure 2 supplement 3 not checked**. *Vitrella* is the organism that could be pre-empted. Also unread: Lee 2018 (*Gracilariopsis* genome), 2018 *Micractinium* genome, 2009 *Micromonas* genomes paper — all paywalled.
- File: `results/trackB_fusion_validation.md`, tables `results/p7_*.tsv`.

### p8 — does the fused AdoMetDC half keep the putrescine site? (19 Nov) — mostly NEGATIVE
Eight positions in human 1MSV numbering (13, 15, 111, 174, 176, 178, 256, 285) mapped onto 23 ODC fusions, 78 other fused, 590 standalone. Same TM-align method as p3; pLDDT ≥ 70 or the call is unresolved.

- **Glu178 is kept in all 23** (pLDDT 71–97, Cα offsets 0.27–1.18 Å). But standalone keeps it **514/579 (89%)** anyway, so there is no room for a difference: Fisher p = 0.159, permutation **p = 0.517**.
- ⚠️ **Leu13→Arg (20 of 23) and Phe111→Lys (15 of 23) are NOT a fusion signal.** Standalone proteins do the same: **133 Arg + 58 Lys at position 13** among 590. Fisher p = 5.5e-07 collapses to permutation **p = 0.667** once compared within phylum. Same for Phe111 (0.671), Phe285 (0.883), Asp174 (0.709), Thr176 (0.998). The basic "filled pocket" is a plant/algal trait, not a consequence of fusion.
- 🆕 **One survivor — Glu15.** Standalone **424/459 kept (92%)**; the 23 fusions **13/22 (59%)**, five of the nine losses to Ser. Fisher p = 4.05e-05, **within-phylum permutation p = 0.0055**. The other 78 fused proteins show nothing (p = 0.647), so it is specific to *this* fusion rather than to fusion in general. Holds in the InterPro-only subset (p = 0.043).
- ⚠️ **Lead, not a result.** Only 5 phyla contain both groups, and eight positions were tested — it barely clears correction for eight, and not for everything run that day. Needs confirmation before it is stated anywhere.
- 💡 **For the original project:** of the 65 standalone proteins that lost Glu178, **30 went to Asn** — the diatom "stuck-on" pattern from Thread 5, sitting untouched in the control group.
- Files: `results/p8_putrescine_site.tsv`, `p8_putrescine_site_tests.tsv`, `p8_putrescine_site_report.txt`.

### Gene-model worry — substantially closed, 19 Nov ✅
`q3_fusion_gene_models.py` (lead-3 chat), results in `results_81/q3_*`. Two tests, both prompted by the lead-1 sandwich protein turning out to be a merged gene model.

- ✅ **Junction-intron size, now at genome level: all 15 linker introns across 10 genes are NOT outliers.** Each junction intron was compared with the other introns of the same gene *and* with the genome-wide intron distribution (median, p95, percentile). Percentiles run 6.5 to 75.3 — several junction introns are among the *smallest* in their gene. Examples: *Galdieria* 53 nt (gene median 52, genome median 50, 66th percentile); *Rhodosorus* 73 nt (7th percentile); *Prototheca* 68 nt (6.5th percentile); *Vitrella* 79 and 98 nt (45th and 75th). Compare the sandwich protein: **1,125 and 1,456 nt against a typical 39–62**. The merge signature is absent from every one of the 23.
- ✅ **Independent annotations: where the test could be run, the fusion reproduces.** 11 species have ≥2 ENA assemblies, and 4 fusion proteins had a second *annotated UniProt proteome* to scan — `A0A087SFF4` and `A0A1D2AC22` (*Auxenochlorella protothecoides*), `A0A5B8MZ45` (*Chloropicon primus*), `A0A830HJR7` (*Pycnococcus provasolii*). **All 8 proteomes scanned return verdict = fusion.** In 4 the fused sequence is byte-identical; in the other 4 a fusion protein is present but not identical, i.e. an independently assembled and annotated reconstruction of the same architecture — the stronger result of the two.
  - 💡 Note: `A0A1D2AC22` was one of the three **transcriptome-only** proteins set aside as weak. It is now confirmed in genome-derived proteomes.
  - ⚠️ **Coverage limit:** only 4 of the 23 had a second annotated proteome. This is a spot check that passed, not a comprehensive clearance. The remaining species have one proteome each.
- ⛔ This is the test that exposed the nematode AdoMetDC–VTT fusion (*C. remanei*: one annotation with it, one without). Applied here, 4 of 4 pass.
- [ ] **Still the strongest outstanding test:** a transcript or RNA-seq read spanning the junction (tblastn against TSA / EST / SRA; **MMETSP** for the marine algae, **1KP** for the streptophytes). Start with the single-exon genes. One spanning transcript settles a protein outright. Tui runs this on the NCBI BLAST site.
- [ ] Linker conservation between related species (*Apatococcus*/*Myrmecia*/*Micractinium*; the red algae) from `p7_linkers.tsv` — intergenic DNA diverges fast, so a linker conserved at the protein level is real coding sequence. Offered by the lead-3 chat, no downloads needed.
- [ ] Gene order: in genomes where AdoMetDC and ODC are separate, are they neighbours? If scattered, the repeated-merge explanation dies outright. This is the direct test of the "same partner every time" argument.

### Literature search — Marc's ask, DONE 19 Nov ✅
Files: `results/literature_adometdc_odc_fusion.md` (49 kB), `results/literature_query_log.tsv` (every query verbatim with hit counts). Europe PMC REST + PubMed E-utilities; full text from the PMC open-access subset, searched sentence by sentence. Seven sections as specified.

**Outcome: not pre-empted, but the claim must be narrowed and one paper must be cited.**

- 🆕 **Volkov et al. 2016, *eLife* 5:e20198 (PMID 27977001, PMC5201418) — Phillips lab.** In the **legend of Supplementary file 2**, describing an NCBI NR survey done for another purpose, they report AdoMetDC fused to PLP ornithine decarboxylase (PF02784 + PF00278) in **"9 apicomplexan and 3 plant sequences"**. The three are **not named**; Supplementary file 2 itself holds only representative AdoMetDCs, not the fused sequences. "Plant" in NR = Viridiplantae, which includes *Micromonas*, so two of the three are probably Green 2011's. ⚠️ ***Vitrella* is a chromerid and may sit inside their 9 apicomplexan — do not claim it as new without checking.** They also count 35 PF08132 leader fusions in plants; this dataset has 43.
  - Mitigation: the survey was NR as of 2016. *Riccia* (2025), *Micractinium* (2026), *Closterium*, *Cymbomonas*, *Pycnococcus* (2023), *Ceratopteris* (2019) and others post-date it. Red algae, *Percolomonas*, the ferns and the liverwort are outside "plant"/apicomplexan entirely.
  - ➡️ **Required wording change:** no longer "unreported outside *Plasmodium* and *Micromonas*", but "reported in *Plasmodium*, in two *Micromonas* species (Green 2011), and noted without identification in three unnamed plant sequences (Volkov 2016)".
  - 📧 Phillips would know which three. **Marc's call, not Tui's** — asking signals the project to the lab most able to pre-empt it.
- ✅ **Woo et al. 2015 (*eLife*, PMC4501334, chromerid genomes) — cleared.** Main text, appendices, Figure 2 and all 5 supplements, all 16 supplementary tables searched. Data are presence/absence of EC numbers per species with no gene IDs and no domain layout. **The fusion is not shown.** This was the largest pre-emption risk and it is closed.
- ✅ **The 23 accessions, 24 locus tags and 6 protein IDs as literal strings: 0 hits** in Europe PMC and PubMed, except *Micromonas* XP_002502359 → Green 2011. The 11 UniProt-linked papers were read in full where open access (8 of 11); none mentions polyamine enzymes.
- ✅ ***Plasmodium* literature (219 papers, core plus all citers of 11 key papers): none surveys outside apicomplexans.** They call the arrangement unique to *Plasmodium* — Birkholtz 2004 confirms it in at least three species. **This is the statement the work revises.**
- ✅ **Michael AJ 2016, *Biochem J* 473:2315** — PDF obtained and full text searched (19 Nov): no *Micromonas*, no AdoMetDC–ODC fusion, nothing on distribution. Clean.
- 🆕 **Michael AJ 2017, *Biochem J* 474:2277 "Evolution of biosynthetic diversity"** — PDF obtained and full text searched (19 Nov). **This is the paper to quote in the introduction.** In its gene-fusion section Michael states that fusions of eukaryotic polyamine biosynthetic enzymes retaining all constituent functions are **rare**, and gives the eukaryotic ODC–AdoMetDC fusion of parasitic *Plasmodium* as the **best-known example** (refs Müller 2000, Krause 2000, Wrenger 2001). In the "Other genes" paragraph he lists polyamine-enzyme fusions found in a **limited number of species** and, for eukaryotic class-2 AdoMetDC fused to AR-fold ODC, names **only *Micromonas*** (ref Green 2011), describing the set as "similarly phylogenetically constrained".
  - ➡️ So the field's leading authority, in a 2017 review *devoted to gene fusion in polyamine metabolism*, catalogued this fusion as *Plasmodium* + one green alga. That is the framing for the paper's opening, and it is independent confirmation that nobody followed up Volkov's three sequences.
- [ ] **Still to obtain (Conn library, Tui):** Michael AJ 2011, *Methods Mol Biol* 720:139 (abstract names gene fusion — highest remaining risk); Fuell et al. 2010, *Plant Physiol Biochem*, doi 10.1016/j.plaphy.2010.02.008 (plants and algae). Lower priority: Michael AJ 2016 *JBC* (PMC4946907); Phillips 2018 *JBC*; Pegg 2018. Search each for "ornithine", "fusion", "Micromonas".
- [ ] Also unread: Worden 2009 *Science* SOM; Arriola 2018 *Plant J*; the Lee 2018 *MBE* supplementary zip. bioRxiv's own search returned HTTP 403.
- ⛔ **Do not email Michael or Phillips for PDFs.** Their labs are the pre-emption risk; a request signals the project. Library and ILL only.

### 20 Nov — both halves, the expanded search, and the red-algal pattern

**p9: does the fused AdoMetDC half still self-process? ✅ Yes, and it is indistinguishable from a standalone enzyme.**
Marc's five sites mapped onto all 23 fusions and all 590 standalone proteins, p8 method, pLDDT ≥ 70 gate.
- Ser68 22/22 vs 501/510 · His243 22/22 vs 524/533 · Cys82 22/22 vs 528/537 · Glu247 21/21 vs 517/526 — **all Fisher p = 1**.
- Ser229 21/22 (one loss: *Riccia* S229→F) vs 537/547; Fisher p = 0.328, within-phylum permutation p = 0.084. ⚠️ Driven entirely by that single loss. Noise at n = 1, footnote only.
- ➡️ These are working AdoMetDC enzymes. Fusion does not compromise the AdoMetDC side at all, which is what makes the ODC-side degradation stand out.
- ✅ **Cut-off robustness check done.** Six calls sit within 1.5 pLDDT points of 70. Flipping them changes nothing for any degraded protein; only *Closterium* (D361, 68.6) and *Vitrella* (S200, 69.2) move, and both move *toward* intact. The pattern does not depend on where the cut-off sits.
- ⚠️ **Methodological correction applied:** the first pass gated the AdoMetDC half at pLDDT ≥ 70 but not the ODC half. Regated, the "7 gutted ODC halves" became **4 confidently degraded** (*Percolomonas*, *Galdieria*, *Rhodosorus*, *Gracilariopsis*); *Prototheca*, *Auxenochlorella* A0A1D2AC22 and *Apatococcus* A0AAW1STH3 are **truncated or badly modelled, not degraded** (partner 65–254 aa, or 61% of residues below pLDDT 70). Files: `p9_self_processing.md`, `p9_both_halves.tsv`, `p9_odc_sites_gated.tsv`, `p9_self_processing_tests.tsv`.
- ✅ **The four degraded ones are real:** lost residues at pLDDT 70–91, spread over 3–4 of 7 site groups, not clustered in a low-confidence patch.

**p10: the sample is much larger than the 1,000-hit run showed.**
UniParc searched on Pfam features with **no structure filter**: `(PF01536 OR PF02675) AND (PF02784 OR PF00278)` → **311 sequences, 288 new**. ⚠️ **Marc's TM gate was NOT changed; this search uses no structure at all.** That means two tiers of evidence and they must be labelled separately in any write-up: 23 fold-verified (TM ≥ 0.5), the rest on domain annotation only.
- Alveolata 114 (~70 *Plasmodium* strains — redundancy, not news; **~40 dinoflagellates**; *Chromera*; a ciliate) · Viridiplantae 92 · Bacteria 66 (64 in the reverse ODC→AdoMetDC order, a different architecture) · **Rhodophyta 10** · **Glaucophyta 1** (*Gloeochaete*) · **Stramenopiles 1** (brown alga *Heribaudiella*) · **CRuMs 3** · Metazoa 1 (probable parasite contaminant).
- ➡️ The fusion is present in **all three Archaeplastida lineages**, in alveolates well beyond *Plasmodium*, in a stramenopile and in CRuMs. "22 unreported proteins" is no longer the right headline.

**The red-algal pattern — restated precisely and strengthened.**
- The original phrasing ("every green-lineage fusion has an intact ODC half") is **wrong**: three green fusions are truncated/incomplete and four undetermined.
- What holds: **every confident in-place substitution at the 13 ODC functional residues is in a red alga or *Percolomonas*** — 8 of 8 proteins. **None of the 18 green-lineage fusions has one.**
- Red-algal fusions with models are now **7 of 7 degraded**, spanning Cyanidiophyceae, Bangiophyceae, Florideophyceae and Stylonematophyceae.
- **K69→R** (the PLP Schiff-base lysine) is shared across three of those classes. ⚠️ K→R is a common substitution; 13 sites cannot exclude convergence.

**The tree (NJ, Kimura, 1000 bootstraps, 197 sequences, 344 columns) — suggestive, not settled.**
- **Red-algal fusions form a clade** (the only other sequence inside it is another *Rhodosorus*).
- ***Percolomonas* does not group with them** → its degradation is an independent event.
- **Degraded halves are not a clade** → the ODC half has died **at least twice independently**.
- ⚠️ **Do not over-read this.** Bootstrap is 0% on every set tested, including the apicomplexan one, which is certainly monophyletic — so the tree carries no statistical support. And the alignment dropped every sequence with a PLP span under 250 aa, which removed *Galdieria partita* (both copies), *G. sulphuraria* and *Porphyra* — i.e. the filter is biased against exactly the degraded sequences the question is about.
- ⛔ **Decision (Tui, 20 Nov): do not chase a maximum-likelihood tree.** "At least two independent events, weak support" is enough; the mechanism work is the priority. Files: `p10_report.md`, `p10_odc_tree.nwk`, `p10_tree_support.tsv`, `trackB_red_algal_pattern.md`.

**⚠️ Correction: the "*Gracilariopsis* has no working ODC" claim is too strong.**
- Guo et al. 2026 (*Plant Physiol Biochem* 238:111700, PMID 42679707) show ODC activity is **functionally required** in *Pyropia yezoensis*: DFMO blocks conchospore development, *PyODC* siRNA silencing arrests it, exogenous putrescine partially rescues.
- The accession the paper prints for *PyODC*, **KAK1862116.1, is not a decarboxylase at all** (576 aa, PF10017 histidine methyltransferase). The gene NCBI calls ornithine decarboxylase is **KAK1862958.1 = A0ACC3BYM5** (582 aa, PF02784 only, separate chromosome, 39% identical to the fusion's ODC half).
- *P. yezoensis* therefore has a fusion **and** a separate ODC gene — same arrangement as *Galdieria* and *Rhodosorus*.
- 🆕 **But neither carries both catalytic residues:** the fusion's ODC half has **K69→R, C360 kept** (5/13); the standalone has **K69 kept, C360→Asn** (10/13).
- ➡️ Since that organism demonstrably has ODC activity, **a single substitution at a human-ODC catalytic position does not prove the enzyme is dead.** The 13-residue criterion is a human-centric yardstick applied across >1 billion years of divergence.
- ⚠️ *Gracilariopsis*' standalone A0A2V3IYX3 was called non-functional on exactly that basis (C360→Arg, 10/13). **Soften to: "carries a substitution at a position catalytic in human ODC."** The *fused* halves remain a different matter — they carry 5–7 confident substitutions each, not one.

**Prepared, not yet run: the dimer models (Marc's "how does the fusion work").**
`data/p11_dimer_fasta/` holds 7 ready-to-submit AlphaFold Server jobs (homodimers): *Klebsormidium* (ODC half intact), *Gracilariopsis* (degraded, no other ODC), *Galdieria* (degraded, has a standalone ODC), *Micromonas*, *Plasmodium* Q8IJ77, plus standalone ODC and standalone AdoMetDC controls, with a README on what to look for. Premise: ODC's active sites are built across the subunit interface (Fuell 2010), AdoMetDC is also a dimer, so one dimerisation of a fusion protein may serve both halves — and a dead ODC half may still be load-bearing as an interface.

### Marc's directions (verbal, 19 Nov, relayed by Tui)
Phylogeny: yes · look at putrescine sensing nearby, and whether anything else plays that role in other AdoMetDCs · the fusion presumably arose because it helped, so work out how · how common is fusion generally · **use a fusion structure as the reference structure** for a new comparison · look at the important parts of the structure.

### Next on this thread
- [ ] **Thorough literature search (Marc's ask).** The fusion as a concept rather than organism by organism; Michael's complete publication list and his pathway-diversity reviews; the 23 accessions as literal strings; UniProt cross-references; preprints; Woo 2015 figure supplements; abstracts and supplementary data of the three paywalled genome papers. Report the query list, not just the findings.
- [ ] **Does the AdoMetDC half still self-process?** Ser68, Ser229, His243, Cys82, Glu247 across all 23. If it is dead, the whole story reads differently. Cheap — reuses the p8 alignments.
- [ ] **Geometric putrescine-site scan** (Marc's "sensing nearby"). Any acidic side chain within a fixed distance of the 1MSV putrescine position after superposition, **regardless of sequence position** — the test p8 could not do. *T. cruzi* is the positive control (shifted site). Fix the cutoff before looking.
- [ ] **Fusion base rate** (Marc's "how common"). Same census on ODC, spermidine synthase and SAM synthetase. Without it, "13% fused" has nothing to compare against.
- [ ] **Inter-domain PAE from AFDB first**, then — only where PAE supports it — distance from the partner's Schiff-base lysine to the putrescine position, buried interface area, and whether the orientation is the same across all 23. Same-orientation vs different answers the one-fusion-or-many question structurally.
- [ ] **Phylogeny of each half separately** (Marc's ask); compare with the orientation result above.
- [ ] **New FOLDTRACE run with a fusion as the query** — its own run with its own settings file, **never** under the locked aim1 protocol.
- ➡️ **Lead 1 (the sandwich protein) and lead 3 (the 81 two-PF01536 proteins) handed to a separate chat:** `1_reports/HANDOFF_81_double_domain.md`. Prefix `q`, results in `9_dataset_survey/results_81/`. Never run two Claude Code jobs against the repo at once. Whatever it concludes about the 81 comes back here, since they are their own class in the Track A statistics.
  - ✅ **Lead 1 (sandwich) CLOSED 19 Nov — gene-prediction artefact.** A complete, normal-length AdoMetDC (Met236–Glu601, 366 aa, 75% identical to *C. elegans* SMD-1, start and end on exon boundaries) sits inside the 865-aa model; the HMG-box and VTT hang off across the two **largest** introns in the model. No other *Caenorhabditis* has the HMG-box; nothing in the literature; all five sites kept at pLDDT 92–98. Report: `9_dataset_survey/results_81/q1_sandwich_report.md`.
  - 🔄 **Lead 3 (the 81) closing, 19 Nov (test sample, 12 of 87).** Most are **one domain that Pfam cut in two**: the two matches cover opposite halves of the human protein (e.g. 6–164 + 169–326) with the cut near human ~165–172, consistent with the internal symmetry of the AdoMetDC protomer reported by Ekstrom 1999. Only *Dillenia* (793 aa) has two genuinely separate units, and it is **not** a prozyme — both units keep all five sites, the halves are 94% identical and do not touch, and inter-unit PAE is 27 Å so the arrangement is unsupported. The 81 are not all pairs: 76 have two Pfam matches, 5 have three. No single-chain enzyme+prozyme found in Europe PMC.
  - ⚠️ **Knock-on for Track B:** the three nematode **AdoMetDC–VTT** fusions are probably the same artefact — *Caenorhabditis remanei* has one genome annotation carrying the fusion and another that does not. Withdrawn from the non-polyamine partner list in `trackB_fusion_partners.md`. **That double-annotation test has not been run on the 23 ODC fusions and should be.**

**Status:** ACTIVE. Track B is the finding; Track A refuted; p8 negative with one lead. The literature search gates everything else.

---

### 23 Nov — p14 tblastn, p17 standalone decarboxylases, and three corrections

**p14 tblastn, query 06 (*Percolomonas cosmopolitus*, 191 aa junction query), run manually on blast.ncbi.nlm.nih.gov against TSA.** `HBGD01001067.1` (MMETSP0759, 2,690 nt) aligns across all 191 residues at 100% identity, positions 938–1510, one unbroken block, E = 7e-135. The junction is carried on a single RNA molecule, and the query construction is correct.

⚠️ **The hit is circular and was always going to be.** MMETSP0759 is the assembly the UniProt entry was built from, so 100% identity was the only possible outcome. The same applies to queries 07 and 08 (*Tetraselmis*, *Auxenochlorella*) — transcriptome-only species whose entries derive from the very TSA assemblies tblastn would search. Queries 01–05 are single-exon genes, already immune to the gene-merge artefact by construction. **Conclusion: the p14 set cannot produce independent evidence for any of its eight proteins and should not be completed.** The one exception is a genus-level search returning a hit in a *different species*, which is a genuinely independent organism; `02_Chloropicon` is being run that way because TSA appears to hold no *C. primus* data.

The informative part of the run was incidental. Four other *Percolomonas* transcripts appear as partial hits — three aligning only to the ODC half at 30–36% identity, one only to the AdoMetDC half. The three ODC-like ones all read `YAVKCNPNK`, human ODC's PLP Schiff-base motif with K69 present, while the fusion reads `YTVSCN` (K69→S, p9 residue 411, pLDDT 82, offset 1.3 Å). Two sit in transcripts of 1,521 and 1,383 nt, too short to hold both domains, so they are standalone decarboxylases; one is from a second library, MMETSP0758.

**p17 (`results/p17_standalone_odc.tsv`, `p17_report.md`): the standalone decarboxylases of the four degraded-fusion species.** Sequence alignment to human ODC P11926 — *not* the p9 structural test — kept in its own file and not merged into any p9 table.

| species | standalone ODC-family copies | best copy |
|---|---|---|
| *Rhodosorus marinus* | 3 intact + 1 truncated + 1 AZIN-like dead | 13/13, TM 0.90, nothing lost |
| *Galdieria partita* | 2 | 11/13, TM 0.91, loses only Y331 |
| *Percolomonas cosmopolitus* | 2 | 12/13 |
| *Gracilariopsis chorda* | 1 | 10/13, **C360→Arg**, TM 0.87 |

Identity to human ODC is bimodal and separates the families cleanly: ~42–48% is ODC-family, ~25–31% is diaminopimelate decarboxylase. `p7_same_species.tsv` counts the wider PF02784/PF00278 category, which is why its "standalone decarboxylase" column overstates the ODC count.

**⚠️ This weakens the broad pseudoenzyme claim.** *Rhodosorus*, *Galdieria* and *Percolomonas* each carry one to three fully intact standalone ODCs. Against that background a dead fused ODC half is explained by redundancy-driven relaxed selection, the default fate of a dispensable duplicate. *Rhodosorus* also carries a second unrelated dead ODC-like paralog (nearest match human AZIN1, 3/13), so free decay of ODC copies is demonstrably common in that genome. The earlier argument that "only one half rots, so the gene is under selection" does not survive: the AdoMetDC half is the functional gene and the ODC half is a passenger, and a passenger rotting under a functional host is exactly what the null predicts.

**What survives is narrower: *Gracilariopsis chorda*.** It is the only one of the four with no intact ODC anywhere. Its single standalone ODC-family protein (A0A2V3IYX3) has lost C360 (→Arg, structurally confirmed, TM 0.87, 10/1/2). Its fusion (A0A2V3IZY9) has lost K69 (→Arg, pLDDT 86) but **kept** C360 (pLDDT 88, CA offset 3.5 Å) and D361. In a real ODC dimer K69 comes from one chain and C360 from the other, so the organism's two broken copies are broken at opposite half-sites.

Caveats that must travel with that observation: n = 1 species; the fusion's C360 sits at 3.5 Å against a 4.0 Å gate; two separate gene products would have to be co-expressed and heterodimerise, for which there is no evidence; and *Gracilariopsis* carries an agmatine iminohydrolase (A0A2V3IZE8) plus nitrilase-family NCPAH candidates — the arginine route to putrescine — so it may not need ODC at all, and both copies may be decaying for the ordinary reason.

**The natural experiment this sets up.** Same fusion architecture, same lineage, one species with three spare ODCs and one with none. dN/dS on the ODC half, *Gracilariopsis* against *Rhodosorus* and *Galdieria*, is the test that separates maintenance from decay. Nucleotide work only — no structures, no AlphaFold. Queued as p18.

A second job worth more than the remaining p14 queries: a ColabFold **heterodimer** of the *Gracilariopsis* fusion A0A2V3IZY9 with its standalone ODC A0A2V3IYX3 (1,163 residues), asking whether K69 of the standalone and C360 of the fusion come together as they do in 1D7K.

**Corrections made today:**

1. ❌ Claude claimed the *Percolomonas* TSA search "filled in a blank UniProt had lost." It did not. `p10_species_enzymes.tsv` already listed two standalone ODC-family proteins for *Percolomonas* (A0A6U1HKK3 12/13, A0A6U0M0T4 12/13) plus a standalone AdoMetDC. The BLAST confirmed results already in the repo.
2. ❌ p17 was largely redundant with p10, which already held the ODC-vs-DAPDC classification and the structural reads. **Read `p10_species_enzymes.tsv` before commissioning any further question of this kind.**
3. ✅ The *Gracilariopsis* conflict is resolved **in favour of the earlier note**: one standalone ODC-family protein, carrying C360→Arg. The apparent contradiction with `p7_same_species.tsv` was a category difference, not an error.

---

## ⛔ Dead ends & resolved — do NOT re-chase
- ⛔ **p14 tblastn junction searches (23 Nov)** — cannot produce independent evidence. Queries 06–08 are transcriptome-only species whose UniProt entries derive from the same TSA assemblies being searched (06 returned a 100% self-hit); queries 01–05 are single-exon genes already immune to the gene-merge artefact. Only a genus-level hit in a *different species* would be informative. Do not complete the set.
- ⛔ **"Only one half of the fusion rots, so the gene must be under selection"** — does not follow. The AdoMetDC half is the functional gene; the ODC half is a passenger, and a passenger decaying under a functional host is what relaxed selection predicts. Superseded by the redundancy finding in p17.
- ⛔ **"The degraded ODC domain keeps its fold, therefore it is maintained"** — near-unfalsifiable from AlphaFold DB models, which fold family members into the family fold regardless. Do not use as evidence without a dN/dS or experimental argument.
- ⛔ **"1I7B chains A and B are the dimer"** — they are not; they are the **α and β chains of one processed protomer**. The real dimer is 1MSV chains A+B, with 1JEN as the apo processed-form check. Claude made this α/β mix-up twice (see the earlier entry below); check it every time a processed AdoMetDC structure is described by chain ID.
- ⛔ **Overclaiming Bodo as fully novel** → corrected to honest incremental framing. (Settled; keep the honest line.)
- ⛔ **Full Pcon sequences via UniProt/NCBI** → not there; TriTrypDB only.
- ⛔ **Running AlphaFold / ColabFold / homology modeling in this workspace** → no GPU, no modeller/pymol, folding APIs blocked. Use AF Server.
- ⛔ **My α/β pyruvoyl-chain mix-up** → fixed; verified absent from reports and GitHub.
- ⛔ **Fusion → catalytic-site loss (Track A, 17 Nov):** the p = 0.004 result came from a length-based proxy for "fused". With the Pfam definition it is p = 0.785 (permutation 0.717). Settled; see Thread 6.
- ⛔ **WOW #3, eukaryote-type (class-2) AdoMetDC in bacteria:** already published (Li et al. 2019 *Biochem J*, PMID 31467246). Found 5 Nov in chunk 1A. Don't pursue as new.
- 🔻 **Low-priority WOW ideas** (not dead, but not lead projects): *Palpitomonas* single-point engine loss (risky); HGT trio; invariant Gly83/Pro246 (supporting figure only).

## 📂 Where everything lives
- **Delivered to you** (`outputs/`): `AdoMetDC_everything_we_know.md`, `dead_paralog_literature_review.md`, `activation_interface_mapping_Bodo.md`, `AlphaFold_modeling_package.md` (+ `alphafold_inputs/`), `WOW_findings_publication_assessment.md`, `RESEARCH_PACKAGE_for_Thursday.md`, `Marc_factcheck_CHANGELOG.md`, `corrected_reports/`, and this log.
- **My working notes** (session workspace, not in your folder): structural notes, PDB table, lit-review notes, interface-mapping results + scripts, alignments.
- **Paywalled papers — status checked 5 Nov via Europe PMC:**
  - ✅ **Both now read (PDFs from you, 5 Nov):** Michael 2011, *FEBS Lett* 585:2627–2634 (open access after all), and Flegontov et al. 2013, *Curr Biol* 23:1787–1793. Neither pre-empts us (details in Threads 4 and 5).
  - Nice-to-have, lower priority: Ribeiro 2019 (Sci Signal); Murphy 2017 (Curr Opin Struct Biol); Willert & Phillips 2012 (Trends Parasitol); Kim & Barrett 2013; Gupta 1987 & Pösö 1975 (old putrescine-sensitivity biochemistry, for validation).
  - ⛔ **Not paywalled after all:** Jackson 2016 (Bodo genome) — open access, PMC4728078, already read.
