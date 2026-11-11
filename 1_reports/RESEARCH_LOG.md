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
  - [ ] 1D/1E: screen citers of Michael 2011 and Bale 2008 (+ Bennett 2002, Stanley 1994, Ndjonka 2003)
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

## ⛔ Dead ends & resolved — do NOT re-chase
- ⛔ **Overclaiming Bodo as fully novel** → corrected to honest incremental framing. (Settled; keep the honest line.)
- ⛔ **Full Pcon sequences via UniProt/NCBI** → not there; TriTrypDB only.
- ⛔ **Running AlphaFold / ColabFold / homology modeling in this workspace** → no GPU, no modeller/pymol, folding APIs blocked. Use AF Server.
- ⛔ **My α/β pyruvoyl-chain mix-up** → fixed; verified absent from reports and GitHub.
- ⛔ **WOW #3, eukaryote-type (class-2) AdoMetDC in bacteria:** already published (Li et al. 2019 *Biochem J*, PMID 31467246). Found 5 Nov in chunk 1A. Don't pursue as new.
- 🔻 **Low-priority WOW ideas** (not dead, but not lead projects): *Palpitomonas* single-point engine loss (risky); HGT trio; invariant Gly83/Pro246 (supporting figure only).

## 📂 Where everything lives
- **Delivered to you** (`outputs/`): `AdoMetDC_everything_we_know.md`, `dead_paralog_literature_review.md`, `activation_interface_mapping_Bodo.md`, `AlphaFold_modeling_package.md` (+ `alphafold_inputs/`), `WOW_findings_publication_assessment.md`, `RESEARCH_PACKAGE_for_Thursday.md`, `Marc_factcheck_CHANGELOG.md`, `corrected_reports/`, and this log.
- **My working notes** (session workspace, not in your folder): structural notes, PDB table, lit-review notes, interface-mapping results + scripts, alignments.
- **Paywalled papers — status checked 5 Nov via Europe PMC:**
  - ✅ **Both now read (PDFs from you, 5 Nov):** Michael 2011, *FEBS Lett* 585:2627–2634 (open access after all), and Flegontov et al. 2013, *Curr Biol* 23:1787–1793. Neither pre-empts us (details in Threads 4 and 5).
  - Nice-to-have, lower priority: Ribeiro 2019 (Sci Signal); Murphy 2017 (Curr Opin Struct Biol); Willert & Phillips 2012 (Trends Parasitol); Kim & Barrett 2013; Gupta 1987 & Pösö 1975 (old putrescine-sensitivity biochemistry, for validation).
  - ⛔ **Not paywalled after all:** Jackson 2016 (Bodo genome) — open access, PMC4728078, already read.
