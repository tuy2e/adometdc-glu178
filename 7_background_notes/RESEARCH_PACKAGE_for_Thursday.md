# AdoMetDC / polyamine-pseudoenzyme research package
**Prepared overnight 3–4 Nov 2026. For Tuy's review before the Thursday meeting with Prof. Zimmer.**
*Everything here is for you to check and decide on. Nothing was sent to Marc, nothing on GitHub/the repo was touched.*

---

## 0. One-paragraph summary (what you can say Thursday)
Using FOLDTRACE + database/sequence analysis, I mapped two polyamine "pseudoenzyme" systems — the **AdoMetDC prozyme** and the **deoxyhypusine synthase (DHS) dead paralog** — across the kinetoplastids. Both are present not only in the parasitic trypanosomatids (already known) but in the **free-living *Bodo saltans*** and the **most basal trypanosomatid *Paratrypanosoma confusum***. Phylogenetics shows each dead paralog has a **single ancient origin** (dead clades 100/100, *Bodo* basal), so these regulatory systems **predate parasitism**. The dead copies are **conserved/maintained** (not decaying), which raises a clean, testable question: *what do they do in a free-living organism that has no parasitic lifestyle to regulate?*

---

## 1. WHAT'S SOLID (verified, reproducible — Marc can check every accession himself)

### AdoMetDC active + dead ("prozyme") pairs
- ***Bodo saltans*** (free-living): **A0A0S4JHX8** (active, 381 aa, full catalytic set) + **A0A0S4JRG8** (dead, 331 aa). Consecutive gene IDs (BSAL_30285 / 30280).
- ***Paratrypanosoma confusum*** (most basal trypanosomatid): **PCON_0025430** (active, 390 aa) + **PCON_0025440** (dead, 339 aa). **Syntenic** — same scaffold (Pcon_scaffold00006_3), same strand, ~1.7 kb apart; both `is_pseudo=false`.
- Trypanosomatid parasites (known): active + short "proenzyme-like" paralog in T. brucei, T. cruzi, Novymonas, Leptomonas (×2), Porcisia, Angomonas, Strigomonas.

### The dead copies really are catalytically dead (definitive residue table, human numbering)
- All **16 active** copies (incl. 6 eukaryotic outgroups): keep E8, E11, Ser68 (pyruvoyl), **Cys82** (catalytic), Ser229, **His243** (base), Glu247.
- All **9 dead** copies: **lost the catalytic Cys82** (→A/D/S, never C) AND **His243** (8/9 → Arg). Shared loss of Glu8 too.

### Single ancient origin (IQ-TREE2, model LG+G4, 1000 UFBoot + 1000 SH-aLRT)
- AdoMetDC dead/prozyme clade = **monophyletic, 100/100**. Branching order **Bodo (basal) → Paratrypanosoma → parasites** = tracks the species tree. Long-branch artifact removed by adding eukaryotic outgroups.
- DHS: the "inactive DHS" dead clade (Bodo + T. brucei DHSp + Leishmania DHSL20) = **monophyletic 100/100, Bodo basal.** Same pattern, second enzyme.

### The dead copies are maintained, not decaying
- Dead-vs-dead mean pairwise identity **48%** (34–80%), 40% of columns conserved ≥70% — far above random decay (~10–15%). → kept under selection → likely functional.

**Figures already delivered in chat:** `prozyme_IQTREE.png`, `DHS_tree.png`, `prozyme_phylogeny.png`, `adometdc_family_tree.png`.

---

## 2. HONEST NOVELTY LINE (read this before you say anything to Marc)

**Already published (NOT new):**
- The *concept* that prozyme/regulatory machinery predates parasitism and exists in free-living *Bodo* — shown for the **PRMT (arginine methyltransferase) prozyme** (bioRxiv 2019 "Structural Basis of PRMT…"; Lorenzon 2022, *ACS Infect. Dis.*), and the *antioxidant* system (Beltrame-Botelho 2016).

**Appears UNREPORTED (checked; likely still yours):**
- The **AdoMetDC prozyme** *and* the **DHS dead paralog** *specifically* in *Bodo* and *Paratrypanosoma*, with the **dated single-origin phylogeny**.
- Papers I read that do NOT report it: Opperdoes 2016 (full text — no AdoMetDC at all), Marchese 2018 (describes the T. brucei prozyme, cites Bodo metabolism, does not say Bodo has it), Beltrame-Botelho 2016 (antioxidant genes only).

**Still to confirm (paywalled — ILL):** Jackson 2016 (*Bodo* genome, *Curr. Biol.*; its indexed gene-family story is "bodonin" surface proteins, "prozyme" not indexed — low risk) and Flegontov 2013 (*Paratrypanosoma* genome).

**Honest framing for Marc:** *not* "I discovered prozymes are ancient" (known). Rather: *"Building on the PRMT-prozyme precedent, I used FOLDTRACE to show the AdoMetDC and DHS pseudoenzyme pairs show the same pattern — present in free-living Bodo and basal Paratrypanosoma, single ancient origin by phylogeny — and the dead copies are conserved, raising the question of their function in a free-living organism."*

---

## 3. PATH TO PUBLICATION (research gaps, ranked)

**Realistic home:** co-authorship on Marc's FOLDTRACE paper (this is an ideal case study — the method catching + dating two real pseudoenzyme systems), or a short standalone evolutionary note. Not a solo high-impact paper — but a real, legitimate contribution.

**Gap 1 — The free-living pseudoenzyme-function paradox (strongest; your question).**
Open: the dead paralog's *function* has only been tested in *parasites* (boosts the enzyme ~1000×). What does it do in free-living *Bodo*? (a) still activates → the "weak-enzyme-needs-dead-partner" logic is ancestral; or (b) *Bodo*'s enzyme works alone → the dead copy does something else.
Computational: check if *Bodo*'s dead copy keeps the activation interface (from the solved T. brucei complex). Lab: express *Bodo* enzyme ± dead paralog, measure activity → **enzyme-kinetics lab**.

**Gap 2 — What does the dead protein conserve, and why.**
The conserved surface of a dead enzyme reveals its real job. Map it (I started: dead copies are 48% conserved — maintained). Lab: mutate the conserved non-catalytic surface, test loss of function → **structural/mutagenesis lab**.

**Gap 3 — A pseudoenzyme-activator "rule" across the polyamine pathway?**
PRMT, AdoMetDC, DHS all recruited dead-paralog activators. Is it a repeated structural strategy (e.g., obligate hetero-oligomers)? Comparative fold analysis across ODC/SpdS too.

**Gap 4 — Moonshot:** FOLDTRACE scan of a fold family for an unreported lineage-specific active-site loss.

---

## 4. READING LIST (for ILL — mark these to loan)

**Open-access (I've read / can read — no ILL needed):**
- Opperdoes et al. 2016, J Eukaryot Microbiol 63:657 (you already have the PDF).
- Beltrame-Botelho et al. 2016, Evol Bioinform 12:263 (PMC5100842).
- Marchese et al. 2018, Pathogens 7:36 (PMC6027508).
- Willert & Phillips 2007 PNAS (PMC1895940); 2008 PLoS Pathog (PMC2562514); 2009 Mol Biochem Parasitol (PMC2730992).
- Volkov et al. 2016 eLife (PMC5201418) — the prozyme complex structure.
- Patel 2018 PLoS Pathog (PMC6221367); Nguyen 2013 JBC (PMC3663545, the DHS dead-paralog paper); Xiao 2013 (PMC3771639); Velez 2013 (PMC3576127).

**Needs ILL (paywalled — please loan):**
1. **Jackson et al. 2016, *Current Biology* 26:161–172** — *Bodo saltans* genome ("Kinetoplastid Phylogenomics…"). *Priority — novelty check.*
2. **Flegontov et al. 2013, *Current Biology*** — *Paratrypanosoma confusum* genome. *Priority — novelty check.*
3. Willert & Phillips 2012, *Trends Parasitol* 28:66–72 (review).
4. Kim & Barrett 2013, *Mol Microbiol* 88:841 (commentary).
5. Lorenzon et al. 2022, *ACS Infect. Dis.* (PRMT prozyme incl. Bodo) — check access.
6. Tanifuji et al. 2017, *Sci Rep* (bodonid genome; may be open — check).

---

## 5. WHAT I DID NOT / COULD NOT DO (honest)
- Did **not** confirm *Bodo*'s two AdoMetDC genes are physically adjacent (only inferred from gene IDs). TriTrypDB coordinate lookup needed — **you can do this** (search BSAL_30280 & BSAL_30285, check scaffold/positions).
- Used **partial** *Paratrypanosoma* sequences (reconstructed from your BLAST). Full sequences (download from TriTrypDB) would let me finalize the residue table.
- dN/dS (definitive "under selection" test) needs coding DNA — not yet run.
- Did not read the two paywalled genome papers (ILL).
