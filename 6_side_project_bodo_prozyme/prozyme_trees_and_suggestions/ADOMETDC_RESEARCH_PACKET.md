# AdoMetDC research packet — for Thursday + your own publication
Prepared overnight 3–4 Nov 2026. Everything here is either verified against primary data (marked ✓) or flagged as open/uncertain. Read the "HONEST novelty line" first so you don't overclaim to Marc.

---

## 0. The one thing to read first — HONEST novelty line

**The general idea is NOT new.** That a prozyme / dead-paralog regulatory system exists in the free-living kinetoplastid *Bodo* and predates the parasites is already published — for a *different* enzyme, the **PRMT (arginine methyltransferase) prozyme** (bioRxiv 2019 "Structural Basis of PRMT…"; Lorenzon 2022, *ACS Infect. Dis.* — both explicitly include *Bodo saltans*). So "prozymes are ancient / predate parasitism" is established.

**What appears unreported (your contribution, pending 2 checks):** the same thing shown *specifically* for **AdoMetDC** and **DHS**, extended to **free-living *Bodo* and basal *Paratrypanosoma***, with a **dated single-origin phylogeny** and the **conserved-surface analysis**. I verified Opperdoes 2016 (the *Bodo* metabolism review) never mentions the AdoMetDC prozyme, and Marchese 2018 discusses only the *T. brucei* prozyme. So nobody seems to have done AdoMetDC/DHS the way the PRMT people did their enzyme.

**Honest framing to use with Marc (defensible, can't be knocked down):**
> "Building on the PRMT-prozyme precedent, I used FOLDTRACE to show the AdoMetDC and DHS pseudoenzyme pairs are also ancient — present in free-living *Bodo* and basal *Paratrypanosoma*, single origin by phylogeny — and I mapped the conserved surface the dead copies still protect."

That is an **incremental extension of known biology**, which is exactly what makes a strong *methods-paper case study* (ideal for Marc's FOLDTRACE paper) and a legitimate contribution. It is **not** a solo "nobody knew this" bombshell — don't pitch it that way.

**2 novelty checks still open (need you / full-text):**
1. The trypanothione/polyamine paper stating *"all genes identified in Bodo… prior to parasitic lifestyles"* — confirm whether it covered AdoMetDC/DHS or only antioxidant genes.
2. The *Bodo* (Jackson 2016) and *Paratrypanosoma* (Flegontov 2013) genome-paper gene tables.

---

## 1. What is SOLID (verified this session)

All of this is from primary data (UniProt, TriTrypDB BLAST, your own runs) + reproducible alignments/trees. True regardless of novelty.

**AdoMetDC active + dead (prozyme-type) pairs:**
- ✓ *Bodo saltans* (FREE-LIVING): active A0A0S4JHX8 / BSAL_30285 (381 aa, full engine) + dead A0A0S4JRG8 / BSAL_30280 (331 aa). Consecutive gene IDs (tandem — coordinates still to confirm, item #1).
- ✓ *Paratrypanosoma confusum* (BASAL trypanosomatid bridge): active PCON_0025430 (390 aa) + dead PCON_0025440 (339 aa). **Syntenic & confirmed**: same scaffold (Pcon_scaffold00006_3), same strand, ~1.7 kb apart; both `is_pseudo=false`.
- ✓ Parasitic trypanosomatids (Leptomonas ×2, Novymonas, Porcisia, Angomonas, Strigomonas, T. brucei, T. cruzi): all have the pair (known).

**Dead-engine signature (definitive residue table, human numbering):**
- ✓ All 16 active copies: E8, E11, S68(pyruvoyl), C82, S229, H243, E247 — intact.
- ✓ All 9 dead copies: lost catalytic **Cys82** (→A/D/S, never C) and **His243** (8/9 → Arg). Shared Glu8→Ser. The His243→Arg is a shared derived "kill-shot" across independent species = one origin.

**Phylogeny (IQ-TREE2, LG+G4, 1000 UFBoot + 1000 SH-aLRT):**
- ✓ AdoMetDC: dead clade **100/100 monophyletic**; branching order **Bodo (basal) → Paratrypanosoma → parasites** = matches species tree. Long-branch artifact fixed with 6 eukaryotic outgroups. → single ancient origin, predating trypanosomatids.
- ✓ DHS (second enzyme): dead-DHS clade (Bodo + T. brucei DHSp + Leishmania DHSL20) **100/100 monophyletic**, Bodo basal. UniProt itself labels the trypanosomatid copies "Inactive deoxyhypusine synthase." → the DHS dead-paralog system is ALSO ancient and in free-living Bodo.

**Two independent polyamine pseudoenzyme systems (AdoMetDC + DHS), both ancient, both in free-living Bodo = an ancestral regulatory strategy, not a parasite invention.** (PRMT is a known third one.)

**Figures delivered to chat:** `prozyme_IQTREE.png`, `DHS_tree.png`, `adometdc_family_tree.png` (earlier), plus `prozyme_phylogeny.*`.

---

## 2. Gap 2 result (new, computational) — "why keep a dead gene?"

Across the 9 dead paralogs, **66 non-catalytic positions are conserved ≥85%**, and at many of them the dead copies are *more* conserved than the active enzymes (e.g., human-pos 93 Pro = 100% in dead vs 25% in active; pos 134 Tyr 100% vs 44%; pos 138 Lys 100% vs 38%; tight conserved blocks ~109–178 and ~255–327).

**Meaning:** a dead gene drifting as junk would not do this. **Selection is actively maintaining a specific surface** → the dead paralog has a function. This is the honest, computable answer to your question. *Limit:* it's conservation evidence, not a structural interface map — calling these "interface residues" needs the complex structure + mutagenesis.

---

## 3. Publication roadmap — ranked gaps (your own angle)

**Gap 1 — The free-living prozyme paradox [sharpest; your question].**
Open: the dead paralog's *function* was only ever tested in *parasites* (it boosts the enzyme ~1000×). Untested in a *free-living* kinetoplastid. Two new outcomes: (a) it still activates → the "weak-enzyme-needs-dead-partner" logic is ancestral; (b) *Bodo*'s enzyme works alone → the dead paralog does something else.
- Compute (done/possible): dead copies conserve a specific surface (Gap 2 ✓); next, map those onto the *T. brucei* prozyme-complex structure to see if they are the activation interface (needs PDB 5TVM).
- **Lab test (for Marc to shop):** express *Bodo* AdoMetDC ± dead paralog; measure activity/activation (standard kinetics). → enzyme-kinetics lab.

**Gap 2 — "What does the dead protein protect?" [evidence in hand].**
The conserved non-catalytic surface (Section 2) is the functional hypothesis.
- **Lab test:** mutate the conserved patch; test loss of activation. → mutagenesis/structural lab.

**Gap 3 — A pseudoenzyme-activator "rule" across the polyamine pathway.**
PRMT + AdoMetDC + DHS all recruited dead-paralog activators. Is there a shared structural logic (all homo-oligomeric enzymes)? Does ODC or spermidine synthase have one? FOLDTRACE comparison across the pathway.

**Gap 4 — Moonshot:** FOLDTRACE scan of a fold family for an unreported lineage-specific active-site loss (fold kept). High ceiling, high risk.

---

## 4. Literature landscape + interlibrary-loan (ILL) list

**Open-access, already have (PMC):** Willert&Phillips 2008 (PMC2562514), 2009 (PMC2730992); Volkov 2016 eLife (PMC5201418); Patel 2018 (PMC6221367); Nguyen 2013 DHS (PMC3663545); Xiao 2013 (PMC3771639); Velez 2013 (PMC3576127); Barker 2009 (PMC2681509); Birkholtz 2011 review (PMC3253432); Marchese 2018 (PMC6027508); Bale&Ealick 2010 review (PMC2847442); Bale 2008 (PMC2646671).

**REQUEST VIA ILL (paywalled / no PMC) — your loan list:**
1. **Willert et al. 2007, PNAS 104:8275** — the original prozyme paper (full text; the kcat/activation numbers).
2. **Opperdoes & Lukeš 2016, J Eukaryot Microbiol 63:657** — you already pulled this ✓ (confirmed: no AdoMetDC mention).
3. **Willert & Phillips 2012, Trends Parasitol 28:66** — review.
4. **Jackson et al. 2016, Curr Biol** — *Bodo saltans* genome (check its gene tables for AdoMetDC/prozyme — novelty check #2).
5. **Flegontov et al. 2013, Curr Biol** — *Paratrypanosoma confusum* genome.
6. **bioRxiv 2019 "Structural Basis of PRMT…"** (likely free on bioRxiv) + **Lorenzon 2022, ACS Infect. Dis.** — the PRMT-prozyme-in-Bodo papers (these are the novelty precedent — read them).
7. The trypanothione/polyamine paper that said *"all genes in Bodo… prior to parasitic lifestyles"* (identify exact ref via its GspS figure; novelty check #1).

**Recent field (2024–2026):** AdoMetDC work is now mostly plants (stress tolerance) and cancer (AMD1); **nothing on kinetoplastid prozyme evolution** — your niche is quiet. One possibly relevant: "Spermine/thermospermine synthases emerged multiple times in eukaryote evolution" (J Biol Chem 2026) — polyamine-enzyme evolution.

---

## 5. Thursday — what to put in front of Marc

1. **The re-run is done** (1MSV/5-site zip) — already delivered.
2. **A real, verified comparative finding:** AdoMetDC + DHS dead-paralog pairs are ancient (free-living *Bodo* + basal *Paratrypanosoma*), single origin, 100/100 trees — with figures.
3. **The honest novelty line** (Section 0) — frame as *extending* the PRMT precedent, not a new discovery.
4. **A concrete lab-testable question** (Gap 1) so he can approach a kinetoplastid enzyme-kinetics lab.
5. **The loan list** (Section 4) — shows you've scoped the literature.

---

## 6. Still open / needs YOU (nothing forced)
- **#1 Bodo tandem coordinates:** look up BSAL_30280 & BSAL_30285 in TriTrypDB (same scaffold/adjacent?) — or BLAST Bodo-active vs *Bodo* and read hit coordinates.
- **#3 full Paratrypanosoma sequences:** download PCON_0025430 & PCON_0025440 full FASTA from TriTrypDB (I used BLAST-partial; consistent, but full is cleaner).
- **Novelty checks #1 & #2** (Section 0).
- **Fact-check fixes:** corrected report/note files are staged; see FIXES_STATUS below / commit when your computer is reachable.

*All analysis files + trees are in this folder. Nothing here was written to Marc's repo or GitHub.*
