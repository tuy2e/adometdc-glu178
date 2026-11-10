# AdoMetDC — everything we know (a deep dive)
### 💡 Tuy's research compilation, 4 Nov 2026 — NOT Marc's protocol, NOT for the project folder

*You asked: what do we know about AdoMetDC across every structure, what's exceptional, what's exceptional about the residues, what's new, does it have metals, how do animals use it — all of it. Here it is. **Every fact is checked against a primary source this session** (RCSB PDB, UniProt P17707, PubMed/PMC, and the Bale & Ealick 2010 structural review read in full). Where a number comes only from a review and I haven't opened the original kinetics paper, I say so. Sources are listed at the end.*

---

## 0. The one-paragraph version
AdoMetDC (S-adenosylmethionine decarboxylase; human gene **AMD1**, UniProt **P17707**, EC 4.1.1.50) makes the aminopropyl donor for spermidine and spermine — it's the committed, rate-limiting gate of polyamine synthesis. It is one of the small club of **pyruvoyl-dependent** enzymes: instead of using vitamin-B6 (PLP) like most amino-acid decarboxylases, it **builds its own cofactor out of its own backbone** by cutting itself in two. It recognises its substrate's positive sulfonium charge by **sandwiching it between two phenylalanine rings** (a cation-π interaction), it is switched on allosterically by **putrescine**, and in parasites it is switched on by a **dead copy of itself** (the "prozyme"). Most AdoMetDCs use **no metal** (humans included); only the bacterial/archaeal class-1 enzymes need Mg²⁺. In animals it's tied to stem-cell self-renewal and is emerging as a **cancer metabolic vulnerability**, especially in the prostate.

---

## 1. What it does (verified: UniProt P17707)
- **Reaction:** S-adenosyl-L-methionine (SAM) + H⁺ → S-adenosyl 3-(methylsulfanyl)propylamine (**decarboxylated SAM, "dcSAM"**) + CO₂.
- dcSAM is the **aminopropyl donor** that spermidine synthase and spermine synthase use to build the higher polyamines from putrescine. This is the **only** job of the enzyme and it's a one-step committed pathway (UniProt: "step 1/1").
- **Function line (verbatim from UniProt):** "Essential for biosynthesis of the polyamines spermidine and spermine. Promotes maintenance and self-renewal of embryonic stem cells, by maintaining spermine levels." → so beyond housekeeping, it has a developmental/stem-cell role.

---

## 2. What is genuinely *exceptional* about this enzyme

**(a) It makes its own cofactor by cutting itself in half — "non-hydrolytic serinolysis."**
AdoMetDC is born as a single inactive proenzyme. It then **self-cleaves** at Glu67–Ser68 with no help from any protease. The trick (UniProt, verbatim mechanism): the side-chain –OH of Ser68 donates its oxygen to cap the new C-terminus of one fragment, while the rest of Ser68 undergoes oxidative deamination to spit out ammonia and leave a **pyruvoyl group** behind. That pyruvoyl group *is* the catalytic cofactor. Only a handful of enzymes in all of biology do pyruvoyl self-processing — this is textbook-exceptional.
- Result: a **β chain (residues 1–67)** and an **α chain (residues 68–334)**, with the pyruvoyl cofactor at the **N-terminus of the α chain**. The mature human enzyme is a **heterotetramer, (αβ)₂**.
- The processing has been caught in the act crystallographically: the **H243A mutant (PDB 1JL0, 1.5 Å)** is frozen at the **ester intermediate** of the reaction.

**(b) It reads a +charge with an aromatic sandwich (cation-π).**
The substrate SAM carries a positively charged **sulfonium** (S⁺). AdoMetDC grips it between two phenylalanine rings, **Phe7 and Phe223** — a cation-π interaction. Bale & Ealick's review puts this at **~4.5 kcal/mol** of binding energy. Both residues were mutated and crystallised (Tolbert et al. 2001; F223A = PDB 3DZ3), confirming they matter. Recognising a sulfonium this way is an unusual and elegant recognition strategy.

**(c) It's turned on by putrescine — allostery from 16–20 Å away.**
Human AdoMetDC has a buried acidic pocket that binds **putrescine**, and putrescine stimulates **both** the self-processing step **and** the decarboxylation step (UniProt "activity regulation"). That's a feedback loop: putrescine (made upstream by ODC) accelerates the enzyme that feeds spermidine/spermine synthesis. The pocket residues were dissected by making mutants and solving them with and without putrescine (PDB 3EP3–3EPB).

**(d) It's turned on by a dead copy of itself — the "prozyme" (parasites).**
In trypanosomes, the active AdoMetDC pairs with a **catalytically dead paralog ("prozyme")** that boosts activity enormously. The mechanism was solved in 2016: **"Relief of autoinhibition by conformational switch explains enzyme activation by a catalytically dead paralog"** (Volkov et al., eLife; PDB 5TVF/5TVM/5TVO). *(This is directly relevant to your own project question about why the dead paralog exists and what it does — the answer in T. brucei is that it flips the active enzyme out of an autoinhibited conformation.)*

**(e) It is an extraordinarily faithful catalyst — and Cys82 is why.**
Reviews report the wild-type enzyme runs **~15,000 catalytic cycles** before a rare transamination side-reaction kills it, and that the **C82A mutant collapses to ~143** cycles — i.e. Cys82 enforces correct re-protonation and protects the enzyme from self-inactivation. *(These exact numbers are from the Bale & Ealick 2010 review; the primary kinetics is Xiong, Stanley & Pegg 1999, Biochemistry — worth pulling the original if you want to quote the numbers in writing.)*

---

## 3. The exceptional residues (the cast of characters, human numbering)
*All verified via UniProt mutagenesis records and/or dedicated PDB structures.*

| Residue | Role | How we know |
|---|---|---|
| **Ser68** | Becomes the **pyruvoyl cofactor** on self-cleavage | Mutagenesis S68A blocks processing (Stanley & Pegg); PDB **1MSV = S68A** non-processing proenzyme (**this is Marc's reference structure**) |
| **Cys82** | **Catalytic acid** — reprotonates Cα, enforces fidelity (~15,000 vs ~143 turnovers) | C82A kinetics (Xiong et al. 1999); inhibited by iodoacetic acid (UniProt) |
| **His243** | **Catalytic base**; also drives self-processing | H243A traps the ester intermediate (PDB **1JL0**); mutagenesis (Xiong & Pegg 1999) |
| **Ser229** | Assists self-processing | Mutagenesis (Xiong & Pegg 1999) |
| **Phe7 + Phe223** | **cation-π sandwich** on the substrate sulfonium (~4.5 kcal/mol) | Mutagenesis (Tolbert 2001); F223A = PDB 3DZ3 |
| **Glu247** | Orients the adenine/substrate | Bale & Ealick review (substrate-binding set) |
| **Asp174, Glu178, Glu256** | **Putrescine-activation pocket** | Mutant structures ± putrescine: PDB 3EP3–3EPB; mutagenesis (Stanley & Pegg 1991) |

Plant twist worth knowing: the **potato** enzyme doesn't use putrescine at all — it has **built-in arginines (Arg18, Arg114)** that structurally mimic the putrescine charge, so it's constitutively "activated" (Bale & Ealick review). A nice example of the same job solved two different ways.

---

## 4. Does it have metals? — the honest answer
**Mostly no.** I checked every one of the **40 experimental PDB entries** for bound metal:
- **Only ONE** has catalytic-looking metal: **2III** (*Aquifex aeolicus*, a class-1b bacterial enzyme) — bound **Mg²⁺ and Ca²⁺**.
- **Every human, potato, Thermotoga and Trypanosoma structure is metal-free** — their only extra atoms are buffer/cryo components (tris, glycerol, glycols, Cl⁻/Na⁺) or functional ligands (putrescine, inhibitors, the pyruvoyl group).
- The textbook statement that **bacterial/archaeal class-1a enzymes (the E. coli type) require Mg²⁺** rests on **biochemistry, not a crystal structure** — there is **no E. coli AdoMetDC (SpeD) structure in the PDB at all** (I confirmed this; it also fixes an error that was in one of your research notes).

**Bottom line:** metal dependence is a **bacterial/archaeal class-1 feature**; eukaryotic AdoMetDC (humans, plants) is metal-free and uses putrescine or built-in arginines instead.

---

## 5. Every known AdoMetDC structure (all 40, verified from RCSB, newest first)

### The newest work — "what's new" structurally (2025)
Three entries released **17 Sep 2025**, all **human**, a **temperature series**:
- **9P1H** — 100 K, 1.81 Å · **9P7Q** — 273 K, 2.21 Å · **9PBB** — 293 K, 2.17 Å (all + putrescine)
- Paper: **Patel, Bonzon, Bakht, Fagbohun, Clinger, "Multi-Temperature Crystallography of S-Adenosylmethionine Decarboxylase Observes Dynamic Loop Motions," *Biomolecules* 2025** (PMID 41008581). Multi-temperature crystallography is done to see **flexible loops / motions** that are frozen out at cryogenic 100 K — i.e. how the enzyme actually *moves* near body temperature. This is the frontier of AdoMetDC structural work right now, and it's about **dynamics**, not a new drug.

### The parasite prozyme heterodimer (the drug-target structures)
*These post-date the 2010 review, which said no such structure existed.*
- **6BM7** (2018, 2.98 Å) — T. brucei AdoMetDC/**prozyme** + pyrimidineamine inhibitor UTSAM568 — *Volkov et al., J Med Chem 2018* (PMID 29271204).
- **5TVF / 5TVM** (2017, ~2.4 Å) — T. brucei AdoMetDC/**prozyme** heterodimer (+inhibitor CGP 40215 in 5TVF) — *Volkov et al., eLife 2016* (PMID 27977001).
- **5TVO** (2016, 1.48 Å) — T. brucei active subunit (Δ26) monomer.

### Prokaryotic / archaeal (class 1)
- *Thermotoga maritima* (class 1b): **1VR7** (1.2 Å, proenzyme), **1TLU** (1.55 Å), **1TMI** (S63A non-processing mutant), **3IWB** (processed), **3IWC** (+SAM methyl ester), **3IWD** (+dimethyl-thioadenosine).
- *Aquifex aeolicus* (class 1b): **2III** (2.3 Å) — **the one with Mg²⁺ + Ca²⁺**.

### Plant (class 2a)
- **1MHM** — potato (*Solanum tuberosum*), 2.3 Å. Uses built-in Arg; no putrescine, no metal.

### Human (class 2b-I) — the large mutant/inhibitor collection
- **Foundational:** 1JEN (1999, first human), 1I7B (+SAM methyl ester), 1I72/1I79/1I7M (inhibitor complexes), **1I7C (+MGBG / mitoguazone — classic anticancer inhibitor)**, 1I7M (+ an amidinoindanone of the SAM486A class).
- **Mechanism:** **1MSV (S68A proenzyme — Marc's reference)**, **1JL0 (H243A ester intermediate)**.
- **Putrescine-site dissection (Bale et al. 2008):** 3EP3/3EP6 (D174N), 3EP4/3EP7/3EPB (E256Q), 3EP5/3EP8/3EPA (E178Q), 3EP9 (apo, no putrescine).
- **Sulfonium / cation-π + 8-substituted-adenosine inhibitors (Ealick lab 2009):** 3DZ3 (**F223A**), 3DZ2/3DZ4/3DZ5/3DZ6/3DZ7, 3H0V/3H0W.

*(Full machine-readable table with organism/method/resolution/ligands for all 40 is in my working notes — happy to hand it over if useful.)*

---

## 6. How animals / humans use it — and how it's controlled

**Where & what for.** AMD1 sits on **human chromosome 6**, 334-aa proenzyme, cytosolic. It's broadly expressed (UniProt/Bgee: >200 tissues) but **tissue-enhanced in prostate** (Human Protein Atlas). Its spermidine/spermine output is needed for cell growth and proliferation generally, and specifically for **embryonic stem-cell self-renewal via spermine** (UniProt).

**It's one of the most tightly regulated enzymes in the cell — at several levels:**
1. **Allosteric activation by putrescine** (above) — a built-in feed-forward from the upstream pathway.
2. **Self-processing** as a committed activation step (the proenzyme is inert until it cleaves itself).
3. **Translational control.** AdoMetDC is a classic example of **polyamine-responsive translational regulation** via its mRNA — there's a deep literature on this (52 papers in PubMed on AdoMetDC translational control). **Brand-new (2026):** Maldosevic et al., **"Mechanism of ribosome stalling by the AMD1 C-terminal tail arrest peptide," *Science Advances* 2026** (PMID 41894501) — a newly worked-out **nascent-peptide / ribosome-stalling** mechanism controlling AMD1 expression. This is genuinely current, high-profile mechanistic work on the human enzyme.
   *(Note: ODC — the other polyamine gatekeeper — is destroyed by "antizyme"; AdoMetDC's control is dominated by putrescine + translational mechanisms, not the antizyme-degradation route, so don't conflate the two.)*

**Disease & drug angle (modern).**
- **Cancer vulnerability — prostate especially (2026):** Kumar et al., **"Androgen Receptor Drives Polyamine Synthesis, Creating a Vulnerability for Prostate Cancer," *Cancer Res* 2026** (PMID 41270152), with a companion commentary (PMID 41766411). Fits the prostate-enhanced expression.
- **AMD1 across disease (2026):** ferroptosis in rheumatoid-arthritis synoviocytes (PMID 42685934), kidney tubular repair/senescence (PMID 42289383), age-related muscle decline (PMID 42602328), CD8⁺ T-cell antitumor function (PMID 42503479). AMD1 shows up as **depleted/essential in many CRISPR screens** (UniProt/BioGRID-ORCS: 275 hits across 1,156 screens) — a dependency signal.
- **Drugs that target it:** **MGBG / mitoguazone** (DB12967), **SAM486A / sardomozide / CGP48664** (DB20551) — the latter still being tested in 2026 combinations (e.g. with the KRAS-G12C drug adagrasib in lung cancer, PMID 41715254). These are the amidinohydrazone-class inhibitors you see in the 1I7C/1I7M structures.
- **Polyamine-imbalance disease context (2023):** DFMO rebalancing polyamines in Snyder-Robinson syndrome (PMID 37702369) — shows the pathway (not AdoMetDC alone) is clinically actionable.

---

## 7. How it differs across life (the five-class map, Bale & Ealick 2010)
| Class | Who | Oligomer | Activator / metal |
|---|---|---|---|
| **1a** | Gram-neg bacteria, some archaea (**E. coli** SpeD) | tetramer | **Mg²⁺ (metal-dependent)** — biochemistry only, no PDB structure |
| **1b** | Gram-pos bacteria, some archaea (*T. maritima, A. aeolicus*) | dimer | none known (2III has adventitious/structural Mg+Ca) |
| **2a** | plants (**potato**) | — | built-in **Arg18/Arg114** (no putrescine, no metal) |
| **2b-I** | **animals / humans** | (αβ)₂ | **putrescine** (no metal) |
| **2b-II** | trypanosomatid parasites | heterodimer | **prozyme** (dead paralog) + putrescine |

The eukaryotic fold is a **four-layer αββα sandwich with internal pseudo-two-fold symmetry** — it arose by **gene duplication + fusion** of an ancestral single-domain (prokaryotic-type) unit. That duplication-and-repurposing theme is exactly the kind of thing your prozyme/pseudoenzyme project is about.

---

## 8. Threads that connect straight to *your* project
- **1MSV** (Marc's new reference) = the human **S68A** non-processing proenzyme — i.e. the enzyme deliberately frozen *before* it makes its pyruvoyl cofactor. Good to be able to say why that mutant exists.
- The **prozyme** is the live research frontier here: the 2016 eLife mechanism ("relief of autoinhibition by a conformational switch") is the published answer to "what does the dead paralog *do*" **in parasites** — which is exactly the contrast that makes your free-living-*Bodo* question interesting (does the ancient dead copy do the same job where there's no parasitic lifestyle to regulate?).
- The 2018 Patel/Phillips paper (PMID 30365568) — "a dual regulatory circuit… controls expression of the paralogous activator prozyme in T. brucei" — is worth reading for how the active enzyme and its product jointly tune the prozyme. Relevant to any "what maintains the dead copy" argument.

---

## 9. Sources (all opened/verified this session)
- **UniProt P17707** (DCAM_HUMAN) — function, reaction, pyruvoyl cofactor, serinolysis mechanism, α/β chains, (αβ)₂, putrescine regulation, residue mutagenesis, PDB list, drugs, prostate-enhanced, CRISPR-ORCS. https://www.uniprot.org/uniprotkb/P17707
- **RCSB PDB** — all 40 EC 4.1.1.50 entries (organism/method/resolution/ligands/release dates) + primary citations, via search + data GraphQL APIs.
- **Bale & Ealick 2010**, "Structural biology of S-adenosylmethionine decarboxylase," *Amino Acids* 38:451–460 (PMC2847442) — classification, mechanism, turnover numbers, cation-π, class table (read in full).
- **PubMed/PMC** (E-utilities) — recent-literature sweep (103 hits 2023–2026; 80 AdoMetDC-relevant) and the specific papers cited above:
  - Patel et al. 2025, *Biomolecules* — multi-temperature crystallography (PMID 41008581).
  - Volkov et al. 2016, *eLife* — prozyme autoinhibition switch (PMID 27977001).
  - Volkov et al. 2018, *J Med Chem* — species-selective inhibitors (PMID 29271204).
  - Patel et al. 2018, *PLoS Pathog* — dual regulatory circuit / prozyme (PMID 30365568).
  - Maldosevic et al. 2026, *Sci Adv* — AMD1 arrest-peptide ribosome stalling (PMID 41894501).
  - Kumar et al. 2026, *Cancer Res* — AR drives polyamine synthesis in prostate cancer (PMID 41270152).
  - plus AMD1-disease papers PMIDs 42685934, 42289383, 42602328, 42503479; SAM486A combo PMID 41715254.

### Flagged as NOT re-verified to the original this session (safe to use, but confirm before quoting a number)
- The exact **~15,000 vs ~143 turnover** numbers (from the 2010 review; primary = Xiong, Stanley & Pegg 1999, *Biochemistry*, PMID 10029540).
- The **~4.5 kcal/mol** cation-π figure (from the 2010 review).
