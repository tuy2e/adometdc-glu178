# The dead copy: what it does, and why it's been kept
### Literature review — catalytically dead enzyme paralogs ("prozymes" / pseudoenzymes), with focus on the AdoMetDC prozyme
*Tuy, 4 Nov 2026. Every paper below was pulled from PubMed and read at the abstract level (several at full text) this session. **Papers are tagged OA (open-access — read now) or PAYWALL (put through interlibrary loan).** Full reading list with links at the end.*

---

## Bottom line up front
Two questions: **what does the dead copy do**, and **why has it survived evolution if it's "dead"?** The literature answers both cleanly:

1. **Function:** the dead copy is not junk — it's a **dedicated allosteric activator**. In trypanosome AdoMetDC it binds the active enzyme and physically flips it out of a self-blocked ("autoinhibited") shape, raising activity **up to ~1,200-fold**. It adds *zero* catalytic residues; all the chemistry is still done by the active partner. The dead copy is a **conformational trigger**, not a second catalyst.

2. **Why conserved:** because a dead paralog that is **maintained under selection is, by definition, doing a job.** Pseudoenzymes make up **~5–10% of most proteomes** and are kept because the enzyme fold is a superb **protein–protein interaction / allostery module**. In trypanosomes specifically there's an extra reason: these parasites **lack transcriptional regulation**, so they lean on protein-level tricks — and dead-paralog activators have evolved **at least three separate times** in their polyamine/translation machinery (AdoMetDC, deoxyhypusine synthase, and PRMT1). That recurrence is the strongest possible evidence these are real, selected functions.

---

## PART 1 — What the dead copy DOES (function)

### The AdoMetDC prozyme — the whole story, step by step

**Discovery (Willert, Fitzpatrick & Phillips 2007, PNAS — OA).** Trypanosomatids uniquely carry *both* a working AdoMetDC *and* a paralog that lost catalytic activity, which they named **"prozyme."** The prozyme binds AdoMetDC in a high-affinity heterodimer and stimulates its activity **1,200-fold**. Their phrase for it: a **"catalytically dead but allosterically active subunit"** — an example of how regulators of multi-subunit enzymes evolve by *gene duplication followed by mutational drift.* This is the paper that named the phenomenon.

**It's conserved across the parasites (Willert & Phillips 2009, Mol Biochem Parasitol — OA).** The same AdoMetDC/prozyme pairing works in *T. cruzi* (110-fold activation) and even in cross-species mix-and-match heterodimers. Mechanistically, **prozyme raises k_cat** (how fast each enzyme turns over) while **putrescine lowers K_m** (how tightly substrate binds) — two independent control knobs on the same enzyme.

**Where it grabs the enzyme (Velez, Brautigam & Phillips 2013, JBC — OA).** The activation works through a **12-residue stretch at the N-terminus of AdoMetDC** (key residues Leu8, Leu10, Met11, Met13). Those residues exist in trypanosomatids but are **absent from human/other eukaryotic AdoMetDCs that have no prozyme** — i.e. the enzyme and its dead partner **co-evolved**. Mutating them doesn't stop the two proteins binding; it stops the *activation* — proving the dead copy works by triggering a shape change, not just by sticking on.

**The actual mechanism (Volkov et al. 2016, eLife — OA) — the most important paper.** X-ray structures showed that **AdoMetDC on its own is autoinhibited**: its own N-terminal sequence sits *inside* its active site and blocks it. When the prozyme docks on, it **displaces that blocking sequence** through a strikingly intricate switch — a **cis→trans proline isomerization**, a **β-sheet reorganization**, and **insertion of the enzyme's N-terminal α-helix into the heterodimer interface** — which opens the active site and turns the enzyme on. Their key evolutionary point: this whole elaborate mechanism was made possible by **one step — the acquisition of a dimerization domain** — which is enough to explain how different regulatory schemes diverged across the AdoMetDC family. *This is the paper that literally answers "what does the dead copy do": it is a conformational key that unlocks an autoinhibited enzyme.*

**The amount of prozyme is itself regulated (Xiao 2013, Mol Microbiol — OA; Patel 2018, PLoS Pathog — OA).** The cell tunes how much prozyme it makes by **translational control**: prozyme protein goes up when AdoMetDC activity drops, and prozyme level is **inversely tied to the pathway product dcAdoMet**. Patel 2018 showed a **dual feedback circuit** — both AdoMetDC protein *and* dcAdoMet suppress prozyme translation, and AdoMetDC does this through an **activity-independent function** (even a catalytically dead AdoMetDC still regulates prozyme). The authors explicitly call this "a regulatory paradigm employing **enzymes and pseudoenzymes that may have broad implications in biology.**"

### Same trick, other enzymes — the dead-paralog strategy is recurrent
This is the part that makes the whole thing look like a *principle*, not a one-off:

- **Deoxyhypusine synthase (DHS) — Nguyen et al. 2013, JBC (OA).** DHS makes the hypusine modification on the translation factor eIF5A. *T. brucei* has **two DHS paralogs** — one badly impaired, one fully dead — and together they form a heterotetramer with a **3,000-fold activity increase**. Both are essential for growth and infection. The authors' key line: *"the **recurrent evolution of paralogous, catalytically dead enzyme-based activating mechanisms** may be a consequence of the unusual gene expression in the parasites, which lack transcriptional regulation."* That sentence is the answer to "why do these keep showing up."
- **PRMT1 (protein arginine methyltransferase) — Kafková et al. 2017, JBC (OA) + Hashimoto et al. 2020, J Mol Biol (OA).** The major trypanosome arginine methyltransferase also runs as an **enzyme–prozyme heterotetramer**; the dead paralog is *required* for the active one to work. The 2020 structure showed the active subunit **can't even fold/dimerize on its own** (a steric clash), so the dead copy acts as a **chaperone** *and* helps bind substrate. The dead copy lost cofactor binding (missing a helix in the Rossmann fold) — that's *why* it's dead — but gained an essential structural role.
- **Plants did it too (Kaltenegger et al. 2021, FEBS Open Bio — OA).** DHS was **duplicated repeatedly in plant lineages** and the copies **recruited to new jobs** (pyrrolizidine-alkaloid biosynthesis) — independent evidence that this enzyme family readily spins off duplicates with new functions.

**So across three different enzymes in the same parasites — AdoMetDC, DHS, PRMT1 — a dead paralog was independently kept and turned into an obligate activator.** That's the single strongest argument that the dead copy has real, selected function.

---

## PART 2 — Why it's been conserved for so long

### The general principle (pseudoenzyme field)
- **Ribeiro et al. 2019, Science Signaling (PAYWALL — ILL this)** — the field-defining review. **5–10% of enzyme families across bacteria, archaea and eukaryotes are pseudoenzymes**, with "notable expansions and contractions associated with specific signaling niches" (i.e. they're gained and lost *for reasons*). Documented functions: **allosterically activate canonical enzymes, act as scaffolds that control assembly and localization of complexes, serve as molecular switches, and regulate networks by sequestering substrates or enzymes.** Some even turn out to retain cryptic activity. They're druggable.
- **Murphy, Mace & Eyers 2017, Curr Opin Struct Biol (PAYWALL — ILL this)** — "Live and let die." Same core list from a structural angle: pseudoenzymes as **allosteric regulators, molecular switches/integrators, assembly hubs, and competitors for substrate/holoenzyme assembly.**
- **Mace & Murphy 2021, JBC (OA)** — "There's more to death than life." Focused on pseudokinases (the best-studied class) but the logic is general: a dead enzyme is kept because the fold is a **versatile protein-interaction module** and can **switch between states** — exactly the properties you need for an allosteric regulator. Catalytic activity is "dispensable for biological function."

**The conservation logic in one sentence:** a gene with no function decays and is lost (pseudogene); a *dead protein that stays intact, expressed, and under purifying selection across millions of years is being kept because it does something* — and for pseudoenzymes that "something" is almost always a non-catalytic job (allostery, scaffolding, switching, sequestration).

### Why trypanosomes in particular keep doing this
Trypanosomatids **transcribe their genes polycistronically and don't regulate transcription of individual genes** — so they can't dial enzymes up and down the way we do. They compensate with **post-transcriptional and post-translational** control. A dead-paralog activator is exactly that kind of control: you regulate an enzyme by regulating how much of its dead partner you make (which, as Xiao 2013 and Patel 2018 showed, is itself feedback-tuned by the pathway's product). This is why the strategy **recurs** in these organisms (Nguyen 2013; Kafková 2017; Phillips 2018 review, OA).

### How this maps onto your own data
Your Bodo/Paratrypanosoma dead copies being **~48% identical to each other and conserved (not decaying)** is precisely the molecular signature the pseudoenzyme field uses to call something "maintained under selection → functional." The published principle predicts your conserved dead copies are doing a job; your data are consistent with that. That's a defensible, literature-backed statement.

---

## PART 3 — The open question that's actually yours
Here's where the literature *stops*, which is where your project can stand:

Every functional study above — the ~1,200× (AdoMetDC), 3,000× (DHS), the PRMT chaperone role — was done in **parasites**, where there's an obvious reason to have evolved tight polyamine control. The published "ancient origin / predates parasitism" idea has so far been shown mainly for the **PRMT prozyme** (including in free-living *Bodo*, established in prior work — this is why your angle is an *extension*, not a first). **What has not been pinned down** is the **function** of the AdoMetDC (and DHS) dead copy in a **free-living** organism that has no parasitic lifestyle to regulate. Two clean, testable possibilities:
- (a) Bodo's dead copy still activates its enzyme → the "weak-enzyme-needs-a-dead-partner" logic is **ancestral**, older than parasitism; or
- (b) Bodo's enzyme works fine alone → the dead copy is doing **something else** there, and was later *co-opted* for activation in the parasite lineage.

Either answer is publishable, and the experiment is concrete (express Bodo AdoMetDC ± its dead paralog, measure activity — the exact assay Willert 2007 and Nguyen 2013 used). The computational precursor (does Bodo's dead copy keep the activation interface mapped in Volkov 2016?) you can do before any bench work.

---

## PART 4 — Reading list

### ✅ Open access — read these now, no ILL needed
| Paper | What it gives you | Link |
|---|---|---|
| Willert, Fitzpatrick & Phillips 2007, *PNAS* 104:8275 | **Discovery** of the prozyme; 1,200× activation | PMC1895940 |
| Willert & Phillips 2009, *Mol Biochem Parasitol* 168:1 | Conserved across trypanosomatids; k_cat vs K_m | PMC2730992 |
| Velez, Brautigam & Phillips 2013, *JBC* 288:5232 | The N-terminal activation residues; co-evolution | PMC3576127 |
| **Volkov et al. 2016, *eLife* 5:e20198** | **THE mechanism — autoinhibition relief** | PMC5201418 |
| Xiao et al. 2013, *Mol Microbiol* 88:846 | Translational control of prozyme; dcAdoMet feedback | PMC3771639 |
| Patel et al. 2018, *PLoS Pathog* 14:e1007404 | Dual feedback circuit; activity-independent function | PMC6221367 |
| Nguyen et al. 2013, *JBC* 288:15256 | **DHS dead paralog** (3,000×); "recurrent evolution" line | PMC3663545 |
| Kafková et al. 2017, *JBC* 292:2089 | **PRMT1 enzyme–prozyme** complex | PMC5313084 |
| Hashimoto et al. 2020, *J Mol Biol* 432:410 | PRMT prozyme **structure**; chaperone role | PMC6995776 |
| Phillips 2018, *JBC* 293:18746 (minireview) | Polyamines in protozoan pathogens — good overview | PMC6290135 |
| Mace & Murphy 2021, *JBC* 296:100705 | Noncatalytic functions of (pseudo)kinases — general principles | PMC8141879 |
| Kaltenegger et al. 2021, *FEBS Open Bio* 11:10 | Plant DHS duplication / neofunctionalization | PMC7780104 |
| Aroonsri et al. 2019, *PeerJ* 7:e6713 | *Plasmodium* DHS as a drug target (pathway importance) | PMC6475138 |

### 🔒 Paywalled — put these through interlibrary loan
| Paper | Why you want it | ID for ILL |
|---|---|---|
| **Ribeiro et al. 2019, *Science Signaling* 12:eaat9797** — "Emerging concepts in pseudoenzyme classification, evolution, and signaling" | **The** field review; your "why conserved / 5–10% of proteomes / expansions & contractions" citation | PMID 31409758; doi 10.1126/scisignal.aat9797 |
| **Murphy, Mace & Eyers 2017, *Curr Opin Struct Biol* 47:95** — "Live and let die: insights into pseudoenzyme mechanisms from structure" | The structural catalogue of pseudoenzyme functions (allostery/switch/scaffold/sequestration) | PMID 28787627; doi 10.1016/j.sbi.2017.07.004 |
| **Willert & Phillips 2012, *Trends Parasitol* 28:66** — review | Phillips-lab review of the prozyme regulatory strategy | PMID 22192816; doi 10.1016/j.pt.2011.11.001 |
| **Kim & Barrett 2013, *Mol Microbiol* 88:841** — commentary | Commentary on the dcAdoMet/prozyme feedback story | PMID 23668674; doi 10.1111/mmi.12243 |
| Carrizo et al. 2026, *Nature* 651:763 — "The transition from monocyte to tissue-resident macrophage requires DHPS" | Shows the *active* hypusine/DHS pathway matters in mammals (context, not dead-copy) | PMID 41565804; doi 10.1038/s41586-025-09972-2 |

*(Check your library first — a few "paywalled" ones, especially Ribeiro 2019, are sometimes free via Europe PMC or the author's site. If any OA link above fails to open, tell me and I'll find an alternate.)*

---

## One honest caveat (so you're not blindsided)
The **general idea** that prozyme/pseudoenzyme regulation is ancient and present in free-living kinetoplastids is already in the literature for the **PRMT prozyme** (including *Bodo*). So your contribution is the **AdoMetDC + DHS dead-paralog mapping into free-living/basal lineages with a dated phylogeny**, plus the **free-living function question** above — an extension of a known theme into new systems, which is a legitimate, publishable unit of work (and the kind of "not previously explored" result you're aiming for). Confirming the two paywalled genome papers (Jackson 2016, Flegontov 2013) is still the step that locks down exactly how new the mapping is.
