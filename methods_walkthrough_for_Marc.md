# Methods walkthrough — the prozyme mapping (point 3) and the diatom lead (point 4)
*For showing Marc. This is the actual step-by-step, with the honest limits at each stage.*

---

# POINT 3 — Do free-living kinetoplastids keep the prozyme activation machinery?

**The question:** In trypanosome parasites, AdoMetDC is switched on by a *catalytically dead paralog* ("prozyme"). Is that whole system a parasite invention, or does it already exist — and still work — in a free-living relative like *Bodo saltans*?

### Step 1 — Ground the biology in primary literature (not memory)
Read the key open-access papers and pulled out exactly what the dead copy does:
- Willert 2007 (PNAS): the prozyme is a dead AdoMetDC copy that heterodimerizes with the active one and boosts activity **~1,200×**.
- Volkov 2016 (eLife): the *mechanism* — the active enzyme is **autoinhibited** by its own N-terminal tail (held by a *cis* proline, **P31**); the prozyme binding flips that proline *cis→trans* and ejects the tail, opening the active site. Mutating two prozyme residues (**M148′, Y152′**) abolishes activation.
- Nguyen 2013 (JBC): the *same* dead-paralog trick independently runs a second enzyme, deoxyhypusine synthase (DHS), at ~3,000× — the active copy keeps the catalytic lysine, the dead copy loses it but keeps substrate-binding residues.

### Step 2 — Establish that Bodo/Paratrypanosoma have the active+dead pairs
- Collected AdoMetDC (and DHS) sequences for *Bodo*, *Paratrypanosoma*, the trypanosome parasites, plus eukaryotic outgroups (human, mouse, fly, yeast, plants).
- Aligned with **MAFFT**; built a tree with **IQ-TREE2** (model LG+G4, 1000 ultrafast bootstraps + 1000 SH-aLRT).
- Result: the **dead/prozyme clade is monophyletic with 100/100 support**, branching **Bodo (basal) → Paratrypanosoma → parasites** — i.e. a *single ancient origin* that tracks the species tree, not independent losses. The DHS dead clade shows the same pattern (100/100, Bodo basal).
- Confirmed every dead copy has lost the catalytic residues (Cys82 and His243 in AdoMetDC numbering), and that the dead copies are ~48% identical to each other (conserved, not decaying).

### Step 3 — The new test: map the *activation interface* residue-by-residue
This is the step beyond "it's dead" — it asks whether the dead copy keeps the surface it uses to *activate* its partner.
1. From Volkov 2016's full text, I listed the exact residues that build the activation mechanism, in T. brucei numbering — on the **active** enzyme (P31 switch; the N-terminal helix residues L8/L10/M11/M13; autoinhibitory F28; D169/H172) and on the **dead** copy (M148′, Y152′, K156′, D171′).
2. **Validated the numbering** against our own alignment before reading anything off — P31 landed on a Pro, C100 on a Cys, the prozyme's M148′/Y152′ on Met/Tyr, etc. (so the positions are trustworthy).
3. For each interface residue, read off what *Bodo* and *Paratrypanosoma* have at the aligned column.

**What came out:**
- **Active enzyme:** Bodo and Paratrypanosoma keep the "activatable" architecture — including **P31, the conformational-switch proline, which humans don't have.** The N-terminal activation residues and D169/H172 are retained; catalytic Cys is intact.
- **Dead copy:** Bodo and Paratrypanosoma retain the **two experimentally-validated activation residues (T. brucei M148′/Y152′ → Bodo M150/Y154)**, plus the salt-bridge/H-bond partners and the prozyme's own putrescine pocket. Only the peripheral edge of the interface varies.
- **DHS (same approach, human DHS numbering):** Bodo's active DHS keeps the catalytic lysine (K329); Bodo's dead DHS has **lost K329 but kept the substrate-binding residues** — the exact signature of the validated T. brucei dead DHS.

**Interpretation:** both halves of the activation machine are present in a *free-living* kinetoplastid → the mechanism most likely **predates the parasites**. This also refines Volkov 2016, who called these features "trypanosomatid-specific" (they never looked at Bodo).

### Honest limits (stated up front)
- Sequence conservation **predicts** function; it doesn't prove it. The clean confirmation is biochemistry (express the Bodo enzyme ± its dead copy, measure activation).
- The Paratrypanosoma sequences are partial (need the full ones from TriTrypDB).
- Novelty: the *general* "this predates parasitism in Bodo" idea is **already published for the PRMT prozyme**; the AdoMetDC+DHS mapping with a dated phylogeny looks unreported but needs two paywalled genome papers to confirm.

### Next steps (ready or defined)
1. **AlphaFold-Multimer** model of the Bodo heterodimer (inputs + analysis checklist already prepared) → test in 3D whether M150/Y154 sit at the interface and the N-terminal helix docks in.
2. Pull full Paratrypanosoma sequences (TriTrypDB).
3. Read the two genome papers to lock novelty.
4. **Bench test** (needs a lab): reconstitute Bodo active ± dead copy, measure activity.

---

# POINT 4 — The diatom "stuck-on" AdoMetDC (an alternative, more independent lead)

**The question:** Human AdoMetDC needs the small molecule putrescine to switch on. Which organisms have rewired that switch — and how?

### Step 1 — Where it came from (a FOLDTRACE survey)
- Took the top ~1,000 structural homologs of human AdoMetDC (Foldseek), kept the 617 that pass the fold gate, and mapped the **putrescine-pocket residues** across all of them.

### Step 2 — The pattern
- A specific pocket rewiring — **Glu178→Asn, Phe111→Gly, Leu13→Thr** — is shared by **6 of 8 diatoms and all 3 brown algae** that passed, and *every one* of the 15 proteins with that exact combination belongs to **Ochrophyta** (the stramenopile group = diatoms + brown algae).
- Plants reach putrescine-independence a *different* way (a built-in Arg/Lys in the pocket). So Ochrophyta look like a **second, independent evolutionary solution** to the same problem.

### Step 3 — Why it's testable (a built-in prediction)
- The human **E178Q** mutant is already known to be **"stuck-on"** — ~4× more active *without* putrescine (Bale 2008). The diatom enzymes carry the natural equivalent (Asn at 178), so the prediction is sharp: **diatom AdoMetDC should be constitutively active / putrescine-independent.**
- Functional logic: diatoms pour polyamines into their **silica cell walls**, so a constitutively-on AdoMetDC fits a high-flux demand.

### Step 4 — Novelty check (done, first-pass)
- PubMed: in all of it, **only one** paper links AdoMetDC to diatoms/brown algae — Michael 2011 (FEBS Lett).
- Checked what it's about: it describes **different** proteins — giant multidomain *fusion* machines that build long-chain polyamines for silica — not the canonical enzyme's regulation.
- Confirmed on UniProt that the WOW-list diatom proteins are the **canonical standalone AdoMetDC** (~320–390 aa, gene "SAMDC"), not those fusions. So the regulation of the *housekeeping* diatom enzyme looks open.

### Why this could be the stronger independent project
- Crisp testable prediction (human E178Q precedent), a **genetically tractable model organism** (*Phaeodactylum tricornutum*), a clear functional story (silica polyamines), and cleaner novelty than the prozyme angle — *provided* Michael 2011 is read in full and carefully distinguished.

### Next steps
1. Pull the exact Ochrophyta "stuck-on" list; build a cross-lineage pocket table + tree to show the convergence is real and single-origin.
2. Model a diatom pocket to show the "on" conformation structurally.
3. Read Michael 2011 (ILL) and draft the novelty framing.
4. **Bench test** (needs a lab): express *Phaeodactylum* AdoMetDC, measure activity ± putrescine.

---

*Both point 3 and point 4 are FOLDTRACE-native (comparative active-/regulatory-site mapping), both end in the same kind of concrete wet-lab experiment, and both have their novelty limits stated honestly.*
