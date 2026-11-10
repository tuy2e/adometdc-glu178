# Does the free-living dead copy keep the activation machinery?
### Residue-by-residue mapping of the prozyme activation interface onto *Bodo saltans* and *Paratrypanosoma confusum*
*Tuy, 4 Nov 2026. 💡 My computational analysis — not Marc's protocol, not for the project folder. Every residue position was taken from a primary source (Volkov et al. 2016 eLife full text for AdoMetDC; Nguyen et al. 2013 JBC + UniProt P49366 for DHS) and mapped onto our own alignments. T. brucei and human numbering were validated against the papers before reading anything off.*

---

## The question
It's one thing to show the Bodo/Pcon dead copies **lost their catalytic residues** (we did that — they're genuinely dead). The sharper question is: **do they keep the machinery that makes a dead copy useful — the surface it uses to grab and switch on its active partner?** If yes, that predicts the whole activation mechanism is **older than the parasites** and was already running in the free-living ancestor. This is the computational test that could become a real result.

I did it for **both** dead-paralog systems.

---

## RESULT 1 — AdoMetDC: both halves of the activation machine are conserved in Bodo

Volkov 2016 solved exactly how the T. brucei dead copy works: the active enzyme is **auto-inhibited** (its own N-terminal tail, held by a *cis* proline **P31**, blocks the active site); the dead copy binds and triggers a conformational switch (P31 flips *cis→trans*, the N-terminal helix docks into the interface) that opens the active site, up to ~1,200× activation. So the machine has two halves — I checked both.

### Half A — does the Bodo/Pcon *active* enzyme still have the "I can be switched on" architecture? **YES.**
| Residue (T. brucei #) | What it does | Human | *T. brucei* | **Bodo** | **Pcon** |
|---|---|---|---|---|---|
| **P31** | the *cis→trans* **switch proline** | absent (T) | P | **P ✓** | **P ✓** |
| L8 / L10 / M11 / M13 | N-terminal activation helix | absent | L/L/M/M | I/L/M/M ✓ | V/L/M/M ✓ |
| F28 (+S27/E29/G30) | autoinhibitory segment | partial | S·F·E·G | S·F·E·G ✓ | S·F·E·G ✓ |
| C100 | catalytic Cys | C | C | **C ✓** | **C ✓** |
| W137 / M146 | connector (mutating these kills activation) | S/F | W/M | W/M ✓ | W/M ✓ |
| D169 / H172 | reach across to the dead copy / putrescine site | N/C | D/H | D/H ✓ | D/H ✓ |

The single most telling one: **P31, the switch proline, is present in Bodo and Paratrypanosoma but absent in humans.** The free-living enzyme carries the exact conformational switch the parasite enzyme uses.

### Half B — does the Bodo/Pcon *dead* copy keep the activator surface? **Mostly YES.**
| Dead-copy residue (T. brucei #) | What it does | *T. brucei* | **Bodo** | **Pcon** |
|---|---|---|---|---|
| **M148′ / Y152′** | **the two residues whose mutation ABOLISHES activation** (Volkov Table 1) | M/Y | **M/Y ✓✓** | **M/Y ✓✓** |
| D171′ | H-bonds the active enzyme's R154 | D | D ✓ | D ✓ |
| K156′ | salt-bridge to the active enzyme's D169 | K | K ✓ | K ✓ |
| E42′/W126′/D195′/D282′ | the dead copy's own putrescine pocket | E/W/D/D | E/W/D/D ✓ | E/W/D/D ✓ |
| R305′/H309′ | peripheral β-sheet edge | R/H | C/E (variable) | R/S (variable) |

The headline: **the two experimentally-proven "activation" residues on the dead copy, M148′ and Y152′, are conserved in free-living Bodo and basal Paratrypanosoma.** Only the peripheral edge of the interface is variable.

**→ Both halves of the Volkov-2016 activation machine are present in a free-living kinetoplastid.**

---

## RESULT 2 — DHS: the Bodo dead copy has the same "dead-but-equipped" signature

For DHS there's **no crystal structure of the complex**, so the interface isn't residue-defined the way AdoMetDC's is. But Nguyen 2013 gave a clear, checkable signature: the active copy keeps the **catalytic lysine (human K329)**; the dead copy **loses it but keeps the substrate-binding residues** (that's how it still contributes to a shared active site). I checked that signature in Bodo.

| Residue (human #) | Role | Human | **Bodo active** | **Bodo DEAD** | *T.b.* DHSc | *T.b.* DHSp (dead) |
|---|---|---|---|---|---|---|
| **K329** | **catalytic nucleophile** | K | **K ✓** | **lost (→Y)** | K | lost (→L) |
| K287 | covalent-intermediate Lys | K | K | **kept ✓** | K | K |
| H288 | spermidine/NAD binding | H | H | **kept ✓** | H | H |
| N106 | NAD/spermidine binding | N | N | **kept ✓** | N | N |
| D238 | NAD binding | D | D | **kept ✓** | D | D |
| W327 | spermidine/covalent | W | W | lost (→H) | W | lost (→N) |

**→ Bodo's DHS dead copy has lost the catalytic lysine K329 but retained the substrate-binding residues (N106, D238, K287, H288) — the exact "dead but still equipped to bind" pattern that the experimentally-validated T. brucei DHSp shows.** It's a bona fide prozyme-type dead copy, not a decaying gene.

---

## What this means (and what it doesn't)

**What it supports.** In *both* systems, the free-living Bodo dead copies aren't just dead — they retain the specific residues their parasite counterparts use for the activation partnership. For AdoMetDC the match is strongest because we could map the actual structure-defined interface (M148′/Y152′ on the dead copy; P31/D169/H172 on the active one). This is concrete evidence that **the dead-paralog activation mechanism predates the trypanosomatid parasites** — i.e. it was likely already working in the free-living ancestor, not invented for the parasitic lifestyle. That directly refines Volkov 2016, who described these features as "trypanosomatid-specific" (they didn't look at Bodo); they're actually kinetoplastid-wide.

This is your "option (a)" from the function review — the weak-enzyme-needs-a-dead-partner logic looks **ancestral** — and it's now backed by a residue-level map, not just a phylogeny.

**What it does NOT prove (the honest limits).**
1. **Sequence conservation predicts function; it doesn't prove it.** The clean test is still the bench experiment: express Bodo AdoMetDC alone vs. with its dead copy and measure activity (the exact assay Willert 2007 / Nguyen 2013 used). My analysis says that experiment is *very likely to work* — which is a strong reason to run it.
2. **P31 being present ≠ the cis-autoinhibition operates identically** — that needs a Bodo structure or a model.
3. **Paratrypanosoma sequences are partial** (reconstructed from BLAST), so the Pcon calls are weaker than the Bodo ones. Pulling the full Pcon sequences from TriTrypDB would firm them up.
4. **DHS has no complex structure**, so that analysis is the catalytic-death + binding-retention signature, not a true interface map.

---

## Why this is useful for you
You now have a **residue-level prediction** to put in front of Marc: *"The Bodo and Paratrypanosoma dead copies retain the experimentally-validated activation-interface residues of the trypanosome prozymes (AdoMetDC: M148′, Y152′, plus the P31 switch and D169/H172 on the active enzyme; DHS: the dead-but-substrate-binding signature). This predicts the prozyme activation mechanism is ancestral to the parasites and should be biochemically reconstitutable from the free-living enzyme pair."* That's a hypothesis with a built-in experiment — exactly the kind of thing that turns an observation into a project.

### Files
- Full result tables + validation: in my working notes (`interface_mapping_RESULTS.md`).
- Residue definitions pulled from the papers: `volkov2016_activation_residues.md`.
- Scripts (reproducible): `map_interface.py`, `map_dhs.py`.
- Alignments used: `pcon_check.aln` (AdoMetDC), `dhs7.aln` (DHS).

### Next concrete steps if you want to push it
1. Download the **full Paratrypanosoma** AdoMetDC/prozyme + DHS sequences from TriTrypDB to upgrade the Pcon calls from "partial" to solid.
2. Build an **AlphaFold / AlphaFold-Multimer model of the Bodo AdoMetDC + dead-copy heterodimer** and check that M148′/Y152′ and the P31/h1 region actually form the predicted interface in 3D (a real figure for a paper).
3. The bench experiment (for Marc to place in a lab): reconstitute Bodo AdoMetDC ± dead copy, measure activation.
