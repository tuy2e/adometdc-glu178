# AlphaFold-Multimer package: does the Bodo dead copy build the activation interface in 3D?
*Tuy, 4 Nov 2026. 💡 My prep — not Marc's protocol. This is everything you need to run the structure prediction and know exactly what to look for. The folding itself needs a GPU service (AlphaFold Server or ColabFold) — it can't run in this workspace — so I've made the inputs submission-ready and written the analysis plan so the output answers a specific question instead of being a pretty picture.*

## Why run this
The sequence mapping showed the Bodo dead copies **keep** the activation-interface residues (AdoMetDC dead copy: M150/Y154 = the Bodo equivalents of T. brucei's mutagenesis-validated M148′/Y152′; active enzyme: the P37 switch proline + N-terminal helix). Sequence conservation *predicts* those residues sit at the heterodimer interface. **A predicted structure tests that directly in 3D** — if the model puts M150/Y154 against the active subunit and docks the active enzyme's N-terminal helix into the interface the way Volkov 2016 saw in T. brucei, that's a real, figure-ready result supporting "the activation mechanism is ancestral/pre-parasitic." It's the strongest thing short of the bench experiment.

## The two jobs (inputs ready)
Files are in `alphafold_inputs/`:
1. **`Bodo_AdoMetDC_heterodimer.fasta`** — Bodo active AdoMetDC (A0A0S4JHX8, 381 aa) + Bodo dead copy (A0A0S4JRG8, 331 aa). **Run this first — it's the headline.**
2. **`Bodo_DHS_heterodimer.fasta`** — Bodo active DHS (390 aa) + Bodo dead DHS (321 aa). Secondary.

Both are full-length, start with Met, no gaps/X. (~710 and ~710 residues total — well within AlphaFold Server limits.)

## How to run it

### Option A — AlphaFold Server (easiest, free, recommended)
1. Go to **alphafoldserver.com** and sign in with a Google account (you do this — I won't log into your account).
2. "New job" → add **two protein entities** → paste the two sequences from the FASTA (active as chain 1, dead as chain 2). Set copies = 1 each.
3. Submit. A heterodimer this size takes ~minutes–tens of minutes. It returns an mmCIF model + a **PAE plot** + confidence scores.
4. Repeat for the DHS pair.
*(If you're already signed in on the desktop browser, tell me and I can fill the sequences in for you with you watching — I just won't sign in or submit on your account without your say-so.)*

### Option B — ColabFold (more control; AF2-Multimer; gives pLDDT + PAE)
Open the "ColabFold: AlphaFold2 using MMseqs2" notebook (the official `ColabFold_batch`), paste each FASTA, set `--templates` on if you want it to see 5TVM, run on the free Colab GPU.

## What to look for — the analysis checklist (this is the point)
Inspect the model against these **specific Bodo residues** (numbers from our alignment, so you can jump straight to them):

**1. Is it even a confident complex?** (do this first)
- **ipTM** (interface pTM): > ~0.6 = a confidently predicted interface; > 0.8 = strong. AF Server reports this.
- **Interface PAE:** in the PAE plot, the off-diagonal blocks (active-vs-dead) should be **low/dark** = the two chains are confidently placed relative to each other. High (light) off-diagonal = don't trust the interface.
- If ipTM is low and PAE off-diagonal is high → the model can't place the complex; report that honestly (still informative, just not a positive).

**2. Does the DEAD copy present its activator residues to the interface?** (the decisive test)
- **Bodo dead M150 and Y154** (= T. brucei M148′/Y152′, the residues whose mutation abolishes activation) — are they **buried at the interface, contacting the active subunit**? If yes → the dead copy is built to activate, in 3D.
- Dead-copy **K158** and **D173** — do they reach across to the active subunit (salt bridge / H-bond)?

**3. Does the ACTIVE subunit adopt the "activated" arrangement?**
- **N-terminal helix (I11/L13/M14/M16)** — docked **into the heterodimer interface** (activated), not folded back blocking its own active site (autoinhibited)?
- **P37** (the switch proline) — positioned at the hinge, consistent with the trans/active state?
- **D175 ↔ dead-copy K158** contact across the interface (mirror of T. brucei D169↔K156′)?
- Active-site **C106** (catalytic Cys) — in an open, competent active-site pocket?

**4. Same architecture as the real parasite complex?**
- Superpose the Bodo model onto **T. brucei 5TVM** (the Volkov heterodimer). Same chain arrangement + the equivalent residues landing in the same interface = strong structural corroboration.

## What each outcome means (honest)
- **Confident interface + M150/Y154 at it + N-helix docked** → strong 3D support that the Bodo free-living enzyme uses the same prozyme activation mechanism; i.e. the mechanism predates the parasites. Figure-ready.
- **Confident interface but residues NOT at it** → interesting — the dead copy may do something other than activation in Bodo (your "option b"). Still publishable.
- **Low-confidence interface** → AF can't resolve it; fall back to the biochemistry. Not a negative result, just inconclusive.

**In all cases: AlphaFold predicts, it doesn't prove.** The confirmation is still expressing Bodo active ± dead copy and measuring activity. But a confident model is a strong figure and a strong reason for Marc to green-light that experiment.

## Loose end: Paratrypanosoma sequences
I tried to upgrade the partial Pcon sequences — they're **not in UniProt or NCBI** (Paratrypanosoma proteins live only in **TriTrypDB**). When you're logged into TriTrypDB, the one-click fix is: open genes **PCON_0025430** (active) and **PCON_0025440** (dead) → Download → protein FASTA. Send me those and I'll drop the "partial" caveat and (optionally) run the Pcon pair through AlphaFold too. Note: the truncation only clipped the extreme N-terminus, so your key Pcon calls (P31, M11/M13, D169/H172, and the dead-copy M148′/Y152′) were already clean.

### Files in this package
- `alphafold_inputs/Bodo_AdoMetDC_heterodimer.fasta`
- `alphafold_inputs/Bodo_DHS_heterodimer.fasta`
- Residue map (what to inspect), above, derived from `pcon_check.aln` + Volkov 2016.
