# Literature check: has anyone compared AdoMetDC putrescine-site chemistry across species?

**Done 10 Nov 2026, Stage 1 step 1 of the pocket-chemistry survey**
(plan: `1_reports/HANDOFF_pocket_chemistry_survey.md`).

**Answer: no.** Nothing found computes or compares pocket chemistry — electrostatics, pKa,
cavity, site charge — across AdoMetDC from different species.

## How it was searched

Europe PMC, title+abstract fields only (`TITLE_ABS:`), via the desktop browser. The cloud
workspace cannot reach Europe PMC, PubMed, NCBI, UniProt or RCSB (the proxy refuses those
hosts), and web search is switched off for the organisation, so every search in this file was
run by hand in the browser. Searches are reproducible by pasting the query into
<https://europepmc.org>.

| # | Query (title/abstract) | Hits | What they were |
|---|---|---|---|
| 1 | ("adenosylmethionine decarboxylase") AND (electrostatic OR "Poisson-Boltzmann" OR pKa OR "molecular dynamics" OR docking OR "homology model") | 21 | inhibitor-design docking (Ai 2025, Ai 2016, Brooks 2007, Grover 2012, Zhao 2025), single-structure papers, one *Leishmania* AdoMetDC-like protein (Singh 2013), *Pf*ODC modelling (Birkholtz 2003). **No cross-species pocket comparison.** |
| 2 | (AdoMetDC OR SAMDC OR "adenosylmethionine decarboxylase") AND putrescine AND (comparative OR evolution OR phylogenetic OR species OR survey) | 40 | plant gene-family expression surveys (pepper, cucumber, peanut, bentgrass), trypanosomatid drug reviews, polyamine-metabolism reviews, Li 2023 on bacterial/archaeal neofunctionalisation. **No structural survey.** |
| 3 | (AdoMetDC OR "adenosylmethionine decarboxylase") AND ("binding site" OR pocket OR cavity OR "regulatory site"), sorted by citations | 16 | exactly the papers already read in full — Bale 2008, Bennett 2002, Clyne 2002, Beswick 2006, Wells 2006 — plus the human crystal papers (Ekstrom 1999, Tolbert 2001) and Makhoba 2016. |

## The nearest existing work (and why it isn't this)

- **Bale S. & Ealick S.E. 2009, *Amino Acids* 38:451–460, PMID 19997761** (open access).
  "Structural biology of S-adenosylmethionine decarboxylase." A review comparing the known
  AdoMetDC structures — human, potato, *Thermotoga*, and the class-1 bacterial enzymes.
  **This is the closest prior art and the paper a reviewer will raise.** It is descriptive:
  it compares folds, active sites and the putrescine site residue by residue, but computes no
  electrostatics, no pKa and no cavity volumes, and covers only the handful of species that had
  crystal structures in 2009 — no stramenopile, cryptophyte or other protist.
  **Not yet read in full. Worth reading before writing up.**
- **Toms A.V. et al. 2004, *J Biol Chem* 279:33837, PMID 15150268.** *Thermotoga maritima*
  AdoMetDC structure, framed as "evolutionary links". Structural evolution of the fold, not
  regulation; prokaryotic.
- **Li B. et al. 2023, *J Biol Chem* 299:105005, PMID 37399976.** Neofunctionalisation of
  AdoMetDC into ornithine/arginine decarboxylases in bacteria and archaea. Sequence and
  biochemistry; a different question, and already noted in the log as Li-group territory.

## What this means for the survey

1. **The gap is real.** No computed pocket chemistry across AdoMetDC species exists, so the
   survey is not duplicating anything.
2. **Novelty should still be stated modestly.** Bale & Ealick 2009 already made the qualitative
   comparison for the enzymes with crystals; our contribution is computing it, and extending it
   to lineages with no structural work at all. The claim is "first systematic computed survey",
   not "first comparison".
3. **Read Bale & Ealick 2009 before writing up** — it is open access, and it is the one paper
   that could pre-empt a framing claim.

## Structures found while checking (RCSB, same session)

Worth recording because the handoff assumed no trypanosome crystal existed:

- **5TVO** — *Trypanosoma brucei* AdoMetDC monomer, 1.48 Å, Volkov et al. 2016 *eLife*.
- **5TVF** — *T. brucei* AdoMetDC/prozyme heterodimer with CGP 40215, 2.42 Å, same paper.
  → A real crystal of a *T. cruzi*-type (activated, no Glu178) enzyme, instead of relying on
  the truncated AlphaFold model that limited the 10 Nov docking round.
- Human crystals with putrescine bound beyond 1MSV: **1JL0** (H243A proenzyme, 1.5 Å),
  **1I7B** (processed, 1.9 Å), **3DZ4** and the rest of the Bale 2009 inhibitor series.
  → Free replicates for estimating the method's own noise (gate criterion G2).
