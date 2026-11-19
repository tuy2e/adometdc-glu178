# Handoff: AdoMetDC pocket-chemistry survey (starts after 10 Nov 2026)

This is for a new chat. It is written so Claude can start without the old conversation.

## The idea in one paragraph

The evolution result stands: every stramenopile lost Glu178, once, at the origin of the group; it is also missing in a few other lineages whose enzymes cluster with them against the species tree. Whether those enzymes are still switched on by putrescine can't be settled without a lab, and Tui can't do lab work.

So the new question is chemical and answerable by computer: **can pocket chemistry computed from structures tell apart the known ways AdoMetDC handles putrescine, and where do the Glu-less lineages fall?** The known ways:

- **Human-type:** activated, Glu178 kept.
- **Plant-type:** pocket filled by built-in Arg/His, always on.
- ***T. cruzi*-type:** activated through a shifted site, no Glu178.
- ***Plasmodium*-type:** not activated, Lys at 174.

The model is Begg, Mason & Zimmer 2026, *J. Chem. Inf. Model.* 66:10266 (Luke and Maddie's paper). It is a systematic structure survey, computational only, validated against known data, with honest limits. Their code: github.com/headsparkle/gfp-barrel-geometry.

## Template: Marc's HEPN paper (follow its methods)

Zimmer M. 2026, "Recurrent catalytic-site change and electrostatic remodeling across the HEPN ribonuclease fold", *J. Mol. Biol.* 439:170055, doi:10.1016/j.jmb.2026.170055. FOLDTRACE v0.2.0, Zenodo doi:10.5281/zenodo.22081699. Tui has the PDF; it is not in the repo (copyright).

His question was what happens to the old active site when evolution removes the catalytic His: does it decay, or is it physicochemically remodelled? **Ours is the same question for a regulatory site:** when Glu178 is lost, does the putrescine pocket decay, fill (plant-like), or get electrostatically remodelled?

| His method | Our version |
|---|---|
| FOLDTRACE fold search + 3D mapping of the chemistry-defining residues (Arg97/His102) | The existing `stuy__adometdc` run (1,000 models); functional centre = 1MSV putrescine centroid; residues 13/15/111/174/176/178/256/285 |
| Pocket retention: fpocket 4.0. The canonical cleft is the fpocket call nearest the functional centre. The distance cutoff is the 95th percentile from intact controls (leave-one-control-out); volume ≥ 80 Å³; lining pLDDT ≥ 70 | The same, with Glu178-retaining enzymes as controls |
| Electrostatics: PDB2PQR 3.6.1 (AMBER, pH 7.0) → APBS 3.4.1 (97³ grid), interpolated at the fpocket alpha-sphere points, reported as the median per pocket | The same |
| Classes: active-like / inverted / decay-like / remodelled-unresolved, using control median ± 3.5 robust SD plus cleft, Arg/Lys and pLDDT gates | human-like / filled (plant-type) / remodelled-shifted (*T. cruzi*-type) / decayed / unresolved |
| Statistics: 50% and 70% identity cluster reduction with bootstrap CIs; one-to-one patristic matching with an exact sign-flip test; spread across 1,000 UFBoot topologies | The same (Glu-less vs nearest Glu-retaining enzyme) |
| Phylogeny: Fitch minimum H/non-H changes, Mk2 model, range across UFBoot trees, 80% taxon jackknife. Reports "repeated changes, range X–Y", never an exact count | Glu/non-Glu changes on the 200-enzyme tree. **This is how to report the unresolved one-vs-two-losses question.** |
| Dimer check: AF-Multimer, 4 seeds. Monomer calls can flip once the dimer is built | Human AdoMetDC is an (αβ)₂ dimer, potato a monomer. Check whether the putrescine site touches the interface; if it does, model dimers for key cases (AF Server, Tui can run it) |
| Discordant-locus check: introns and spliced neighbours to rule out contamination | Ties to to-do 1F (genome vs transcriptome) for the cross-supergroup sequences |

**What's new compared with his paper:**

- The site is allosteric, not catalytic.
- The controls are biochemical (enzymes with measured putrescine responses), not only structural.
- AdoMetDC has known precedents for each outcome: human E178Q, plant filled pocket, *T. cruzi* shifted site, *Plasmodium* Lys-filled pocket.

## Stage 1 = go/no-go (do this first, 1–2 sessions)

1. **Literature check.** Has anyone compared AdoMetDC putrescine-site electrostatics across species? Start with Bale 2008, Bennett 2002, Clyne 2002, Beswick 2006, Wells 2006 / Williams 2011 (notes in `5_literature/paper_notes_activation_mechanisms.md`).
2. **Choose the measurements before looking at any stramenopile** (fixed in advance, like the contamination rule):
   - electrostatic potential at the putrescine position (PDB2PQR + APBS);
   - pKa of the acidic site residues, human numbering 15, 174, 178, 256 (PROPKA);
   - hydrogen-bond partners of the residue-178 side chain;
   - basic residues filling the site (13, 111, 285);
   - free cavity volume;
   - docking in-site fraction (pipeline already built, see below).
3. **Test on structures with known answers:**
   - Crystals: human 1MSV (S68A proenzyme with putrescine), potato 1MHM, and a trypanosome AdoMetDC if one exists (check in RCSB, probably *T. brucei* with prozyme).
   - AlphaFold models of enzymes with measured responses:
     - activated, keep Glu: human, yeast, *Neurospora*, *C. elegans*, *Onchocerca*, *Acanthamoeba*;
     - activated, no Glu: *T. cruzi*, *Euglena*;
     - not stimulated: plants, *Plasmodium*;
     - mutants: human E178Q (4× basal, no stimulation, Bale 2008 Table 4) and *C. elegans* E194Q.
4. **Gate:**
   - **Go:** the measurements separate activated from not-stimulated enzymes, *and* AlphaFold and crystal agree for human and potato.
   - **No-go:** write the evolution paper only, and record why in the log.

**Stage 2 (only after a go):** run all 1,000 FOLDTRACE models, map the results onto the 200-enzyme tree, and run robustness checks. Optional stage 3: short MD on 3–4 enzymes (Marc's specialty).

## What already exists (don't redo)

- **Docking pilot:** `3_data/test_results/RESULTS_docking.md`, scripts `4_scripts/docking_*.py`, all files in `8_full_workspace_backup/workspace_part5_docking.zip`.
  - Method: superimpose AlphaFold models on 1MSV; Vina 1.2.5; box 22 Å and 12 Å centred on the putrescine centroid (61.324, 11.189, 7.513); "in site" = within 2.5 Å.
  - Controls passed: human, *Acanthamoeba* and lamprey dock in site; plant filled; *T. cruzi* shifted 6 Å.
- **Known pitfalls:**
  - AutoDock4 scoring failed its human control, so don't use it.
  - AlphaFold apo side chains can close a pocket (the haptophyte keeps Glu but shows a filled pocket), so robustness checks are needed.
  - The *T. cruzi* model lacks its first ~20 residues.
  - The first superposition script applied the transform twice; that bug is fixed in `docking_prep.py`.
- **Residue table by lineage** (human numbering): `5_literature/paper_notes_activation_mechanisms.md`.
- **Evolution results:** `README.md` (findings 1–8), `1_reports/RESEARCH_LOG.md`, `3_data/test_results/RESULTS_round4.md`.

## Where things are on Tui's computer

- **Repo (read/write):** `C:\Users\Soniya\Documents\GitHub\adometdc-glu178`.
- **FOLDTRACE AlphaFold models (READ-ONLY, 1,000 PDB files named by UniProt accession):** `C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc\models\`

## Ground rules

- **Read-only on Marc's FOLDTRACE project** and the `unlabeled_fp_search` copy. Nothing goes to Marc without Tui.
- **Tui commits and pushes with GitHub Desktop.** Claude only writes files into the repo folder. Commit in a separate step *after* edits, and check `expectedMtimeMs` first.
- **Save usage:** write scripts Tui can run on Tui's own Windows PC (Python is installed) and bring results back; avoid long polling. Install missing tools in the cloud workspace (`autodock-vina` via apt, `openbabel`, `rdkit`, `pdb2pqr`, `propka`, APBS); they don't carry over between chats.
- **Web search is off for the organization.** Use the desktop app's built-in browser for RCSB, PubMed and NCBI; the cloud workspace can't reach them directly.
- No paper PDFs in the repo (copyright). Verify facts against sources. Don't oversell novelty. Predictions stay "plausible", not "proven".
- Plain-language updates to Tui at each step, and a "what failed and what I did" list at the end.
