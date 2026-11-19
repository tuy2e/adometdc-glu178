# Pre-registered measurements and go/no-go rule: pocket-chemistry survey, Stage 1

**Written 10 Nov 2026, before any structure was measured.** Nothing in this file may be changed
after the Stage 1 numbers are seen. If something here turns out to be unworkable, the change is
recorded as a dated amendment at the bottom, with the reason, not by editing the text above it.

This follows the same discipline as the contamination rule in round 4: the decision criteria are
fixed while the answer is still unknown.

---

## 1. The question Stage 1 answers

Not "are stramenopiles activated by putrescine?" — Stage 1 does not touch a single stramenopile.

Stage 1 asks only: **do structure-derived pocket measurements separate AdoMetDCs whose putrescine
response has already been measured in a lab?** If they cannot separate enzymes where the answer is
known, they cannot be trusted where it is unknown, and the survey stops.

## 2. Reference frame

All measurements are made in the frame of **1MSV chain A** (human S68A proenzyme with putrescine
bound). The **putrescine position** is the centroid of the crystal putrescine heavy atoms,
`(61.324, 11.189, 7.513)`, as used in the 10 Nov docking round. Models are superimposed onto 1MSV
by the method already validated in `4_scripts/docking_prep.py` (Cα of human residues within 18 Å of
putrescine, mapped by pairwise MAFFT, core refined twice to atoms under 2 Å).

Residue positions are given in **human numbering** throughout, mapped per structure by the same
pairwise alignment.

## 3. The six measurements

Fixed in advance. No measurement is added, dropped or re-weighted after the results are seen.

| # | Measurement | How | Reported as |
|---|---|---|---|
| M1 | **Electrostatic potential at the putrescine position** | PDB2PQR 3.x (PARSE force field, pH 7.0, `--titration-state-method=propka`) → APBS linearised Poisson–Boltzmann, 150 mM 1:1 salt, solute dielectric 2, solvent 78.54, 298 K | mean potential in kT/e over a 2 Å sphere centred on the putrescine position, and the value at the centroid itself |
| M2 | **pKa of the acidic site residues** | PROPKA3 on each prepared structure | predicted pKa at human positions 15, 174, 178, 256 (blank where the residue is not acidic) |
| M3 | **H-bond partners of the 178 side chain** | heavy-atom donor–acceptor ≤ 3.5 Å with D–H···A angle ≥ 120°, computed on the PDB2PQR-protonated structure | list of partner residues, and the count |
| M4 | **Basic residues filling the site** | any Arg/Lys/His side-chain nitrogen within 8 Å of the putrescine position | count, plus which of human 13, 111, 285 are basic |
| M5 | **Free cavity volume** | 3 Å probe, free volume within 4 Å of the putrescine position (identical to the 10 Nov docking round, so the earlier numbers are directly comparable) | Å³ |
| M6 | **Docking in-site fraction** | AutoDock Vina 1.2.5, putrescine dication, focused 12 Å box on the putrescine position, 9 modes; "in site" = pose heavy-atom centre within 2.5 Å of the crystal putrescine centre | poses in site out of 9, and the best Vina score with its distance from the site |

**Expected directions**, stated now so that a result in the wrong direction counts against the
method rather than being reinterpreted afterwards: a pocket that binds the putrescine dication
should have a **more negative** M1, a **larger** M5, and a **higher** M6 than one that does not.
M2–M4 are descriptive and are not used to pass or fail the gate on their own.

## 4. The test set (fixed before measuring)

Every structure here has a **published, measured** putrescine response, or is a crystal of one.

**Class A — activated by putrescine, Glu178 kept**

| Structure | What | Source |
|---|---|---|
| 1MSV | human S68A proenzyme + putrescine (crystal) | Tolbert 2003 |
| 1JL0 | human H243A proenzyme + putrescine (crystal) | Ekstrom 2001 |
| 1I7B | human processed enzyme + putrescine (crystal) | Tolbert 2001 |
| P17707 | human (AlphaFold) — crystal-vs-model control | — |
| P21182 | *Saccharomyces cerevisiae* (AlphaFold) | — |
| Q1K889 | *Neurospora crassa* (AlphaFold) | — |
| O02655 | *Caenorhabditis elegans* (AlphaFold) | — |
| Q27883 | *Onchocerca volvulus* (AlphaFold) | — |
| L8GRS6 | *Acanthamoeba castellanii* (AlphaFold, already aligned) | — |
| A0AAJ7UEA7 | lamprey *Petromyzon marinus* (AlphaFold, already aligned) | — |

**Class B — activated by putrescine, no Glu178** (the case that proves activation can survive the loss)

| Structure | What |
|---|---|
| 5TVO | *Trypanosoma brucei* AdoMetDC monomer (crystal, 1.48 Å; Volkov 2016) |
| 5TVF | *T. brucei* AdoMetDC/prozyme heterodimer (crystal, 2.42 Å) |
| O76240 | *T. cruzi* full-length (AlphaFold) |
| A0A2V2WHV6 | *T. cruzi* (AlphaFold, already aligned; **truncated**, missing ~20 N-terminal residues) |

**Class C — not stimulated by putrescine**

| Structure | What |
|---|---|
| 1MHM | potato (crystal; Bennett 2002) |
| Q04694 | potato (AlphaFold) — crystal-vs-model control |
| D8S2S4 | *Selaginella* (AlphaFold, already aligned) |
| Q9U8D4 | *Plasmodium falciparum* AdoMetDC–ODC (AlphaFold; AdoMetDC domain only) |

**Class D — reported but NOT used in the gate**

| Structure | What | Why excluded from the gate |
|---|---|---|
| human E178Q | in-silico mutant of 1MSV | 4× basal, no stimulation (Bale 2008 Table 4). Its defect is in **signalling**, not binding, and the 10 Nov docking round showed it docks exactly like wild type. No static measurement can capture it. |
| O02655 E194Q | in-silico *C. elegans* mutant | same reason |

Excluding Class D from the gate is itself a pre-registered decision, so that "the mutants didn't
separate" cannot later be waved away.

**No stramenopile, cryptophyte, oomycete, haptophyte, *Perkinsus* or dinoflagellate structure is
measured in Stage 1.** They are the test set for Stage 2 and are not looked at until the gate has
been decided.

## 5. The gate

**G1 — separation.** On at least one of M1, M5, M6, the Class C (not stimulated) values and the
combined Class A + B (activated) values must **not overlap**, in the direction stated in §3, and the
gap between the two classes must be **larger than the method's own noise**, defined in G2.

**G2 — the method's own noise.** For the two enzymes where both a crystal and an AlphaFold model
exist (human: 1MSV vs P17707; potato: 1MHM vs Q04694), the crystal-to-model difference on each of
M1, M5, M6 defines the method's error bar. The three human crystals (1MSV, 1JL0, 1I7B) and the two
*T. brucei* crystals (5TVO, 5TVF) give a second, crystal-only error bar for the same measurements.
**The G1 gap must exceed the larger of these two.** For M6, agreement means the crystal and the
model fall on the same side of the in-site/not-in-site call.

**Go** = G1 and G2 both pass. Proceed to Stage 2: all 1,000 FOLDTRACE models, mapped onto the
200-enzyme tree, plus robustness checks.

**No-go** = either fails. Write the evolution paper only (core claims: findings 1–3), record the
Stage 1 numbers and the reason in `1_reports/RESEARCH_LOG.md`, and keep the pocket-chemistry result
as a negative methodological note. A no-go is a real result and gets written up as one.

**Partial** = G1 passes but only for the **filled-pocket** contrast (Class C plants vs everything
else), while Class B (*T. brucei*, *T. cruzi*) sits with the activated enzymes on M5/M6 but not on
M1. This is the most likely outcome given the 10 Nov docking round, and it is **not a go**: it would
mean the measurements detect "pocket filled in" but cannot distinguish "not activated" from
"activated through a shifted site", which is exactly the gap the survey was meant to close. In that
case Stage 2 is restricted to the plant-like/filled hypothesis for oomycetes, and the stramenopile
question stays open.

## 6. Known limits, stated now

- AlphaFold models are apo and side chains are unrelaxed; a model pocket can close that is open in
  reality (the haptophyte, 10 Nov). G2 is the explicit check on this.
- Proenzyme and processed forms are mixed in the test set (1MSV/1JL0 proenzyme, 1I7B processed,
  5TVO/5TVF processed). Processing reshapes the site. The human proenzyme-vs-processed comparison
  is itself informative and is reported.
- The *Plasmodium* structure is the AdoMetDC domain of a 1,419-residue fusion; its model confidence
  at the site is reported and the structure is dropped from the gate if the site pLDDT is below 70.
- PROPKA and continuum electrostatics both assume a single rigid conformation; neither is a free
  energy.
- No measurement here tests signalling, only the site. This is the standing limitation of the whole
  approach and belongs in the paper's Limitations section whatever Stage 1 returns.

## Amendments

**Clarification, 10 Nov 2026, written before any structure was measured.** Implementation details
that §3 left implicit, recorded here rather than silently decided at the keyboard:

1. **Everything that is not protein is removed before PDB2PQR**: waters, ions, buffer molecules,
   inhibitors, and the bound putrescine itself in 1MSV, 1JL0 and 1I7B. M1 is therefore the
   potential of the **empty** site. Leaving the dication in would make its own +2 charge dominate
   the number and would not be comparable with the apo models.
2. **Multi-chain structures keep all their chains** for M1, M4 and M5, because a neighbouring
   subunit is part of the site's electrostatic environment. This matters for 1I7B (processed
   α + β chains), 5TVF (AdoMetDC + prozyme) and 5TVO. Superposition uses only the chain with the
   highest identity to human AdoMetDC; which chain that was is recorded per structure.
3. **Protonation at pH 7.0** throughout, via PDB2PQR with PROPKA titration states. His tautomers
   are left to PDB2PQR.
4. The APBS grid is centred on the putrescine position with a fine-grid spacing of ≤0.5 Å, inside
   a coarse grid that encloses the whole solute.

**Amendment 1, 11 Nov 2026, before any Stage 1 structure was measured.** Reason: an overnight
literature check (`BACKGROUND_factors_pocket_chemistry_survey.md`) found two facts that make the
original §5 wording wrong.

1. **5TVO and 5TVF are removed from the G2 noise estimate.** Volkov et al. 2016 (*eLife* 5:e20198)
   show that monomeric *T. brucei* AdoMetDC is autoinhibited by its own N-terminal sequence and that
   prozyme heterodimerisation relieves it. 5TVO is therefore the **autoinhibited** state and 5TVF the
   **activated** state of the same enzyme, not two draws from the same distribution. Treating their
   difference as measurement noise would inflate the error bar with a real conformational change and
   make the gate too easy to pass. They are reported as a **state comparison**. The crystal-to-crystal
   error bar in G2 now comes only from the human crystals.
2. **1JEN is added to the test set** (human AdoMetDC, Ekstrom 1999, 2.25 Å, **no ligand bound**).
   Bale 2008 shows that putrescine binding closes loop 312–320 over the site; 1MSV, 1JL0 and 1I7B are
   all holo, while every AlphaFold model is apo, so "crystal vs model" as originally written also
   measures "closed lid vs open lid". 1MSV-vs-1JEN isolates the lid effect and 1JEN-vs-P17707 gives a
   clean apo crystal-to-model comparison.

**Not amended.** The six measurements and the gate thresholds are unchanged. The signalling-relay and
aromatic-lid positions (human 11, 80, 229, 243 and 285, 315, 318, 320) found informative in that same
overnight work are deliberately **not** added to the gate, because adding a measurement after seeing
data is what this file exists to prevent. They are recorded as candidate Stage 2 measurements.
