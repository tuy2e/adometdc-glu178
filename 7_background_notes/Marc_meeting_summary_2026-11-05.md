# Summary for Marc — 5 Nov 2026
*What I've done since the Aim 1 update, what's solid, what's pending your call, and proposed next steps.*

## 1. Aim 1 re-run (your protocol) — done
- Re-ran the full six-step pipeline on the new reference **1MSV with 5 sites** (Ser68, Ser229, His243, Cys82, Glu247).
- **Result is a null:** assembly-based loss ~2.3% vs chemistry-based ~1.9%, with overlapping confidence intervals — the hypothesis isn't supported by this run. One nuance: Glu247 came out "altered."
- Output packaged (`aim1_stuy_adometdc_20261028.zip`); kept the old 1I7B run and zip intact. **I haven't sent it to you yet — tell me if you want it.**

## 2. Fact-check of the project files — done, a few real errors found
- Checked the reports and research notes against primary sources (NCBI/UniProt/PMC/RCSB).
- Fixed verified errors in the reports (e.g. a wrong table citation; restored the correct K80A kinetic numbers).
- Found several **citation errors in the research notes** (wrong PMIDs/PMC IDs/DOIs, a non-existent E. coli structure claim, a mis-assigned enzyme) — documented with the corrections, **not yet applied** pending your OK.
- Corrected report files are staged; **nothing committed to the repo** — your call on whether I commit or you drop them in.

## 3. Research direction — polyamine "prozymes" (dead-enzyme activators)
This is the thread I went deep on. Three parts:

**(a) Literature grounding.** Did a full review of what the catalytically dead paralog actually *does* and why it's conserved. Short version: in trypanosomes the dead copy is a dedicated **allosteric activator** — it binds the active enzyme and flips it out of an autoinhibited state (Volkov 2016), raising activity ~1,200×. The same dead-paralog trick evolved independently for **three** enzymes in these parasites (AdoMetDC, DHS, PRMT1), and pseudoenzymes are a recognized ~5–10% of proteomes kept for exactly these regulatory roles.

**(b) Free-living angle + a concrete new result.** I mapped the finding that free-living *Bodo saltans* and basal *Paratrypanosoma* each have an active AdoMetDC **plus a dead paralog** (and a parallel dead DHS). Then I did a residue-level test: **do the free-living dead copies keep the activation machinery?** They largely do — the Bodo dead copy retains the two experimentally-validated interface residues (T. brucei M148′/Y152′), and the Bodo active enzyme keeps the P31 conformational switch that even humans lack. For DHS, the Bodo dead copy has lost the catalytic lysine but kept the substrate-binding residues — the same signature as the validated T. brucei DHSp. **This predicts the activation mechanism predates the parasites.**
  - *Honest caveat:* the general idea that this machinery predates parasitism in Bodo is **already partly published** (for the PRMT prozyme). My specific AdoMetDC + DHS mapping with a dated phylogeny looks unreported, but that needs confirming against two paywalled genome papers (Jackson 2016; Flegontov 2013) — I've requested them.
  - *Ready to go:* an AlphaFold-Multimer run of the Bodo heterodimer is prepped (inputs + an analysis checklist) to test the interface in 3D.

**(c) A possible alternative lead — diatoms.** Separately, FOLDTRACE flagged that diatoms + brown algae (Ochrophyta) share a convergent **"stuck-on" AdoMetDC** — a *second, independent* route to putrescine-independence, different from how plants do it. It has a crisp testable prediction (human E178Q is a known stuck-on mutant), a genetically tractable model organism (*Phaeodactylum*), and — on first check — cleaner novelty (only one adjacent paper, on a different protein). I'd flag this as a candidate for a more independent project.

## 4. What I'd like to ask you
- Do you want the Aim 1 null-result package, and how should I frame the Glu247 "altered" nuance?
- OK to apply the research-note citation fixes and to commit the corrected reports?
- Of the two directions — (b) the **Bodo/DHS prozyme mapping** (strongest as a case study / co-authorship on the FOLDTRACE paper) or (c) the **diatom stuck-on convergence** (more of an independent project) — which do you think is worth pushing, and toward what venue?
- For either, the real confirmation is a bench experiment (express the enzyme ± its partner / ± putrescine). Is there a lab you'd point me to for that?

*(Full details and files for any of the above are ready if you want them.)*
