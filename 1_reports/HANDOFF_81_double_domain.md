# Handoff: the 81 double-domain proteins, and one sandwich protein

This is for a fresh chat. Read it first, then the files listed in §2. Do not assume
any knowledge from other conversations.

Student: Tui (Soniya Hossain Tuy), sophomore at Connecticut College, doing
independent computational research with Prof. Marc Zimmer.

## 1. One paragraph of background

AdoMetDC makes dcSAM for polyamine synthesis. It is made as a proenzyme that cuts
itself in two (a pyruvoyl group forms at Ser68) and works as an (alpha-beta)2 dimer.
In trypanosomatids it is regulated by **prozyme**: a second, catalytically dead
AdoMetDC paralogue that has no activity of its own but switches the real enzyme on
by pairing with it. Prozyme is a separate protein in those organisms.

Tui's project started from the putrescine regulatory site (Glu178 in human
numbering) and has since grown into a broad survey of a 1,000-hit FOLDTRACE run.
That survey is **lead 2** and is being worked in a different chat. Do not touch it.

## 2. Read these before doing anything

In `C:\Users\Soniya\Documents\GitHub\adometdc-glu178`:

- `9_dataset_survey/README.md` — what every script does
- `9_dataset_survey/results/architecture_census_summary.txt` — the architecture census
- `9_dataset_survey/results/architecture_per_protein.tsv` — per-protein domain layout
- `9_dataset_survey/results/p5_loss_vs_fusion_report.txt` — site loss by architecture class
- `9_dataset_survey/scripts/p3_residue_grid.py` and `p5_loss_vs_fusion.py` — the
  existing mapping and statistics methods. Reuse them rather than writing new ones.
- `1_reports/RESEARCH_LOG.md` — the running record

## 3. Where the data is

**Repo (read/write):** `C:\Users\Soniya\Documents\GitHub\adometdc-glu178`

**Marc's FOLDTRACE run (READ ONLY — never write, move or delete):**
`C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc\`
Contains `models/` (1,000 AlphaFold PDB files named by UniProt accession),
`foldtrace_calls.tsv`, `foldseek_hits.m8`, `sites_used.tsv`, `prespecification.json`,
`aim1_settings_used.json`, `REPORT.md`, `1MSV_A.pdb`.

Of the 1,000 hits, **773 pass the fold gate** (TM-score >= 0.5 to human 1MSV
chain A). Everything downstream uses those.

## 4. Ground rules — these are not negotiable

- **`aim1_settings_used.json` is LOCKED** (protocol aim1-v1.0). Do not change the
  TM gate (0.5), the offset threshold (4.0 A), or any other setting. A new search
  with a different query must be filed as its own run, never under Aim 1.
- **Never run `git commit` or `git push`.** Write files into the repo folder and
  stop. Tui commits with GitHub Desktop herself.
- **Nothing goes to Marc without Tui.**
- **No paper PDFs in the repo** (copyright). Notes and citations only.
- **pLDDT >= 70** for any residue-level claim. Below that, call it unresolved, not
  present or absent. This is the rule FOLDTRACE itself uses.
- Verify facts against sources. Do not oversell novelty. Predictions stay
  "plausible", never "proven". Say plainly when numbers are too small to conclude.
- Tui prefers to run long jobs herself in a separate PowerShell window. Write the
  script, tell her the command, let her run it and paste back the output. Claude
  Code's own command timeout is 10 minutes.

## 5. File naming — important, two chats share this repo

Another chat is actively working on lead 2 and uses the `p1_`...`p8_` prefix and
`9_dataset_survey/results/`.

**Use `q1_`, `q2_`, ... and write results to `9_dataset_survey/results_81/`.**

Never run a job while the other chat is running one. One at a time.

## 6. Job 1 — the 81 double-domain proteins

Among the 773 fold-passers, **81 have two separate PF01536 (AdoMetDC) matches in
one chain**. Nobody has checked what that means.

Two possibilities:
- Pfam split one domain into two matches — an annotation artefact, nothing there.
- They really are two AdoMetDC units fused into one polypeptide.

If it is the second, the obvious question is whether one unit is dead. That would
be a **single-chain prozyme**: the active enzyme and its dead regulator joined into
one protein, instead of the two separate proteins seen in trypanosomatids. As far
as Tui and this project know, no such thing has been described. Check that claim
before repeating it.

What to measure:

1. Are they two folded units? Superpose each PF01536 region separately onto
   1MSV chain A and report TM-score and coverage for each. One domain split in two
   gives two partial matches; two real units give two good ones.
2. Does each unit keep the processing and catalytic machinery? Map Marc's five
   prespecified sites onto each unit separately: **Ser68** (pyruvoyl self-cleavage),
   **Ser229** and **His243** (processing proton acceptors), **Cys82** (catalytic
   proton donor), **Glu247** (substrate positioning). Report per unit, with pLDDT.
   The interesting pattern is one complete unit and one degraded unit in the same
   chain.
3. Do the two units pack against each other the way the real AdoMetDC dimer does?
   Superpose the pair onto the processed human dimer (1I7B, chains A and B) and
   report how well it fits.
4. Check PAE from AlphaFold DB between the two units before making any claim about
   their relative arrangement. AlphaFold is unreliable on domain orientation across
   a floppy linker. If PAE is high, say the arrangement is unsupported and stop.
5. Taxonomy: which organisms, and does any of them also have a separate standalone
   AdoMetDC?

**Known confound, state it in the report.** This class shows the highest site-loss
rate in the dataset — 13 of 81 (16%) lost at least one prespecified site, against
4.2% for standalone proteins. That is probably circular: a degraded domain is more
likely to get split into two Pfam matches in the first place. Do not report the
16% as a finding without addressing this.

## 7. Job 2 — the sandwich protein (small, do it second)

**A0A8S1GM37**, *Caenorhabditis auricularia*. Architecture `HMG-box > AdoMetDC > VTT`.
The only protein in all 773 with AdoMetDC in the middle of two other domains.
A domain in the middle needs either two fusion events or an insertion into an
existing gene, which is why it was flagged.

This is **one protein**, so it cannot carry a project. Spend an hour, not a week:

1. Is the gene model real? Genome or transcriptome assembly, exon structure, any
   intron spanning either junction. Does *C. auricularia* also have a normal
   standalone AdoMetDC? Do other *Caenorhabditis* species have this architecture?
2. Does the AdoMetDC part keep the five prespecified sites, and at what pLDDT?
3. Quick literature check: has anyone reported an HMG-box or VTT domain fused to
   AdoMetDC, in any organism?

If it is a single transcriptome-only gene model with no relatives, say so and close
it. That is a perfectly good outcome.

## 8. Reporting

Write a markdown report per job into `9_dataset_survey/results_81/`, plus the TSVs.
End each with a short "what failed and what I did instead" list. Keep Tui updated in
plain language as you go — she is not a programmer and dislikes jargon and padding.
