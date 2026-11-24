# ColabFold inputs — homodimer predictions

Marc asked whether this works on Colab. It does: **ColabFold handles multimers**, and a
homodimer is written as the same sequence twice, joined by a colon.

Each `*_COLAB.txt` file holds exactly that — `SEQUENCE:SEQUENCE` on one line, ready to
paste. Nothing to edit.

## The notebook

https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

## Settings

| field | value |
|---|---|
| `query_sequence` | paste the whole contents of one `_COLAB.txt` file |
| `jobname` | the file's number and organism, e.g. `07_human_AdoMetDC` |
| `num_relax` | 0 |
| `template_mode` | none |
| `msa_mode` | mmseqs2_uniref_env (the default) |
| `pair_mode` | **unpaired_paired** (the default; this is what makes the multimer work) |
| `model_type` | **alphafold2_multimer_v3** |
| `num_recycles` | 3 |
| `num_seeds` | **4** — under Advanced settings. This is the whole point of running it here rather than on AlphaFold Server |
| `save_all` | **tick it** — otherwise the PAE arrays are not written out |

Then Runtime → Run all.

## Order to run them — smallest first

Free Colab gives a T4 GPU, and these are large complexes. Size is the main risk.

| run | job | total residues | expectation |
|---|---|---|---|
| 1 | 07 human AdoMetDC (control) | 668 | fine |
| 2 | 06 *Rhodosorus* ODC (control) | 912 | fine |
| 3 | 08 human ODC (control) | 922 | fine |
| 4 | 02 *Gracilariopsis* (ODC half degraded) | 1460 | should be fine |
| 5 | 03 *Galdieria* (ODC half degraded) | 1524 | should be fine |
| 6 | 01 *Klebsormidium* (ODC half intact) | 1728 | probably fine |
| 7 | 04 *Micromonas* | 2078 | borderline |
| 8 | 05 *Plasmodium* | 2868 | likely too big for a free T4 |

**Run the three controls first.** If human ODC does not come back as a proper dimer,
something is wrong with the setup and there is no point running the rest. If it does,
the method is working and everything after it is interpretable.

The comparison that matters is **01 against 02 and 03** — intact ODC half versus
degraded. Get those done before worrying about *Micromonas* or *Plasmodium*.

If a run dies with an out-of-memory error, that is the GPU, not a mistake. Note it and
move on; Colab Pro or Marc's cluster can pick those up.

## What to keep

Each finished run offers a zip download. **Take the whole zip, do not cherry-pick the
best PDB.** It contains the PAE JSON files, and the inter-chain PAE is the number the
entire analysis depends on — a ranked model without it is useless here.

Drop the zips into `9_dataset_survey/data/p16_colab_results/`.

## One methods note for the write-up

ColabFold builds its alignments with MMseqs2 rather than the full AlphaFold database
pipeline. The model weights are the same AlphaFold2-multimer weights. This is standard
and widely published, but it is a difference from AlphaFold-Multimer as run in
Zimmer 2026 and should be stated in the methods rather than glossed.
