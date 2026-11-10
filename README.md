# AdoMetDC Glu178 project: research workspace

Tui's computational research in Prof. Marc Zimmer's lab (Connecticut College), done with Claude (Cowork).
**Last updated: 10 Nov 2026.**

This repository is deliberately **separate from Marc's FOLDTRACE project**. Nothing here changes his project. His reports, FOLDTRACE results and manuscripts stay in his repository, and this one only refers to them.

> **Labels used in these files:** where a file says "Marc's protocol", it follows his instructions. Anything marked 💡 or "my suggestion" is an idea from our side, not part of his protocol.

---

## The question in one paragraph

AdoMetDC (S-adenosylmethionine decarboxylase) makes the aminopropyl donor for spermidine. In humans and most eukaryotes it is switched on by putrescine, and the residue that senses putrescine is **Glu178**. We found that **every stramenopile AdoMetDC has lost Glu178**, and so have cryptophytes, kinetoplastids/euglenids and one dinoflagellate/*Perkinsus* copy. Those Glu-less enzymes sit in unexpected places on the family tree. The main question is how and when the sensor was lost, and what that predicts about how these enzymes are regulated.

## Main findings so far

1. **97 of 97 stramenopile AdoMetDCs lack Glu178** (mostly Asn; His in oomycetes). Cryptophytes, kinetoplastids/euglenids and one dinoflagellate/*Perkinsus* copy also lack it.
2. **Two Glu-less groups break the species tree:** (a) stramenopiles + the dinoflagellate/*Perkinsus* Asp copy, and (b) kinetoplastids/euglenids + cryptophytes. AU tests reject moving any group back to its species-tree place, under both a standard model and a C60 mixture model.
3. **Not settled:** whether (a) and (b) are one lineage with one loss (the best tree) or two groups with two losses (cannot be rejected; this is the arrangement in Li et al. 2019).
4. **Biochemistry:** all 8 sequenced putrescine-activated enzymes outside trypanosomes keep Glu178. Two exceptions, *T. cruzi* and *Euglena gracilis*, lack Glu178 but are still activated. Our ancestral reconstruction shows they inherited the same pocket (Arg13 / Asp174 / Asn178) from a shared euglenozoan ancestor. **The clean prediction ("not switched on by putrescine") is therefore for stramenopiles only; cryptophytes are uncertain.**
5. **Structure:** the diatom pocket is intact and acidic, with Asn where the Glu was, and a Trp next to it (human position 109) that is shared by most Glu-less enzymes.
6. **Novelty:** no paper read so far reports this (about 20 read cover to cover; see `5_literature/`).

The full write-up is `1_reports/stramenopile_Glu178_steps2-6_results.md`. The running record of every step, dead end and decision is `1_reports/RESEARCH_LOG.md`.

---

## Folder guide

| Folder | What's in it |
|---|---|
| `1_reports/` | The main results write-ups (step 1, steps 2–6) and the research log. **Start here.** |
| `2_figures/` | Family-tree figures and the human vs. stramenopile pocket figure. |
| `3_data/sequences/` | The 187 AdoMetDCs in the main tree, the 142 spermidine synthases, *Euglena* sequences from our BLAST searches, and the BLAST query files. |
| `3_data/alignments/` | The main trimmed alignment, plus the robustness versions (pocket columns removed, kinetoplastids removed) and the spermidine synthase alignment. |
| `3_data/trees/` | The main tree (Q.pfam+R7), the C60 mixture-model tree, the constrained trees used in the AU tests, and the spermidine synthase tree. Newick format; open in iTOL or FigTree. |
| `3_data/test_results/` | IQ-TREE reports (models, support values, AU test tables) and the ancestral reconstruction summaries, plus two plain-language results files. |
| `3_data/annotations/` | Which group each sequence belongs to, and its residues at the key positions. |
| `4_scripts/` | The Python scripts used for mapping residues, building constraint trees, ancestral reconstruction and figures. Named by the folder they came from. |
| `5_literature/` | The reading list with notes on every paper read in full, and the literature-search record. Citations and DOIs only; no PDFs. |
| `6_side_project_bodo_prozyme/` | The earlier side thread on a possible AdoMetDC prozyme (dead helper copy) in free-living *Bodo*. |
| `7_background_notes/` | Background on AdoMetDC, a methods walkthrough, meeting notes and the Thursday research package. |
| `8_full_workspace_backup/` | Two zip files holding the complete working directory, so nothing is lost if the session workspace is wiped. See the README inside. |

## What is NOT in this repository, and why

- **Marc's project files** (his reports, FOLDTRACE run outputs, unpublished manuscripts): they stay in his repository.
- **Paper PDFs and extracted paper text:** copyright. The reading list gives full citations and DOIs.
- **IQ-TREE checkpoint and distance files:** large, and regenerated automatically by rerunning the commands in the `.iqtree` reports.

## Live documents (on claude.ai, shared separately)

- **Team doc**, "AdoMetDC project: what I'm doing and where I am": https://claude.ai/artifact/8iaf7suLAcn8Wvfwsx3BKL
- **Project Map** (flow map of what held, changed or was dropped): https://claude.ai/artifact/NxFhcawF69HgBQ5KxdrtUW

Both are private until shared from their Share menu.

## Tools and versions

- IQ-TREE 2.0.7 (ModelFinder, UFBoot, SH-aLRT, AU test with `-zb 10000 -au`, ancestral reconstruction with `-asr`, PMSF C60)
- MAFFT L-INS-i; ClipKIT (gappy 0.5)
- Python 3 with Biopython
- Sequences from UniProt, UniParc and NCBI (protein and TSA via BLAST)
