# Dataset survey: all 1,000 FOLDTRACE candidates

Exploratory survey of every AdoMetDC-like candidate from the `stuy__adometdc` FOLDTRACE run
(reference 1MSV chain A), not only the Glu178 question. Computational only.

The FOLDTRACE run itself is Marc's and is only read from, never changed:
`C:\Users\Soniya\Downloads\unlabeled_fp_search\unlabeled-fp-search\aim1_workflow\runs\stuy__adometdc\`

| Folder | Contents |
|---|---|
| `scripts/` | Python scripts, one per step (`p1_` = Phase 1, `p2_` = Phase 2, ...) |
| `data/` | Downloaded annotation (UniProt, InterPro) |
| `results/` | Output tables |
| `figures/` | Plots |
| `chimerax/` | ChimeraX scripts for the most interesting structures |

## Phase 1 (audit), 15 Nov 2026

- 1,000 candidates; the same accessions in `foldseek_hits.m8`, `foldtrace_calls.tsv` and `models/`.
  Every model sequence equals the Foldseek target sequence.
- 773 pass the fold gate (TM >= 0.5, normalised by the reference); 227 fail and keep `fold_ok = false`.
- **The `foldseek_*` columns of `foldtrace_calls.tsv` are mislabelled** (the web-server m8 was read as
  standard m8: `foldseek_evalue` = target start, `foldseek_bits` = target end, `foldseek_qcov` = E-value,
  `foldseek_alntmscore` = query length). All Foldseek values here come from the raw m8. TM-align
  values, the fold gate and the site calls are not affected.
- Reference 1MSV chain A: residues 4-328, numbering = UniProt P17707, only mutation S68A
  (SEQADV), residues 22-27 have no coordinates. In 1I7B, chain A is the alpha chain only (69-328);
  residues 4-67 are chain B.
- "Unresolved" is not loss: `*_unresolved_reason` says whether no residue aligned, or an aligned
  residue was more than 4 A away, or the protein failed the fold gate.

`results/master_candidates.tsv`: one row per candidate (made by `scripts/p1_build_master_table.py`).

## Phase 2: a third of the hit list no longer resolves in UniProt (17 Nov 2026)

The Foldseek search (AFDB50, October 2026) returned 1,000 AlphaFold models. In UniProt release
2026_04, **326 of those 1,000 accessions (33%) are inactive**, including **211 of the 773
fold-passers (27%)**. All 326 were **deleted**; none were merged or demerged, so UniProt names no
replacement accession. Stated reasons: "Not part of a reference proteome" 319, "Redundant
sequence" 4, "Deleted from sequence source (ENSEMBL)" 3. AlphaFold DB still serves models for
accessions UniProt has since removed, so any AFDB-based hit list will contain them.

The sequences themselves are not lost. All 326 are still in UniParc with exactly the sequence of the
AlphaFold model. `scripts/p2b_resolve_inactive.py` recovers what it can (`results/inactive_resolution.tsv`):
- 12 (10 fold-passers) have an active UniProt entry with the identical sequence (same UniParc ID).
  This is a same-sequence twin, not a UniProt merge, and is labelled `same_sequence_active_entry`.
  Some twins carry an updated species name (e.g. *Shewanella* sp. YLB-09 -> *S. eurypsychrophilus*).
- 266 have domain positions from UniParc sequence features (Pfam, Gene3D, SUPERFAMILY, PANTHER, ...),
  in `results/domain_hits_uniparc.tsv`.
- 48 have no annotation at all; all 48 fail the fold gate.
Every fold-passer therefore has at least domain-level annotation. Inactive accessions are kept in all
tables; nothing has been dropped.

**What the deletions mean.** 319 of the 326 deletions are "Not part of a reference proteome": UniProt
removes unreviewed entries from genomes it does not choose as the reference proteome for that species.
These are real translated proteins, not bad data. But the inactive set is therefore **biased toward
under-studied and less-curated species** (non-reference strains, draft genomes, environmental
genomes), and any pattern carried mostly by inactive accessions is also a pattern of those organisms.

**Rule for headline results.** Every headline result is also reported for the InterPro-only set (the
674 active accessions annotated directly from the InterPro API), and any pattern that appears only in
the replacement-record or UniParc-feature data is flagged.

**Excluded from position-based work:** A0A8C5BDJ7 (UniProt sequence 227 aa, AlphaFold model 245 aa, so
annotation positions do not map onto the model).

## Phase 2c: architecture census (17 Nov 2026; census only, no interpretation)

`scripts/p2c_architectures.py` merges the three annotation sources into `results/domain_hits_merged.tsv`
(column `source` = interpro / replacement_record / uniparc_features; one source per protein, never mixed).
The architecture unit is **Pfam**, the one database reported with the same accessions by InterPro and
UniParc. Regions of >= 40 aa outside every Pfam domain but matched by CATH-Gene3D, SUPERFAMILY, CDD,
SMART or a PROSITE profile are listed as **uncertain segments**, not as domains.

- `results/architecture_per_protein.tsv`: Pfam domains in order with positions, linker lengths,
  AdoMetDC position (only / first / last / middle / repeated), partner domains, uncertain segments,
  and whether the layout is unique.
- `results/fusion_architectures.tsv`: one row per distinct Pfam architecture, with counts split by source
  and fold gate, `only_in_non_interpro_sources`, and organisms for rare (<= 5) architectures.
- `results/architecture_census_summary.txt`: the same counts for all proteins, the InterPro-only set,
  fold-passers, and fold-passers in the InterPro-only set.

Caveats on the census itself: "PF01536 > PF01536" means Pfam reported two separate PF01536 matches in one
chain; whether that is one split domain or two domains is not decided here. UniParc features include
databases InterPro results here do not (FunFam, PIRSR), which is why the architecture unit is Pfam only.

## Scripts, in run order

| Script | What it does | Main outputs |
|---|---|---|
| `p1_build_master_table.py` | one row per candidate from the FOLDTRACE run | `master_candidates.tsv` |
| `p2_fetch_annotations.py` | UniProt + InterPro download (resumable, cached in `data/`) | `annotation_uniprot.tsv`, `domain_hits.tsv` |
| `p2b_resolve_inactive.py` | why inactive accessions are inactive; UniParc + same-sequence recovery | `inactive_resolution.tsv`, `domain_hits_uniparc.tsv` |
| `p2c_architectures.py` | merge the 3 annotation sources; Pfam architecture census | `domain_hits_merged.tsv`, `architecture_per_protein.tsv`, `fusion_architectures.tsv` |
| `p5_loss_vs_fusion.py` | Track A: site loss vs architecture, with lineage controls, pLDDT, per-site and InterPro-only versions | `p5_loss_vs_fusion_report.txt`, `p5_loss_vs_fusion_proteins.tsv` |
| `p6_decarboxylase_partner.py` | Track B q5: nearest reference for PLP-decarboxylase partners + ODC catalytic residues | `p6_decarboxylase_partner.tsv`; literature table in `trackB_fusion_partners.md` |
| `p7_fusion_validation.py` | Track B checks: gene-model source and exon structure (ENA), same-species AdoMetDC/decarboxylase genes, linkers, partner-domain integrity vs human ODC | `p7_*.tsv`; summary in `trackB_fusion_validation.md` |
| `p8_putrescine_site.py` | Track B step 1: putrescine-site positions (13, 15, 111, 174, 176, 178, 256, 285) mapped onto 1MSV as in p3; 23 ODC fusions / all 101 fused / 78 other fused vs 590 standalone; Fisher + within-phylum permutation, also InterPro-only | `p8_putrescine_site.tsv`, `p8_putrescine_site_tests.tsv`, `p8_putrescine_site_report.txt` |
| `p3_residue_grid.py` | TM-align every model onto 1MSV (same aligner as FOLDTRACE); one call per human position | `p3_residue_grid.tsv`, `p3_indel_events.tsv`, `p3_protein_summary.tsv` |
| `p4_variation.py` | where proteins differ from human; pLDDT-filtered rankings, indel hotspots, unique changes, co-changing groups | `p4_*.tsv`, `figures/p4_*.png` |

pLDDT rule: a difference counts only if its region has pLDDT >= 70 in that model; lower ones are
counted in separate "low_conf_excluded" columns and never used in rankings.

## Reference set (agreed 15 Nov 2026)

| Structure | Role | State | Coverage (human numbering) |
|---|---|---|---|
| 1MSV chain A | main reference (same as Marc's run) | holo (putrescine), unprocessed proenzyme, S68A | 4-328; no coordinates for 22-27 |
| AF-P17707-F1 (v6) | unmutated full-length check | apo model, unprocessed, wild type | 1-334 |
| 1JEN chains A+B | apo crystal, for comparison with the apo AlphaFold models | apo (no putrescine), processed (PYR68), no mutations | B 4-67 (no 21-27), A 69-329 (no 165-171, 293-298); second copy in C+D |
| 1I7B chains A+B | processed-form check | holo (putrescine + SMM), processed (PYR68) | B 4-67 (no 24-26), A 69-328 (no 167-170, 293-298); A146G vs UniProt, not yet explained |

In both processed crystals, residue 68 is a pyruvoyl HETATM, so Ser68 itself cannot be mapped from them.
Local copies: `data/references/1JEN.pdb` (downloaded from RCSB).
