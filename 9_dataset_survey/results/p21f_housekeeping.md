# p21 F: housekeeping

## trackB_fusion_validation.md: intron counts (done)

Section 1's two bullets now count the **19 unreported** genome gene models only:

- **9 run through the junction unbroken**: 5 single-exon genes (*Adiantum*, *Chloropicon*, *Gracilariopsis*,
  *Pycnococcus*, *Ceratopteris*) and 4 with the junction inside an exon (*Riccia*, *Apatococcus lobatus*,
  *Cymbomonas*, *Closterium*).
- **10 rely on a predicted intron in the linker** (*Galdieria*, *Klebsormidium*, *Rhodosorus*, *Myrmecia*,
  *Micractinium*, *Auxenochlorella* A0A087SFF4, *Prototheca*, *Apatococcus fuscideae* x2, *Vitrella*).

*Micromonas* C1E5E7 (Green 2011) is stated separately: it is also unbroken across the junction, which is where the
old "10 of the 20" came from. The table in that file was already right and was not changed; only the two bullets
were rewritten, with a line giving 9 + 10 = 19 (20 with *Micromonas*, 23 with the 3 transcriptome-only proteins).

## Thread 5 status line (NOT done)

Thread 5 exists only in `RESEARCH_LOG.md` (the copies at the repository root, in `1_reports/` and in
`Claude outputs/`). The task's standing rule says not to touch RESEARCH_LOG.md, so it was not edited. The two
instructions conflict; tell me which one wins and I will make the change, or edit the status line yourself.

## IQ-TREE (for part C)

No IQ-TREE binary is on PATH or in the usual install folders on this machine. The repository holds output of
earlier IQ-TREE runs (`3_data/test_results/*round4*.iqtree`, and Thread 5 step 4 in RESEARCH_LOG.md), so it was
presumably run elsewhere. `scripts/p21c_tree.py --iqtree <path>` will use a binary wherever it is.
