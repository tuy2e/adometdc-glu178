# AlphaFold Server dimer results

Drop the downloaded `.zip` files from AlphaFold Server here, one per job.
Do not unzip them — `p13` handles that.

Expected, seven jobs:

| # | protein | organism | why |
|---|---|---|---|
| 01 | A0A1Y1HLG4 | *Klebsormidium nitens* | fusion, ODC half intact (13/13) |
| 02 | A0A2V3IZY9 | *Gracilariopsis chorda* | fusion, ODC half degraded; species has no other working ODC |
| 03 | A0A9C7PQ51 | *Galdieria partita* | fusion, ODC half degraded; species does have a standalone ODC |
| 04 | C1E5E7 | *Micromonas commoda* | the fusion already reported (Green 2011) |
| 05 | Q8IJ77 | *Plasmodium falciparum* | the known working bifunctional enzyme |
| 06 | A0AAV8UQI0 | *Rhodosorus marinus* | control: standalone ODC, 13/13 |
| 07 | A0A9D4U7K2 | *Adiantum capillus-veneris* | control: standalone AdoMetDC |

All submitted as homodimers (Copies = 2) on alphafoldserver.com, seed on Auto,
one seed per job. These are a pilot — the publication version would be four seeds
per complex on Marc's cluster, as in Zimmer 2026 JMB.

## The question

ODC's active sites are built across the interface between two subunits, and
AdoMetDC also works as a dimer. So one dimerisation of a fusion protein may serve
both halves at once. And if a degraded ODC half still forms the interface while
having lost its catalytic residues, it is being kept as a structural element
rather than as an enzyme.

## Read PAE first

Inter-chain PAE decides whether the predicted arrangement means anything at all.
If it is high, the geometry is a guess and must not be interpreted, however
confident the per-residue pLDDT looks.
