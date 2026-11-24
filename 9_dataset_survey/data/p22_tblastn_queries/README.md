# p22 tblastn junction queries: the intron-dependent fusions

Built by `scripts/p22_build_inputs.py` from Marc's models (read only; sequence from CA atoms, checked against UniProt) and the domain boundaries in `results/p7_linkers.tsv`. Same recipe as `data/p14_tblastn_queries/`: the **last 60 residues of the AdoMetDC Pfam match**, the **whole linker**, and the **first 60 residues of the PLP-decarboxylase Pfam match**. Positions below are in the full protein.

## Why these two

p21f lists 10 fusions whose gene model puts a predicted intron in the linker. For them, the gene structure says nothing for or against a real fusion, because a gene-prediction merge of two neighbouring genes looks exactly like that. A transcript read across the junction is the test.

| file | protein | type into NCBI's organism box | query length | AdoMetDC / linker / PLP residues | sequence source | genome BioProject(s), from q3 |
|---|---|---|---|---|---|---|
| `01_Galdieria_partita_A0A9C7PQ51.txt` | A0A9C7PQ51 | *Galdieria partita* (taxid 83374) | 234 aa | 260-319 / 320-433 (114 aa) / 434-493 | model sequence identical to cached UniProt | PRJDB12742 |
| `02_Klebsormidium_nitens_A0A1Y1HLG4.txt` | A0A1Y1HLG4 | *Klebsormidium nitens* (taxid 105231) | 225 aa | 275-334 / 335-439 (105 aa) / 440-499 | model sequence identical to cached UniProt | PRJDB718 |

## How to run

blast.ncbi.nlm.nih.gov -> tblastn -> paste one file -> Database: **Transcriptome Shotgun Assembly (TSA)** -> Organism box: the name above (let NCBI autocomplete it) -> BLAST. If nothing comes back, try the genus, then the database **Nucleotide collection (nr/nt)**.

From p15: through the URL API, restricting `tsa_nt` with `ENTREZ_QUERY` left searches waiting for over 50 minutes; the web form restricts TSA through its organism box instead, which is why the web form is suggested here.

## What a positive looks like

- **One subject sequence (a TSA contig or an mRNA) aligning across the whole query in a single HSP**: starting inside the AdoMetDC part, running through the entire linker, and ending inside the PLP part, with no break at the linker. That is one RNA molecule carrying both halves, and here it also shows the predicted linker intron is really spliced out as predicted.
- Hits covering only the AdoMetDC part or only the PLP part are **not** evidence: they can come from separate genes.
- No hit means nothing on its own: the species may have little or no RNA data in TSA (not checked here).
- If the hit's TSA BioProject is the genome's BioProject (table above), the same project built both, so it is **not independent** evidence; record that with the hit.

## Note on the p14 choice

The task calls the p14 choice of single-exon genes backwards. Partly: a single-exon fusion cannot be an intron-spanning merge, so these two intron-dependent genes do need transcripts more. But the single-exon cases are not settled either: `trackB_fusion_validation.md` lists read-through of adjacent genes and assembly errors (a lost stop codon, an indel) as ways an uninterrupted ORF can still be artefactual. Transcript evidence helps both groups.
