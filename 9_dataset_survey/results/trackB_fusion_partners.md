# Track B: fusion partners of AdoMetDC (17 Nov 2026)

Architectures with >= 2 proteins, from `fusion_architectures.tsv`. Partner descriptions are the Pfam
entry text (InterPro API, InterPro 111.0). Literature searched in Europe PMC (abstracts and open full
text) on 17 Nov 2026; "not found" means none of the searches below returned a report, not that none exists.

## 1-2. Partner domains

| Architecture | n (fold-pass / InterPro-source) | Partner | Pfam description / function (one line) | Group |
|---|---|---|---|---|
| PF01536 > PF02784 > PF00278 | 11 (11 / 9) | PF02784 + PF00278 | PLP-dependent decarboxylase (TIM-barrel PLP-binding domain + C-terminal sheet domain) acting on ornithine, lysine, arginine and related substrates | polyamine |
| PF01536 > PF02784 | 7 (7 / 6) | PF02784 | same family, barrel domain only | polyamine |
| PF01536 > PF00278 > PF02784 > PF00278 | 3 (3 / 3) | PF02784 + PF00278 | same family | polyamine |
| PF08132 > PF01536 | 36 (32 / 28) | PF08132 | AdoMetDC leader peptide: product of the conserved small upstream ORF in plant AdoMetDC mRNA 5' leaders; represses AdoMetDC translation | polyamine (AdoMetDC's own regulatory peptide) |
| PF08132 > PF01536 > PF01536 | 2 (2 / 2) | PF08132 | as above | polyamine |
| PF01536 > PF08132 | 2 (1 / 0) | PF08132 | as above, but C-terminal; UniParc-feature data only | polyamine |
| PF00566 > PF08132 > PF01536 | 3 (3 / 3) | PF00566 (+ PF08132) | Rab-GTPase-TBC domain: GTPase-activating protein for Rab-like small GTPases (vesicle traffic) | **no obvious polyamine link** |
| PF01536 > PF09335 | 3 (3 / 3) | PF09335 | VTT domain (VMP1/TMEM41/TVP38): membrane domain linked to vesicular trafficking; formerly "SNARE-associated" | **no obvious polyamine link** |
| PF10609 > PF01536 | 2 (2 / 2) | PF10609 | NUBPL iron-transfer P-loop NTPase: Fe-S cluster assembly factor (NBP35) / plasmid-partition ATPase family | **no obvious polyamine link** |
| PF06441 > PF01536 | 2 (2 / 2) | PF06441 | N-terminal region of eukaryotic epoxide hydrolase (adds water to epoxides) | **no obvious polyamine link** |
| PF01536 > PF01536 (79), PF01536 x3 (5) | 84 | none | two or three separate PF01536 matches in one chain | not a partner |
| PF02675 | 2 (0 / 1) | none | bacterial/archaeal (class 1) AdoMetDC alone | not a partner |
| PF36262 | 7 (0 / 7) | none | TIGR04255 protein, bacterial; Pfam: "structural similarity to AdoMetDC (PF01536), function unknown" | no AdoMetDC Pfam |
| PF33883 | 5 (0 / 5) | none | DUF9824, bacterial (mostly Alphaproteobacteria), function unknown | no AdoMetDC Pfam |
| PF12092 | 2 (0 / 0) | none | DUF3568, bacterial; includes Francisella lipoprotein Flpp3, function unknown | no AdoMetDC Pfam |

Organisms of the non-polyamine-partner fusions:
- TBC > leader > AdoMetDC: *Camellia sinensis* A0A4S4DFF5, *Hevea brasiliensis* A0A6A6MN09, *Rhododendron williamsianum* A0A6A4MJE6
- AdoMetDC > VTT: nematodes *Parelaphostrongylus tenuis* A0AAD5WKP5, *Trichostrongylus colubriformis* A0AAN8F2B4, *Caenorhabditis remanei* A0A261AI60 — ⛔ **probable gene-merge artefacts, see the literature table below (withdrawn 19 Nov)**
- NUBPL > AdoMetDC: Sordariomycetes *Staphylotrichum tortipilum* A0AAN6RUM9, *Emericellopsis cladophorae* A0A9P9XYK1
- Epoxide hydrolase N > AdoMetDC: Dothideomycetes *Hortaea werneckii* A0A1Z5T579, *Stemphylium lycopersici* A0A364NFZ2

## 3-4. Literature status

| Fusion | Status | Source(s) | What the papers say it does |
|---|---|---|---|
| AdoMetDC (N) - ODC (C), *Plasmodium* | **reported** (not in this dataset) | Muller 2000 JBC (PMID 10713131); Krause 2000 Biochem J (11085920); Wrenger 2001 JBC (11390378); Birkholtz 2004 Biochem J (12974675); Hart 2016 Biol Open (27387533) | One ORF, AdoMetDC N-terminal, ODC C-terminal, joined by a hinge; recombinant protein has both activities; ~330 kDa heterotetrameric complex. Wrenger 2001: the two active sites work independently (no domain-domain interaction), half-life > 2 h, putrescine does not stimulate the AdoMetDC domain. Birkholtz 2004: parasite-specific inserts and the hinge are needed for complex formation, and AdoMetDC-ODC interactions are "vital for optimal ODC activity" (this differs from Wrenger 2001). Hart 2016: deletion in *P. yoelii* abolishes male gametogenesis and blocks mosquito transmission. |
| AdoMetDC (N) - ODC (C), *Micromonas* | **reported**; C1E5E7 in this dataset is the reported *Micromonas* sp. RCC299 protein (XP_002502359) | Green 2011 Mol Microbiol (21762220) | Sequence identification only (activity was measured for the bacterial fusions in that paper, not for *Micromonas*). Both domains most similar to green-alga/plant genes, so an origin independent of *Plasmodium*; mRNA has a small upstream ORF. |
| AdoMetDC - PLP decarboxylase, the other 22 proteins (trebouxiophyte, prasinophyte and other green algae, *Klebsormidium*, *Closterium*, ferns, liverwort, red algae, *Percolomonas*, *Vitrella*) | **not found** | searched: AdoMetDC/SAMDC + ODC + fusion/bifunctional (abstracts, 15 hits, all *Plasmodium*/*Micromonas*/unrelated); "AdoMetDC/ODC"-type terms excluding *Plasmodium* (47 hits, none); + chromerid/*Vitrella*/*Perkinsus*/dinoflagellate; + *Chlorella*/fern/*Galdieria*/red alga; Lin & Lin 2018 microalgae review (full text: only AdoMetDC-APT fusions); Fuell 2010 abstract | - |
| Leader peptide (PF08132) + AdoMetDC in one protein | **uncertain** | Franceschetti 2001 Biochem J (11139406); Hanfrey 2005 JBC (16176926); Ivanov 2010 NAR (19920120) | These papers describe PF08132 as the product of a separate small uORF in the 5' leader (overlapping a tiny uORF), translated on its own and repressing translation of the main AdoMetDC ORF in a polyamine-responsive way. No report found of a single protein containing both. Whether the 43 database proteins are real translation products or gene models that include the leader is not established here. |
| Rab-GAP TBC + AdoMetDC | **not found** | Europe PMC full-text co-mention search (8 hits, all unrelated expression/proteomics papers) | - |
| AdoMetDC + VTT | ⛔ **probable annotation artefact — withdrawn 19 Nov** | full-text co-mention search (9 hits, all unrelated); lead-3 chat, `results_81/q1_sandwich_report.md` | Not a literature question. *Caenorhabditis remanei* has **one genome annotation carrying the fusion and another that does not**, which is the signature of gene-prediction software merging neighbouring genes. The same VTT domain appears in the lead-1 sandwich protein (`A0A8S1GM37`), which was independently shown to be a gene merge. Treat all three nematode AdoMetDC–VTT entries as artefacts unless re-examined. |
| NUBPL P-loop NTPase + AdoMetDC | **not found** | full-text co-mention search (262 hits on generic terms; first 25 unrelated) | - |
| Epoxide hydrolase N + AdoMetDC | **not found** | full-text co-mention search (30 hits, all unrelated) | - |
| Michael 2011 FEBS Lett 585:2627 (PMID 21827754) | covers a **different** fusion family | abstract + earlier full reading (see `1_reports/RESEARCH_LOG.md`) | Diatom fusions of **bacterial-type** AdoMetDC with an aminopropyltransferase (+/- histone methyltransferase, Tudor), proposed to make and transfer long-chain polyamines for silica. None of these architectures (bacterial PF02675 fused to an aminopropyltransferase) occurs in this dataset. |

Other Michael-lab papers checked for fusion statements: Green 2011 (full text), Li 2023 JBC (class-2 AdoMetDC arose from a fusion of two class-1b genes; bacterial SpeD-ODC pyruvoyl fusions), Li 2026 JBC (animal spermine synthase derives from a bacterial AdoMetDC-aminopropyltransferase fusion), Li 2023 PNAS (virus ODCs). None mentions any of the non-ODC partners above.

## 5. Which decarboxylase is the PLP-decarboxylase partner? (`p6_decarboxylase_partner.tsv`)

Method: local alignment (BLOSUM62) of the region after the AdoMetDC domain against 14 references of
known activity (ODC: human, *T. brucei*, yeast, tobacco, *Plasmodium* bifunctional; ADC: *Arabidopsis*
ADC1/2, *E. coli* SpeA, chlorella virus PBCV-1 ODC-fold ADC; DAP decarboxylase LysA; CANSDC x2;
antizyme inhibitors AZIN1/2), plus the ODC catalytic residues Lys69 and Cys360 read via human ODC.
23 proteins (the 21 in the three architectures above, plus 2 singletons with the same domains).

| Result | Proteins |
|---|---|
| Closest to an ODC (44-54% identity), ODC Lys69 and Cys360 both present, score margin >= 117 over the best non-ODC family | 16: *Klebsormidium*, *Myrmecia*, *Adiantum*, *Riccia*, *Micractinium*, *Apatococcus lobatus*, *Cymbomonas*, *Tetraselmis*, *Micromonas* C1E5E7, *Closterium*, *Chloropicon*, *Prototheca*, *Apatococcus fuscideae* A0AAW1T1C3, *Pycnococcus*, *Ceratopteris*, *Vitrella* |
| Closest to ODC, Cys360 present, Lys69 position falls in an alignment gap | 1: *Auxenochlorella* A0A087SFF4 |
| Closest to ODC but partner truncated (Cys360 region absent; margin 40-47) | 2: *Auxenochlorella* A0A1D2AC22 (547 aa), *Apatococcus fuscideae* A0AAW1STH3 (411 aa) |
| ODC Lys69 replaced (Ser or Arg); closest references ODC or AZIN1 at 28-32% identity with small margins (14-55); identity not determined | 4: *Percolomonas* A0A7S1KLJ4 (S), *Galdieria* A0A9C7PQ51 (R; Cys360 also replaced by Gly), *Rhodosorus* A0AAV8UZF7 (R), *Gracilariopsis* A0A2V3IZY9 (R) |
| Closest to ADC (plant/bacterial), DAP decarboxylase or CANSDC | 0 |

Caveat stated in the script: similarity places these in the ODC/ODC-fold group, but an ODC-fold enzyme
can prefer arginine (PBCV-1, Shah 2004), so substrate is not established by this comparison. A check
of the residue Shah 2004 names (T. brucei ODC position 332) was attempted and **not used**: UniProt
P07805 has Ala at 332, so their numbering does not match that sequence.
