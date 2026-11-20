# Literature search: AdoMetDC-ODC fusion (19 Nov 2026)

📋 Requested by Marc (via Tuy): an exhaustive search covering seven items, with every query and its hit count recorded.
Databases: Europe PMC REST API (searches titles, abstracts and the open-access full texts it indexes) and PubMed
(NCBI E-utilities). Full texts were read from the PMC open-access XML (NCBI efetch) and searched sentence by sentence.
The query log is at the end of this file and in `literature_query_log.tsv`.

Terms used throughout:
- **searched, nothing found**: the query ran, and none of the returned records (title, abstract, or full text where
  one was available) mentions the subject.
- **could not search**: the text could not be read at all (paywall, captcha, HTTP 403, or no machine-readable text).
  Nothing is concluded about what that text contains.

A full-text search covers only papers in the PMC open-access subset. For every other paper only the title and
abstract were searched. Supplementary files were searched only where named below.

---

## Summary

| # | Item | Status | Result |
|---|---|---|---|
| 1 | Fusion as a concept | searched | Outside *Plasmodium*, the AdoMetDC-ODC fusion is described in **two** papers: **Green et al. 2011** (*Micromonas*, two species) and **Volkov et al. 2016**, where a database survey reports AdoMetDC fused to ODC Pfam domains in "9 apicomplexan and 3 plant sequences". Neither names the other algae, ferns, liverwort, red algae, *Percolomonas* or *Vitrella*. Separately, Green 2011 describes a bacterial ODC-AdoMetDC fusion and Li et al. 2023 describe bacterial pyruvoyl AdoMetDC-ODC (SpeD) fusions; both are different protein families (see 1c). |
| 2 | *Plasmodium* literature | searched (219 papers: core papers plus all papers citing 11 key ones) | The *Plasmodium* papers call the arrangement "unique to Plasmodia" (Birkholtz 2004: "shared by at least three species"). None of them surveys organisms outside apicomplexans. One 2024 review says "Plasmodium and other apicomplexans ... contain bifunctional enzymes". |
| 3 | Author publication lists (Michael, Hanfrey, Fuell, Elliott, Phillips, Ealick, Pegg) | searched (513 polyamine papers). **The key Michael reviews could not be searched in full text.** | The only fusion statements in these authors' papers are the Green 2011 and Volkov 2016 ones above, plus the bacterial fusion papers. The broad reviews (Fuell 2010, Michael 2011, 2016a, 2016b, 2017, 2018; Pegg 2009, 2018; Phillips 2018) could be read only as abstracts, and none of those abstracts mentions an AdoMetDC-ODC fusion. |
| 4 | The 23 accessions as literal strings, plus their UniProt references | searched, nothing found | The 23 accessions, 24 locus tags and 6 protein IDs give 0 hits in Europe PMC and PubMed. The one exception is *Micromonas* XP_002502359, which points to Green 2011. The 11 UniProt-linked papers were read in full where open access (8): none mentions polyamine enzymes. |
| 5 | Preprints | Europe PMC preprints searched, nothing relevant; **bioRxiv's own search could not be searched** (HTTP 403) | 4 preprint hits, none about the fusion. |
| 6 | The three unread genome papers | partly searched / could not search | *Gracilariopsis* (Lee 2018): about 80k characters of main text read, with no terms; the supplementary zip was not opened. *Micractinium* (Arriola 2018) and *Micromonas* (Worden 2009): **could not search** beyond the abstract (publisher HTTP 403). |
| 7 | Woo 2015 eLife, *Vitrella* | searched (main text, appendices, Figure 2 and all 5 supplements, all 16 supplementary tables) | **The fusion is not shown.** The data are presence/absence of EC numbers per species, with no gene IDs and no domain layout. |

---

## 1. The fusion as a concept

The queries combine (AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine
decarboxylase") AND (ODC OR "ornithine decarboxylase") with each of: fusion, bifunctional, chimeric, fused, hybrid,
multidomain, "single polypeptide", "gene fusion". Each combination was run three ways: Europe PMC title/abstract
only, Europe PMC all fields (includes open-access full text), and PubMed [tiab] (with MeSH for the enzymes).

- Title/abstract hits: 25 in Europe PMC and 25 in PubMed across the 8 words; 27 distinct papers in all. They are the
  *Plasmodium* papers (section 2), Green 2011, Li 2023, and papers where the words are unrelated (cell fusion,
  antisense constructs, transgenic mice).
- Full-text hits: 428 distinct records. Of these, 290 open-access full texts were scanned for sentences that name
  AdoMetDC and ODC together with a fusion word. A further 115 have a PMCID but no open-access text, and 23 have no
  PMCID; for those 138 only the abstract was available. The sentences found came from 32 papers. All of them are
  about *Plasmodium*, or are false matches (tables, abbreviation lists, "fused cultures", cancer constructs), except
  for the papers below.

### 1a. Papers that describe a eukaryotic AdoMetDC-ODC fusion outside *Plasmodium*

| Paper | What it says (verbatim) |
|---|---|
| Green R, Hanfrey CC, Elliott KA, McCloskey DE, ... Michael AJ (2011) *Mol Microbiol*. PMID 21762220, PMC3196669 | "In eukaryotes, we found a gene fusion of a eukaryotic AdoMetDC and an ornithine decarboxylase in two species of single-celled prasinophyte green algae, Micromonas pusilla CCMP1545 (EEH58717, 1019 a.a.) and Micromonas sp. [RCC299]". "Although the Micromonas AdoMetDC-ODC gene fusions are arranged similarly to the Plasmodium AdoMetDC-ODC fusion genes (Muller et al., 2000), these fusion genes represent entirely independent eukaryotic gene fusion events because both domains of the Micromonas fusion protein exhibit greatest homology with the individual component genes from green algae and plants and present only distant similarity to the Plasmodium ..." |
| Volkov OA, Kinch L, Ariagno C, ... Phillips MA (2016) *eLife* 5:e20198. PMID 27977001, PMC5201418 | In the legend of Supplementary file 2, describing their survey of AdoMetDC sequences in NCBI NR: "The confirmed AdoMetDC domains were found to be fused to several additional domains in select sequences, including an AdoMetDC leader (pfam08132) in 35 plant sequences, an f-box-like domain (pfam12937) in one fungal sequence, **pyridoxal-dependent ornithine decarboxylase (pfam02784 and pfam00278) in 9 apicomplexan and 3 plant sequences**, and protein prenyltransferase alpha subunit repeat (pfam01239) in one Blastocystis sequence." |

Notes on Volkov 2016:
- The paper does not name the 3 "plant" sequences.
- Supplementary file 2 (an alignment PDF, eLife CDN, read in full) holds only representative AdoMetDCs, not the fused
  sequences.
- "Plant" here is a taxonomy grouping in NR (Viridiplantae would include green algae such as *Micromonas*). Which
  three sequences were meant cannot be determined from the paper.
- It also counts 35 leader-peptide (PF08132) fusions in plants. This dataset has 43 proteins with that architecture
  (see `trackB_fusion_partners.md`).

### 1b. Reviews that mention the *Micromonas* or *Plasmodium* fusion

| Paper | Statement |
|---|---|
| Pérez-Pertejo Y et al. (2024), PMC10820115 | "... Plasmodium and other apicomplexans, which contain bifunctional enzymes capable of synthesizing simultaneously both polyamines [74]". Which other apicomplexans is not stated. |
| Lin HY, Lin HJ (2018) microalgae review, PMC6356823 | Mentions only diatom SAMDC-APT fusions; no AdoMetDC-ODC fusion (already in `trackB_fusion_partners.md`). |

### 1c. Bacterial ODC/AdoMetDC fusions (different protein families; listed so they are not confused with this one)

- **Green 2011**: an actinobacterial ODC-AdoMetDC fusion, 498 aa in total: "the ornithine decarboxylase ORF has fused
  with the immediately downstream AdoMetDC ORF to form a 498 aa ODC-AdoMetDC fusion protein". At 498 aa it is
  shorter than any of the 23 here. Green 2011 also describes AdoMetDC-aminopropyltransferase fusions.
- **Li B, Liang J, Phillips MA, Michael AJ (2023)** *JBC*, PMC10407285: "fusion proteins of bona fide AdoMetDC/SpeD
  with homologous L-ornithine decarboxylases that possess two ... internal protein-derived pyruvoyl cofactors".
  These are bacterial class-1 AdoMetDC fusions (Omnitrophica). The ODC part is pyruvoyl-dependent, not
  PLP-dependent.

---

## 2. *Plasmodium* AdoMetDC/ODC literature

Searches:
- *Plasmodium* + enzyme terms in title/abstract (Europe PMC 16 hits, PubMed 15).
- Enzyme terms + fusion words in all fields (63).
- The written forms "AdoMetDC/ODC", "PfAdoMetDC/ODC", "ODC/AdoMetDC" and similar (Europe PMC 152, PubMed 12).
- (*Plasmodium* OR malaria) + polyamine in title/abstract (Europe PMC 166, PubMed 154).
- Every paper citing 11 key papers (Europe PMC citations; the counts are in the query log).
- Apicomplexan and chromerid polyamine queries.

The union is 219 papers. Of these, 79 open-access full texts and 140 abstracts were scanned for fusion sentences and
for sentences about how widespread the arrangement is.

### 2a. Papers on the bifunctional enzyme itself

| Year | Paper | PMID | Full text read? | Does it look outside *Plasmodium*? What it says about distribution |
|---|---|---|---|---|
| 2000 | Müller S et al., *JBC* (discovery) | 10713131 | abstract only | "a unique bifunctional ODC/AdoMetDC in P. falciparum" |
| 2000 | Krause T et al., *Biochem J* (ODC domain) | 11085920 | abstract | none |
| 2001 | Wrenger C et al., *JBC* (domain independence) | 11390378 | abstract | none |
| 2003 | Birkholtz L et al., *Biochem J* (ODC model) | 12557188 | abstract | none |
| 2004 | Birkholtz LM et al., *Biochem J* (inserts, hinge) | 12974675 | abstract | **"The bifunctional nature of AdoMetDC/ODC is unique to Plasmodia and is shared by at least three species."** |
| 2005 | Das Gupta R et al., *Biochem J* (APA inhibitors) | 15980361 | abstract | none |
| 2006 | Wells GA et al., *Biochem J* (AdoMetDC model) | 16257247 | abstract | none |
| 2007 | Williams M et al. (deletion mutagenesis) | 17519001 | full | "unique bifunctional nature" |
| 2008 | Müller IB et al., review | 18455248 | abstract | "Plasmodium has its own and unique polyamine synthesis ... (AdoMetDC/ODC), on a bifunctional protein" |
| 2008 | Clark K et al. (DFMO transcriptome) | 18163886 | abstract | none |
| 2009 | van Brummelen AC et al., *JBC* (co-inhibition) | 19073607 | full (PMC3783055) | none |
| 2010 | Becker JV et al., *BMC Genomics* | 20385001 | full | "a unique bifunctional protein" |
| 2011 | Williams M et al., *Biochem J* (monofunctional PfAdoMetDC) | 21803076 | abstract | "a unique bifunctional ... (PfAdoMetDC/ODC)" |
| 2011 | Birkholtz LM et al., *Biochem J* review (protozoa) | 21834794 | full | "Possibly the most distinctive feature of polyamine biosynthesis in plasmodia is that a single open reading frame encodes both ODC and AdoMetDC in a unique bifunctional protein" |
| 2014 | le Roux D et al. (AdoMetDC inhibitors) | 24596666 | full | "its unique association with ODC in a heterotetrameric bifunctional protein complex" |
| 2014 | Kronenberger T et al. (B6 enzymes review) | 24524072 | full | "P. falciparum possesses a unique polyamine biosynthesis due to the bifunctional organisation" |
| 2015 | Burger PB et al. | 25651815 | full | "unique and exploitable parasite-specific characteristics such as the association of ... AdoMetDC and ODC, in a heterotetrameric bifunctional protein" |
| 2016 | Makhoba XH et al. (chimeric Hsp70) | 27031344 | full | "the unique bifunctional PfAdoMetDC-ODC" |
| 2016 | Hart RJ et al., *Biol Open* (P. yoelii knockout) | 27387533 | full | Fig. 1B: "Schematic representation of conservation of the bifunctional enzyme AdoMetDC/ODC in rodent and human malaria parasite species" (within *Plasmodium* only) |
| 2017 | Meireles P et al. | 28642498 | full | "the two rate-limiting decarboxylase activities are found in a unique protein" |
| 2018 | van Biljon R et al. | 30409996 | full | none |
| 2023, 2026 | Kaiser A; Makhoba XH reviews | 37238673; 42515022 | full | describe it as *Plasmodium*-specific ("unique bifunctional enzyme structures in the malaria parasite") |
| 2024 | Pérez-Pertejo Y et al. review | 38251386 | full | "Plasmodium and other apicomplexans, which contain bifunctional enzymes" |

Reviews that could be read only as abstracts (no fusion statement in the abstract): Phillips MA 2018 *JBC*
"Polyamines in protozoan pathogens" (PMC6290135, not in the open-access subset); Pegg AE 2018 "Unique characteristics
of the parasite polyamine pathway"; Heby et al. 2003; Müller S et al. 2001; Roberts & Ullman 2017; Kappes et al. 2011.

### 2b. Outcome for section 2

Searched, and no *Plasmodium* paper reports the fusion outside apicomplexans. The *Plasmodium* papers repeatedly call
it unique to *Plasmodium*. Birkholtz 2004 confirms it in at least three *Plasmodium* species, and Hart 2016 within
rodent and human malaria species. The only broader statements are Volkov 2016 (9 apicomplexan sequences) and the 2024
review ("other apicomplexans"). Apicomplexan/chromerid-specific searches (Europe PMC 9 and 32 title/abstract hits)
found no paper describing the fusion in *Toxoplasma*, *Cryptosporidium*, *Eimeria*, piroplasms, *Perkinsus* or
chromerids. Woo 2015 is covered in section 7.

---

## 3. Author publication lists

Queries: AUTH:"Michael AJ" (103 hits), AUTH:"Michael A" + polyamine terms (78), Hanfrey C / CC (17 / 8), Fuell C (8),
Elliott KA (158, of which 6 are on polyamines), Phillips MA (406; 85 on polyamines), Ealick SE (230; 41), Pegg AE (560;
333), each repeated in PubMed. AUTH:"Michael AJ" also returns other people with the same initials (ophthalmology,
seismology, cardiology); their titles were read and set aside. The 513 polyamine-related records from all authors
were scanned: 77 as open-access full text, the rest as abstracts.

Flagged as relevant to pathway diversity, unusual architectures or algal/protist enzymes (titles and abstracts read):

| Paper | Readable? | Relevant content |
|---|---|---|
| Illingworth C ... Michael AJ 2003 *FEBS Lett* "diverse bacterial origins of the Arabidopsis polyamine pathway" | abstract | origin of plant ODC/AdoMetDC; no fusion |
| Lee J, Michael AJ ... 2007 *JBC* "Phylogenetic diversity ... beta/alpha-barrel fold basic amino acid decarboxylases" | abstract | ODC-fold diversity; no fusion |
| **Fuell C, Elliott KA, Hanfrey CC, Franceschetti M, Michael AJ 2010** *Plant Physiol Biochem* "Polyamine biosynthetic diversity in plants and algae" | **abstract only (paywalled)** | abstract: "Single-celled green algae have lost the arginine route and are dependent ... on putrescine biosynthesis from the ornithine"; no fusion mentioned in the abstract |
| Shaw FL ... Michael AJ 2010 *JBC* (homospermidine pathway, horizontal transfer) | abstract | no |
| **Michael AJ 2011** *Methods Mol Biol* 720:139 "Exploring polyamine biosynthetic diversity through comparative and functional genomics" | **abstract only (Springer login)** | abstract: "Diversification of polyamine biosynthesis has depended on gene duplication ..., on gene fusion, and on horizontal gene transfer. Four examples ... are presented". The four examples are not named. |
| Michael AJ 2011 *FEBS Lett* (diatom AdoMetDC-APT fusions) | abstract (read earlier) | different fusion |
| Green R ... Michael AJ 2011 *Mol Microbiol* | **full** | see 1a |
| Li B ... Michael AJ 2015 *Mol Microbiol* (ciliates, AdoMetDC-SpdSyn fusion) | abstract | bacterial AdoMetDC-spermidine synthase fusion in *Tetrahymena*; not ODC |
| **Michael AJ 2016** *JBC* (PMC4946907) "Polyamines in Eukaryotes, Bacteria, and Archaea" | **abstract only** (PMC page behind captcha; not in the open-access subset) | abstract does not mention fusions |
| **Michael AJ 2016** *Biochem J* "Biosynthesis of polyamines and polyamine-containing molecules" | **abstract only (HTTP 403)** | abstract: "patterns of evolutionary reiteration ..."; no fusion named |
| **Michael AJ 2017** *Biochem J* "Evolution of biosynthetic diversity" | **abstract only (HTTP 403)** | abstract: "Fusion of different protein domains results in enzymes with novel properties. This review will consider ... gene fusion" (examples not named) |
| Volkov OA ... Phillips MA 2016 *eLife* | **full** | see 1a |
| **Michael AJ 2018** *JBC* "Polyamine function in archaea and bacteria" | abstract only (captcha) | bacteria/archaea only |
| Li B ... Michael AJ 2021 *JBC* (ancestral ODC/LDC) | abstract | no |
| Li B ... Phillips MA, Michael AJ 2023 *JBC* (SpeD neofunctionalization) | **full** | bacterial SpeD-ODC fusions (see 1c) |
| Li B ... Michael AJ 2026 *JBC* (spermine/thermospermine synthases in eukaryotes) | **full** | AdoMetDC-SpmSyn fusion origin of animal spermine synthase; no AdoMetDC-ODC |
| Bale S, Ealick SE 2010 "Structural biology of S-adenosylmethionine decarboxylase" (PMC2847442) | **full** | only "evolved by gene duplication and fusion" (the AdoMetDC fold itself) |
| Pegg AE 2009 "S-Adenosylmethionine decarboxylase" (PMID 20095968); Pegg 2016 *JBC*; Pegg 2018 parasite review | abstract only | no fusion statement in the abstracts |
| Phillips MA 2018 *JBC* "Polyamines in protozoan pathogens" | abstract only | no fusion statement in the abstract |

Outcome: searched. The only statements about the AdoMetDC-ODC fusion in these authors' work are Green 2011
(*Micromonas*) and Volkov 2016 (9 apicomplexan + 3 plant sequences). **The broad diversity reviews could not be
searched in full text**, so whether any of them catalogues a wider distribution is open. Michael 2011 and Michael
2017 are the most likely places, since both abstracts name gene fusion as a theme.

---

## 4. The 23 accessions and their UniProt references

- Each of the 23 accessions as a quoted string: Europe PMC 0 hits, PubMed 0 hits (46 queries).
- 24 locus tags / ORF names and 6 EMBL or RefSeq protein IDs from the UniProt entries, as quoted strings in Europe PMC:
  0 hits each, except **XP_002502359** (*Micromonas commoda* C1E5E7), which gives 1 hit: Green 2011.
- References attached to each UniProt entry (for the 3 deleted entries, from the last UniSave version):

| Accession | Organism | UniProt reference | Read? | Polyamine terms in text |
|---|---|---|---|---|
| A0A9C7PQ51 | *Galdieria partita* | Hirooka S et al. 2022 *PNAS* (PMID 36194630, PMC9565259) | full | none |
| A0A1Y1HLG4 | *Klebsormidium nitens* | Hori K et al. 2014 *Nat Commun* (PMC4052687) | full | none |
| A0AAV8UZF7 | *Rhodosorus marinus* | Lee Y et al. 2023 *Nat Commun* "minicircular mitochondrial genomes" (PMC10250338) | full | none |
| A0AAW1PYA1, A0AAW1RFQ0, A0AAW1T1C3, A0AAW1STH3 | *Myrmecia*, *Apatococcus* ×3 | Puginier C et al. 2024 *Nat Commun* (PMC11126685) | full | only a buffer recipe (spermidine/spermine in nuclei isolation buffer) |
| A0A2P6V5T4 | *Micractinium* | Arriola MB et al. 2018 *Plant J* (PMID 29178410) | **could not search** (abstract only) | - |
| A0A087SFF4 | *Auxenochlorella* | Gao C et al. 2014 *BMC Genomics* (PMC4111847) | full | none |
| A0AAE0FQH6 | *Cymbomonas* | Burns JA et al. 2015 *GBE* (PMC5741210) | full | none |
| C1E5E7 | *Micromonas commoda* | Worden AZ et al. 2009 *Science* (PMID 19359590) | **could not search** (abstract only) | - |
| A0A9P3J9N9 | *Closterium* (deleted) | Sekimoto H et al. 2022 *New Phytol* (PMID 36533897) | **could not search** (not open access) | - |
| A0A2V3IZY9 | *Gracilariopsis* | Lee J et al. 2018 *MBE* (PMID 29688518) | partial (see section 6) | none in the part read |
| A0A830HJR7 | *Pycnococcus* | submission linked to Makita Y et al. 2021 *Nat Commun* (PMC8209157) | full | none |
| A0A9D4ZGS2, A0ABD3ICG1, A0A5B8MZ45, A0AAD9ING6, A0A8T2T6R3, A0A0G4EYD2 | *Adiantum*, *Riccia*, *Chloropicon*, *Prototheca*, *Ceratopteris*, *Vitrella* | database submissions only, no article | - | (genome papers by title checked earlier, in `trackB_fusion_validation.md`) |
| A0A7S1KLJ4, A0A1D2AC22, A0A061RZV3 | *Percolomonas*, *Auxenochlorella*, *Tetraselmis* (all deleted, TSA) | transcriptome submissions only | - | - |

No PDB, KEGG or literature cross-references beyond these exist on any of the 23 entries.

**Correction to `trackB_fusion_validation.md`:** that file says "no *G. partita* genome paper found". UniProt links
Hirooka 2022 *PNAS*. It has now been read in full, and it contains no polyamine terms.

---

## 5. Preprint servers

- Europe PMC preprint index (covers bioRxiv, medRxiv, Research Square, Preprints.org and others; mostly abstract
  level): the section 1 terms with all fusion words OR'd gave 4 hits; without the fusion words, 11 hits. None is
  about an AdoMetDC-ODC fusion (they are parasite/TB polyamine reviews, a tomato symbiont paper, trypanosome cullins,
  ATP13A3).
- **bioRxiv's own full-text search: could not search.** All four search URLs returned HTTP 403. The bioRxiv API
  responds, but it has no keyword search.

---

## 6. The three genome papers previously unread

| Paper | Abstract | Full text | Supplementary data |
|---|---|---|---|
| Lee J, Yang EC, Graf L, et al. (2018) *Mol Biol Evol* 35:1869 "Analysis of the draft genome of the red seaweed *Gracilariopsis chorda* ..." | read; genome size evolution, no polyamine content | publisher page loaded. About 79k characters read verbatim: no polyamine/ODC/AdoMetDC terms. The rest was seen only as a summary or not at all: **partial** | one zip, `msy081_supp.zip` (Figs S1-S16, Tables S1-S24), on the OUP site; **not opened**. Data: GenBank NBIV00000000 (BioProject PRJNA361418) is the source of the gene model already checked in p7. Raw reads at MAGIC Korea (MN00281). |
| Arriola MB et al. (2018) *Plant J* (doi 10.1111/tpj.13789) *Chlorella sorokiniana* / *Micractinium conductrix* genomes | read; maltose excretion focus, no polyamine content | **could not search** (Wiley, HTTP 403) | could not check; the annotation is in GenBank (LHPF02), already used in p7 |
| Worden AZ et al. (2009) *Science* (doi 10.1126/science.1167222) *Micromonas* genomes | read; no polyamine content | **could not search** (Science, HTTP 403) | could not check whether the Supporting Online Material is free. Green 2011 identifies the *Micromonas* fusion from these genome sequences. |

---

## 7. Woo et al. 2015 *eLife* (chromerid genomes), PMC4501334

What was read:
- The main text and appendices (open-access XML, 230k characters).
- Figure 2 and figure supplements 1-5 (images and legends).
- All 16 supplementary/source-data spreadsheets, downloaded from Europe PMC and searched cell by cell.

Findings:
- **Text.** Only two sentences mention polyamines. One reports "loss of enzymes mediating polyamine biosynthesis in
  all lineages except Plasmodium (Figure 2B and Figure 2—figure supplement 3)". The other lists "pathways shared
  between chromerids and a apicomplexan lineage: with Plasmodium, polyamine synthesis". Neither mentions a fusion or
  a bifunctional enzyme.
- **Figure 2B.** Pie charts of pathway completeness per ancestral node, with "Polyamines" as one slice.
- **Figure 2—figure supplement 3** (metabolic map). Box 14 "Polyamine" shows Orn → Put → Spd → Spn arrows coloured
  by phyletic pattern across *Cryptosporidium*, *Babesia*, *Plasmodium*, *Toxoplasma*, *Chromera* and *Vitrella*.
  Arrows stand for reactions; no gene, domain or fusion is drawn.
- **Figure 2—source data 1** (KEGG EC distribution; spreadsheet `elife06974s003.xlsx`, sheet F2-SD1), "Polyamine
  Metabolism, Ornithine → Putrescine → Spermine" block:
  - EC 4.1.1.17 (ODC) is present in VBRA (*Vitrella*) and CVEL (*Chromera*), absent in Atha.
  - EC 4.1.1.50 (AdoMetDC) is present in VBRA, **absent in CVEL**, and also blank for all *Plasmodium* columns
    (Pfal, Pber, Pyoe and others) in that row.
  - EC 2.5.1.16 (spermidine synthase) is present in both.
  - The table lists only EC numbers per species. It has no gene IDs and no domain information.
- No other supplementary table contains polyamine genes. The tables cover genome statistics, orthogroups,
  photosystem, heme, extracellular proteins, invasion genes, RAP domains, flagella and trafficking.

Outcome: searched. **Woo 2015 does not show the *Vitrella* AdoMetDC-ODC fusion.** It shows only that *Vitrella* has
both activities by KEGG EC assignment. One oddity of that table, recorded but not interpreted: the bifunctional
*Plasmodium* enzyme is not assigned EC 4.1.1.50 in it.

---

## What remains unsearched

These items are **could not search**, listed so they can be read with library access:
1. Michael AJ 2011 *Methods Mol Biol* 720:139. Its abstract names gene fusion among four diversification examples.
2. Michael AJ 2017 *Biochem J* 474:2277 "Evolution of biosynthetic diversity". Its abstract names gene fusion.
3. Michael AJ 2016 *Biochem J* 473:2315 and Michael AJ 2016 *JBC* (PMC4946907).
4. Fuell et al. 2010 *Plant Physiol Biochem* (doi 10.1016/j.plaphy.2010.02.008) (plants and algae).
5. Phillips MA 2018 *JBC* "Polyamines in protozoan pathogens"; Pegg AE 2018 "Unique characteristics of the parasite
   polyamine pathway".
6. Worden 2009 *Science* SOM; Arriola 2018 *Plant J*; the Lee 2018 *MBE* supplementary zip and the remainder of its
   main text.
7. bioRxiv full-text search.
8. Volkov 2016 does not list its 3 "plant" sequences. The corresponding author (M. A. Phillips) would know which they
   are.

💡 Claude's suggestion (not Marc's protocol): items 1-4 can likely be opened through the Connecticut College library.
A text search of each for "ornithine", "fusion" and "Micromonas" would settle whether a review catalogued the
distribution.

---

## Query log (all queries run, verbatim, with hit counts)

Also saved as `literature_query_log.tsv`. Section = item number in the request. "Hits" is the count the database
reported on the date shown. Europe PMC "all fields" searches include the open-access full text that Europe PMC
indexes; "TITLE_ABS" restricts to title and abstract. For citation rows, the query is "papers citing PMID X". Not in
the table:
- Web-page fetches: bioRxiv search ×4, HTTP 403. Publisher pages for Lee 2018 (loaded), Arriola 2018 (403), Worden
  2009 (403). Michael/Fuell/Phillips reviews ×6: captcha, 403 or login.
- Full-text sentence scans: 290 + 79 + 77 open-access texts. Woo 2015 and Volkov 2016 supplementary files.

| # | Section | Database | Query (verbatim) | Hits | Date |
|---|---|---|---|---|---|
| 1 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:fusion` | 2 | 2026-11-19 |
| 2 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND fusion` | 218 | 2026-11-19 |
| 3 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:bifunctional` | 12 | 2026-11-19 |
| 4 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND bifunctional` | 99 | 2026-11-19 |
| 5 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:chimeric` | 1 | 2026-11-19 |
| 6 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND chimeric` | 60 | 2026-11-19 |
| 7 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:fused` | 1 | 2026-11-19 |
| 8 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND fused` | 115 | 2026-11-19 |
| 9 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:hybrid` | 8 | 2026-11-19 |
| 10 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND hybrid` | 132 | 2026-11-19 |
| 11 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:multidomain` | 0 | 2026-11-19 |
| 12 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND multidomain` | 8 | 2026-11-19 |
| 13 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:"single polypeptide"` | 1 | 2026-11-19 |
| 14 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND "single polypeptide"` | 9 | 2026-11-19 |
| 15 | 1 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase") AND TITLE_ABS:"gene fusion"` | 0 | 2026-11-19 |
| 16 | 1 | Europe PMC (all fields incl. OA full text) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND "gene fusion"` | 19 | 2026-11-19 |
| 17 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND fusion[tiab]` | 4 | 2026-11-19 |
| 18 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND bifunctional[tiab]` | 11 | 2026-11-19 |
| 19 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND chimeric[tiab]` | 1 | 2026-11-19 |
| 20 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND fused[tiab]` | 1 | 2026-11-19 |
| 21 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND hybrid[tiab]` | 7 | 2026-11-19 |
| 22 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND multidomain[tiab]` | 0 | 2026-11-19 |
| 23 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND "single polypeptide"[tiab]` | 1 | 2026-11-19 |
| 24 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH]) AND "gene fusion"[tiab]` | 0 | 2026-11-19 |
| 25 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab] OR "Ornithine Decarboxylase"[MeSH])` | 796 | 2026-11-19 |
| 26 | 1 | PubMed | `(AdoMetDC[tiab] OR SAMDC[tiab] OR "S-adenosylmethionine decarboxylase"[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "Adenosylmethionine Decarboxylase"[MeSH])` | 1489 | 2026-11-19 |
| 27 | 2 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:Plasmodium AND TITLE_ABS:(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase" OR "S-adenosyl-L-methionine decarboxylase") AND TITLE_ABS:(ODC OR "ornithine decarboxylase")` | 16 | 2026-11-19 |
| 28 | 2 | Europe PMC (all fields incl. OA full text) | `Plasmodium AND (AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase" OR "S-adenosyl-L-methionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND (bifunctional OR fusion OR fused)` | 63 | 2026-11-19 |
| 29 | 2 | Europe PMC (all fields incl. OA full text) | `"AdoMetDC/ODC" OR "PfAdoMetDC/ODC" OR "ODC/AdoMetDC" OR "PfODC/AdoMetDC" OR "AdoMetDC-ODC" OR "SAMDC/ODC" OR "ODC/SAMDC"` | 152 | 2026-11-19 |
| 30 | 2 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(Plasmodium OR malaria) AND TITLE_ABS:(polyamine OR polyamines)` | 166 | 2026-11-19 |
| 31 | 2 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(Apicomplexa OR apicomplexan OR apicomplexans) AND TITLE_ABS:(polyamine OR polyamines)` | 9 | 2026-11-19 |
| 32 | 2 | Europe PMC, title/abstract only (TITLE_ABS) | `TITLE_ABS:(Toxoplasma OR Cryptosporidium OR Eimeria OR Theileria OR Babesia OR Perkinsus OR chromerid OR Chromera OR Vitrella) AND TITLE_ABS:(polyamine OR polyamines OR "ornithine decarboxylase" OR "S-adenosylmethionine decarboxylase" OR AdoMetDC)` | 32 | 2026-11-19 |
| 33 | 2 | Europe PMC citations endpoint | `papers citing PMID 10713131` | 66 | 2026-11-19 |
| 34 | 2 | Europe PMC citations endpoint | `papers citing PMID 11085920` | 23 | 2026-11-19 |
| 35 | 2 | Europe PMC citations endpoint | `papers citing PMID 11390378` | 41 | 2026-11-19 |
| 36 | 2 | Europe PMC citations endpoint | `papers citing PMID 12974675` | 26 | 2026-11-19 |
| 37 | 2 | Europe PMC citations endpoint | `papers citing PMID 27387533` | 18 | 2026-11-19 |
| 38 | 2 | Europe PMC citations endpoint | `papers citing PMID 21762220` | 20 | 2026-11-19 |
| 39 | 2 | Europe PMC citations endpoint | `papers citing PMID 12557188` | 14 | 2026-11-19 |
| 40 | 2 | Europe PMC citations endpoint | `papers citing PMID 16257247` | 13 | 2026-11-19 |
| 41 | 2 | Europe PMC citations endpoint | `papers citing PMID 17519001` | 7 | 2026-11-19 |
| 42 | 2 | Europe PMC citations endpoint | `papers citing PMID 19073607` | 58 | 2026-11-19 |
| 43 | 2 | Europe PMC citations endpoint | `papers citing PMID 21803076` | 7 | 2026-11-19 |
| 44 | 2 | PubMed | `Plasmodium[tiab] AND (AdoMetDC[tiab] OR SAMDC[tiab] OR "adenosylmethionine decarboxylase"[tiab] OR "S-adenosyl-L-methionine decarboxylase"[tiab]) AND (ODC[tiab] OR "ornithine decarboxylase"[tiab])` | 15 | 2026-11-19 |
| 45 | 2 | PubMed | `(Plasmodium[tiab] OR malaria[tiab]) AND (polyamine[tiab] OR polyamines[tiab])` | 154 | 2026-11-19 |
| 46 | 2 | PubMed | `"AdoMetDC/ODC"[tiab] OR "PfAdoMetDC/ODC"[tiab] OR "ODC/AdoMetDC"[tiab]` | 12 | 2026-11-19 |
| 47 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Michael AJ"` | 103 | 2026-11-19 |
| 48 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Michael AJ" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 66 | 2026-11-19 |
| 49 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Michael A"` | 1565 | 2026-11-19 |
| 50 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Michael A" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 78 | 2026-11-19 |
| 51 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Hanfrey C"` | 17 | 2026-11-19 |
| 52 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Hanfrey C" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 14 | 2026-11-19 |
| 53 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Hanfrey CC"` | 8 | 2026-11-19 |
| 54 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Hanfrey CC" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 7 | 2026-11-19 |
| 55 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Fuell C"` | 8 | 2026-11-19 |
| 56 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Fuell C" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 3 | 2026-11-19 |
| 57 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Elliott KA"` | 158 | 2026-11-19 |
| 58 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Elliott KA" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 6 | 2026-11-19 |
| 59 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Phillips MA"` | 406 | 2026-11-19 |
| 60 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Phillips MA" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 85 | 2026-11-19 |
| 61 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Ealick SE"` | 230 | 2026-11-19 |
| 62 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Ealick SE" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 41 | 2026-11-19 |
| 63 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Pegg AE"` | 560 | 2026-11-19 |
| 64 | 3 | Europe PMC (all fields incl. OA full text) | `AUTH:"Pegg AE" AND (polyamine OR polyamines OR spermidine OR spermine OR putrescine OR decarboxylase OR AdoMetDC OR "aminopropyltransferase" OR agmatine OR norspermidine OR "carboxyspermidine" OR hypusine)` | 333 | 2026-11-19 |
| 65 | 3 | PubMed | `Michael AJ[au] AND (polyamine*[tiab] OR spermidine[tiab] OR decarboxylase[tiab] OR putrescine[tiab] OR spermine[tiab])` | 65 | 2026-11-19 |
| 66 | 3 | PubMed | `Pegg AE[au]` | 560 | 2026-11-19 |
| 67 | 3 | PubMed | `Phillips MA[au] AND (polyamine*[tiab] OR decarboxylase[tiab] OR spermidine[tiab] OR AdoMetDC[tiab] OR prozyme[tiab])` | 71 | 2026-11-19 |
| 68 | 3 | PubMed | `Ealick SE[au] AND (decarboxylase[tiab] OR polyamine*[tiab] OR spermidine[tiab] OR AdoMetDC[tiab])` | 28 | 2026-11-19 |
| 69 | 3 | PubMed | `Hanfrey C[au]` | 17 | 2026-11-19 |
| 70 | 3 | PubMed | `Fuell C[au]` | 8 | 2026-11-19 |
| 71 | 3 | PubMed | `Elliott KA[au] AND (polyamine*[tiab] OR spermidine[tiab] OR decarboxylase[tiab])` | 6 | 2026-11-19 |
| 72 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A7S1KLJ4"` | 0 | 2026-11-19 |
| 73 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A9C7PQ51"` | 0 | 2026-11-19 |
| 74 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A1Y1HLG4"` | 0 | 2026-11-19 |
| 75 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAV8UZF7"` | 0 | 2026-11-19 |
| 76 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAW1PYA1"` | 0 | 2026-11-19 |
| 77 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A9D4ZGS2"` | 0 | 2026-11-19 |
| 78 | 4 | Europe PMC (all fields incl. OA full text) | `"A0ABD3ICG1"` | 0 | 2026-11-19 |
| 79 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A2P6V5T4"` | 0 | 2026-11-19 |
| 80 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A087SFF4"` | 0 | 2026-11-19 |
| 81 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAW1RFQ0"` | 0 | 2026-11-19 |
| 82 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A1D2AC22"` | 0 | 2026-11-19 |
| 83 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAE0FQH6"` | 0 | 2026-11-19 |
| 84 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A061RZV3"` | 0 | 2026-11-19 |
| 85 | 4 | Europe PMC (all fields incl. OA full text) | `"C1E5E7"` | 0 | 2026-11-19 |
| 86 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A9P3J9N9"` | 0 | 2026-11-19 |
| 87 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A5B8MZ45"` | 0 | 2026-11-19 |
| 88 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAD9ING6"` | 0 | 2026-11-19 |
| 89 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAW1T1C3"` | 0 | 2026-11-19 |
| 90 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A2V3IZY9"` | 0 | 2026-11-19 |
| 91 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A830HJR7"` | 0 | 2026-11-19 |
| 92 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A8T2T6R3"` | 0 | 2026-11-19 |
| 93 | 4 | Europe PMC (all fields incl. OA full text) | `"A0A0G4EYD2"` | 0 | 2026-11-19 |
| 94 | 4 | Europe PMC (all fields incl. OA full text) | `"A0AAW1STH3"` | 0 | 2026-11-19 |
| 95 | 4 | Europe PMC (all fields incl. OA full text) | `"GpartN1_g196"` | 0 | 2026-11-19 |
| 96 | 4 | Europe PMC (all fields incl. OA full text) | `"KFL_000310060"` | 0 | 2026-11-19 |
| 97 | 4 | Europe PMC (all fields incl. OA full text) | `"NDN08_007560"` | 0 | 2026-11-19 |
| 98 | 4 | Europe PMC (all fields incl. OA full text) | `"WJX72_004436"` | 0 | 2026-11-19 |
| 99 | 4 | Europe PMC (all fields incl. OA full text) | `"GOP47_0010864"` | 0 | 2026-11-19 |
| 100 | 4 | Europe PMC (all fields incl. OA full text) | `"R1sor_018251"` | 0 | 2026-11-19 |
| 101 | 4 | Europe PMC (all fields incl. OA full text) | `"C2E20_7009"` | 0 | 2026-11-19 |
| 102 | 4 | Europe PMC (all fields incl. OA full text) | `"F751_3176"` | 0 | 2026-11-19 |
| 103 | 4 | Europe PMC (all fields incl. OA full text) | `"WJX74_007903"` | 0 | 2026-11-19 |
| 104 | 4 | Europe PMC (all fields incl. OA full text) | `"CYMTET_27057"` | 0 | 2026-11-19 |
| 105 | 4 | Europe PMC (all fields incl. OA full text) | `"MICPUN_58289"` | 0 | 2026-11-19 |
| 106 | 4 | Europe PMC (all fields incl. OA full text) | `"XP_002502359"` | 1 | 2026-11-19 |
| 107 | 4 | Europe PMC (all fields incl. OA full text) | `"CLOM_g8357"` | 0 | 2026-11-19 |
| 108 | 4 | Europe PMC (all fields incl. OA full text) | `"A3770_19p83610"` | 0 | 2026-11-19 |
| 109 | 4 | Europe PMC (all fields incl. OA full text) | `"QBZ16_000624"` | 0 | 2026-11-19 |
| 110 | 4 | Europe PMC (all fields incl. OA full text) | `"WJX84_010671"` | 0 | 2026-11-19 |
| 111 | 4 | Europe PMC (all fields incl. OA full text) | `"BWQ96_02393"` | 0 | 2026-11-19 |
| 112 | 4 | Europe PMC (all fields incl. OA full text) | `"PPROV_000455800"` | 0 | 2026-11-19 |
| 113 | 4 | Europe PMC (all fields incl. OA full text) | `"KP509_15G068700"` | 0 | 2026-11-19 |
| 114 | 4 | Europe PMC (all fields incl. OA full text) | `"Vbra_8476"` | 0 | 2026-11-19 |
| 115 | 4 | Europe PMC (all fields incl. OA full text) | `"WJX84_010316"` | 0 | 2026-11-19 |
| 116 | 4 | Europe PMC (all fields incl. OA full text) | `"TSPGSL018_17490"` | 0 | 2026-11-19 |
| 117 | 4 | Europe PMC (all fields incl. OA full text) | `"TSPGSL018_21677"` | 0 | 2026-11-19 |
| 118 | 4 | Europe PMC (all fields incl. OA full text) | `"PCOS0759_LOCUS847"` | 0 | 2026-11-19 |
| 119 | 4 | Europe PMC (all fields incl. OA full text) | `"CAD9077615"` | 0 | 2026-11-19 |
| 120 | 4 | Europe PMC (all fields incl. OA full text) | `"JAT76708"` | 0 | 2026-11-19 |
| 121 | 4 | Europe PMC (all fields incl. OA full text) | `"JAC77513"` | 0 | 2026-11-19 |
| 122 | 4 | Europe PMC (all fields incl. OA full text) | `"JAC75932"` | 0 | 2026-11-19 |
| 123 | 4 | Europe PMC (all fields incl. OA full text) | `"GJP49097"` | 0 | 2026-11-19 |
| 124 | 4 | Europe PMC (all fields incl. OA full text) | `"XP_011397346"` | 0 | 2026-11-19 |
| 125 | 4 | Europe PMC (title / ID lookup) | `EXT_ID:36194630 OR EXT_ID:38789482 OR EXT_ID:25012212 OR EXT_ID:26224703 OR EXT_ID:36533897 OR EXT_ID:29178410 OR EXT_ID:24865297 OR EXT_ID:19359590 OR EXT_ID:29688518` | 9 | 2026-11-19 |
| 126 | 4 | Europe PMC (title / ID lookup) | `TITLE:"Origin of minicircular mitochondrial genomes in red algae"` | 1 | 2026-11-19 |
| 127 | 4 | Europe PMC (title / ID lookup) | `TITLE:"Dualchrome1"` | 1 | 2026-11-19 |
| 128 | 4 | Europe PMC (title / ID lookup) | `TITLE:"transcriptome of the halophilic microalga Tetraselmis"` | 0 | 2026-11-19 |
| 129 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Chloropicon AND genome)` | 1 | 2026-11-19 |
| 130 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Adiantum AND genome)` | 7 | 2026-11-19 |
| 131 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Ceratopteris AND genome)` | 5 | 2026-11-19 |
| 132 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Riccia AND genome)` | 3 | 2026-11-19 |
| 133 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Vitrella AND genome)` | 0 | 2026-11-19 |
| 134 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Prototheca AND wickerhamii AND genome)` | 8 | 2026-11-19 |
| 135 | 4 | Europe PMC (all fields incl. OA full text) | `AUTH:"Makita Y" AND (Pycnococcus OR photoreceptor) AND TITLE:(photoreceptor)` | 3 | 2026-11-19 |
| 136 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Pycnococcus)` | 6 | 2026-11-19 |
| 137 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Tetraselmis AND transcriptome)` | 0 | 2026-11-19 |
| 138 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Percolomonas)` | 4 | 2026-11-19 |
| 139 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Gracilariopsis AND genome)` | 8 | 2026-11-19 |
| 140 | 4 | Europe PMC (title / ID lookup) | `TITLE:(Micractinium AND genome)` | 7 | 2026-11-19 |
| 141 | 4 | PubMed | `A0A7S1KLJ4[All Fields]` | 0 | 2026-11-19 |
| 142 | 4 | PubMed | `A0A9C7PQ51[All Fields]` | 0 | 2026-11-19 |
| 143 | 4 | PubMed | `A0A1Y1HLG4[All Fields]` | 0 | 2026-11-19 |
| 144 | 4 | PubMed | `A0AAV8UZF7[All Fields]` | 0 | 2026-11-19 |
| 145 | 4 | PubMed | `A0AAW1PYA1[All Fields]` | 0 | 2026-11-19 |
| 146 | 4 | PubMed | `A0A9D4ZGS2[All Fields]` | 0 | 2026-11-19 |
| 147 | 4 | PubMed | `A0ABD3ICG1[All Fields]` | 0 | 2026-11-19 |
| 148 | 4 | PubMed | `A0A2P6V5T4[All Fields]` | 0 | 2026-11-19 |
| 149 | 4 | PubMed | `A0A087SFF4[All Fields]` | 0 | 2026-11-19 |
| 150 | 4 | PubMed | `A0AAW1RFQ0[All Fields]` | 0 | 2026-11-19 |
| 151 | 4 | PubMed | `A0A1D2AC22[All Fields]` | 0 | 2026-11-19 |
| 152 | 4 | PubMed | `A0AAE0FQH6[All Fields]` | 0 | 2026-11-19 |
| 153 | 4 | PubMed | `A0A061RZV3[All Fields]` | 0 | 2026-11-19 |
| 154 | 4 | PubMed | `C1E5E7[All Fields]` | 0 | 2026-11-19 |
| 155 | 4 | PubMed | `A0A9P3J9N9[All Fields]` | 0 | 2026-11-19 |
| 156 | 4 | PubMed | `A0A5B8MZ45[All Fields]` | 0 | 2026-11-19 |
| 157 | 4 | PubMed | `A0AAD9ING6[All Fields]` | 0 | 2026-11-19 |
| 158 | 4 | PubMed | `A0AAW1T1C3[All Fields]` | 0 | 2026-11-19 |
| 159 | 4 | PubMed | `A0A2V3IZY9[All Fields]` | 0 | 2026-11-19 |
| 160 | 4 | PubMed | `A0A830HJR7[All Fields]` | 0 | 2026-11-19 |
| 161 | 4 | PubMed | `A0A8T2T6R3[All Fields]` | 0 | 2026-11-19 |
| 162 | 4 | PubMed | `A0A0G4EYD2[All Fields]` | 0 | 2026-11-19 |
| 163 | 4 | PubMed | `A0AAW1STH3[All Fields]` | 0 | 2026-11-19 |
| 164 | 5 | Europe PMC, preprints only (SRC:PPR) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase") AND (fusion OR bifunctional OR chimeric OR fused OR hybrid OR multidomain OR "single polypeptide" OR "gene fusion")` | 4 | 2026-11-19 |
| 165 | 5 | Europe PMC, preprints only (SRC:PPR) | `(AdoMetDC OR SAMDC OR "S-adenosylmethionine decarboxylase" OR "adenosylmethionine decarboxylase") AND (ODC OR "ornithine decarboxylase")` | 11 | 2026-11-19 |
