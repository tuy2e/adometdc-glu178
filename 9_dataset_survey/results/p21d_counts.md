# p21 D: the headline counts, restated

Source: `results/p10_fusion_search.tsv` (the 311 UniParc sequences with an AdoMetDC and a PLP-decarboxylase Pfam
feature, which `p10_all_fusions.tsv` also lists). Genus ranks from UniProt taxonomy (all 199 taxids resolved;
cached in `data/p21/taxonomy_ranks.json`, genus membership in `data/p21/genera.json`).

## Organisms and genera

| | count | note |
|---|---|---|
| sequences | 311 | UniParc sequences, not proteins of distinct organisms |
| unique organisms | **199** | 199 organism names = 199 taxids. Your 199 is right. |
| unique genera, first word of the name | 83 | your 83 is right for this method, but it counts non-genera ("uncultured", "unclassified", "bacterium", the class "Actinomycetes") and splits nothing |
| unique genera, UniProt rank "genus" | **79** | plus 7 taxa with no genus rank: uncultured *Pyramimonas*, uncultured *Micromonas*, uncultured *Mantoniella*, uncultured Pyramimonadales, uncultured alveolate, Actinomycetes, "bacterium" |

## Two different fusions: count them separately

| set | sequences | organisms (taxids) | genera (rank genus) | domain order |
|---|---|---|---|---|
| eukaryotic | **245** | 163 | 70 | AdoMetDC (PF01536, eukaryotic class 2) before the PLP domain in all 245: 239 start with PF01536; 4 start with the AdoMetDC leader PF08132 (*Closterium* x3, *Marchantia*); 2 *Chlorella* sequences carry trehalose-6-phosphate synthase/phosphatase domains (PF00982, PF02358) before the AdoMetDC, possibly merged gene models |
| bacterial | **66** | 36 | 9 (+6 sequences with no genus rank) | 64 have the PLP decarboxylase FIRST (PF02784) and a prokaryotic AdoMetDC (PF02675) after it; 1 (*Tistrella*, A0A2E3ETT7) has PF02675 first; 1 ("bacterium" UPI0023698E55) has the eukaryotic PF01536 first |

- The 64 PLP-first bacterial sequences are the fusion Green et al. 2011 already reported (ODC-AdoMetDC in the
  actinobacterium *Streptomyces* sp. AA4, which is itself in this set; see `p21a_literature.md`). Their
  decarboxylase is nearest the *Selenomonas* LDC/ODC (p10/p20). Genera: *Amycolatopsis* 33 sequences,
  *Kibdelosporangium* 17, *Streptomyces* 3, *Actinosynnema* 2, one each *Actinocrispum*, *Actinophytocola*,
  *Pseudonocardia*, *Rhodococcus*; all Actinomycetota apart from *Tistrella*.
- UPI0023698E55 ("bacterium", taxid 1869227, 602 aa) has a eukaryotic class-2 AdoMetDC domain and a PLP domain whose
  nearest reference is a plant LDC2 (50 % identity, p10). That combination looks eukaryotic; it is probably a
  eukaryotic sequence placed in a bacterial bin. Not tested further.
- Headline wording: "245 eukaryotic sequences from 163 organisms in 70 genera (244 / 162 without the primate row
  below), and separately 66 bacterial sequences of a different fusion type (PLP domain first, prokaryotic
  AdoMetDC) from 36 organisms, mostly *Amycolatopsis* and *Kibdelosporangium*".

## Plasmodium

65 sequences = 47 organism names (strains and isolates) of **23 species**. The largest: *P. falciparum* 16
sequences, *P. vivax* 10, *P. ovale* 5, *P. vinckei* 5, *P. chabaudi* 4, *P. berghei* 3, *P. reichenowi* 3;
16 species have one or two. *Plasmodium* is 65 of the 245 eukaryotic sequences (27 %) but 23 of the 163 eukaryotic
organisms.

## The Metazoa row

UPI00108F000D = UniProt A0A8C9GKE0, *Piliocolobus tephrosceles* (Ugandan red colobus monkey), 1,532 aa, PF01536
9-269 > PF01536 713-818 > PF02784 1109-1383. It is an unreviewed entry from the Ensembl annotation of the host
genome (Ensembl ENSPTET00000003839.1, proteome UP000694416).

**Call: contamination by a haemosporidian parasite, not a primate gene.** Evidence:
- Its best match among all 311 is the fusion of ***Hepatocystis* sp. ex *Piliocolobus tephrosceles*** (UPI0012E4CB6E):
  68.6 % identity over 1,452 aligned residues. The next matches are all *Plasmodium* (*P. gaboni*, *P. billcollinsi*,
  *P. adleri*, about 54 %).
- *Hepatocystis* is a malaria-related parasite of this same monkey species; a genome assembled from an infected
  animal can contain parasite contigs, which the host annotation then calls as host genes.
- Its architecture (AdoMetDC split by a long insertion, then ODC, about 1,500 aa) is the haemosporidian one.
- No other metazoan in the set has the fusion.

What is not settled: at 68.6 % it is not the sequenced *Hepatocystis* isolate itself, so the exact source parasite
(another *Hepatocystis* lineage or a co-infecting species) is unknown. It should be removed from organism counts
and labelled as parasite-derived, not counted as a metazoan.
