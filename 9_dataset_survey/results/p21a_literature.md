# p21 A: Green 2011 and Li 2023, re-read

Both papers were read in full text for this note (Green 2011 as the NIH author manuscript from NCBI PMC efetch;
Li 2023 from Europe PMC full-text XML). Quotations are verbatim.

## Verdict

- **Green 2011 is one paper, and both descriptions are right.** It reports two different fusions: a eukaryotic
  AdoMetDC-ODC fusion in *Micromonas* (sequence only), and a bacterial ODC-AdoMetDC fusion in the actinobacterium
  *Streptomyces* sp. AA4. That bacterial organism is in our own p10 set, so Green 2011 already reported the type of
  fusion that makes up most of our bacterial rows.
- **Li 2023 does describe AdoMetDC fusions whose ODC part has lost activity and pyruvoyl formation.** The family is
  different from ours (a pyruvoyl-dependent, SpeD-homologous ODC, not the PLP fold-III ODC), but it is direct
  precedent for "a fused ODC domain degrades and is kept". Our novelty claim has to be narrowed accordingly.

## 1. Green et al. 2011

**Citation.** Green R, Hanfrey CC, Elliott KA, McCloskey DE, Wang X, Kanugula S, Pegg AE, Michael AJ (2011).
"Independent evolutionary origins of functional polyamine biosynthetic enzyme fusions catalyzing de novo diamine to
triamine formation." *Mol Microbiol* 81(4):1109-1124. PMID 21762220,
PMC3196669, doi:10.1111/j.1365-2958.2011.07757.x.

One paper. The question "bacterial fusion or *Micromonas* fusion?" has the answer "both":

**(a) The *Micromonas* fusion (eukaryotic, AdoMetDC N-terminal, ODC C-terminal; identified from sequence only).**

> "In eukaryotes, we found a gene fusion of a eukaryotic AdoMetDC and an ornithine decarboxylase in two species of
> single-celled prasinophyte green algae, Micromonas pusilla CCMP1545 (EEH58717, 1019 a.a.) and Micromonas sp.
> RCC299 (XP_002502359, 1039 a.a.). The AdoMetDC domain is N-terminal and the ODC domain C-terminal [...]"

> "Although the Micromonas AdoMetDC-ODC gene fusions are arranged similarly to the Plasmodium AdoMetDC-ODC fusion
> genes (Muller et al., 2000), these fusion genes represent entirely independent eukaryotic gene fusion events
> because both domains of the Micromonas fusion protein exhibit greatest homology with the individual component genes
> from green algae and plants and present only distant similarity to the Plasmodium sequences."

XP_002502359 is *Micromonas commoda* RCC299 = our C1E5E7. The paper gives no activity data for the *Micromonas*
protein; the abstract says these fusions were "identified".

**(b) A bacterial ODC-AdoMetDC fusion (ODC N-terminal, prokaryotic AdoMetDC C-terminal).**

> "In the genome of the actinobacterium Streptomyces sp. AA4, the ornithine decarboxylase ORF has fused with the
> immediately downstream AdoMetDC ORF to form a 498 aa ODC-AdoMetDC fusion protein while maintaining a separate
> downstream aminopropyltransferase ORF (Fig. 8A)."

> "[...] a cluster of unfused ODC/AdoMetDC/aminopropyltransferase ORFs has given rise to ODC-AdoMetDC,
> AdoMetDC-aminopropyltransferase, and ODC-aminopropyltransferase fusions."

*Streptomyces* sp. AA4 is UniProt taxon 591158, and it is in `results/p10_fusion_search.tsv` with two sequences
(UPI0018F8B29C, 464 aa; UPI00256FE547, 474 aa), both PF02784 first and PF02675 (prokaryotic AdoMetDC) after it.
That is the architecture of 64 of the 66 bacterial rows (p20). The paper does not report activity for this fusion
either; its functional work is on the AdoMetDC-aminopropyltransferase fusions.

**What was wrong in our files.** Nothing is contradictory: `literature_adometdc_odc_fusion.md` already lists Green
2011 under both the *Micromonas* fusion (section 1a) and the bacterial fusion (1c). The phrase "Green 2011 describes
a bacterial ODC-AdoMetDC fusion" is correct but incomplete if read alone; "the Micromonas fusion report" is also
correct but incomplete. What the files miss is that the Green 2011 bacterial organism is in our dataset.

## 2. Li et al. 2023, *J Biol Chem* (PMC10407285)

**Citation.** Li B, Liang J, Phillips MA, Michael AJ (2023). "Neofunctionalization of S-adenosylmethionine
decarboxylase into pyruvoyl-dependent L-ornithine and L-arginine decarboxylases is widespread in bacteria and
archaea." *J Biol Chem* 299(8):105005. PMID 37399976, PMC10407285, doi:10.1016/j.jbc.2023.105005 (full text read).

**Does it describe AdoMetDC fusions in which ODC activity or pyruvoyl formation has been lost? Yes.**

The fusions are SpeD1-SpeD2 (AdoMetDC/SpeD fused to a SpeD-homologous, pyruvoyl-dependent ODC, "SpeDC"). Some are
active for both reactions; others have a degraded ODC domain:

> "We identified fusion proteins of bona fide AdoMetDC/SpeD with homologous L-ornithine decarboxylases that possess
> two, unprecedented internal protein-derived pyruvoyl cofactors."

> "However, many of the Lentisphaerota fusion proteins lack key active site residues in the C-terminal SpeD domain."

> "Omnitrophica SpeD1SpeD2 fusion proteins are fully processed into the predicted α-, β- and γ-subunits, whereas the
> V. vadensis SpeD1/degradedSpeD2 fusion protein is processed into only two subunits."

> "The homologous fusion protein encoded by Lentisphaerota bacterium V. vadensis has lost critical active site
> residues in the C-terminal SpeDC domain and consequently cannot self-cleave and so lacks ODC activity."

> "We detected additional SpeD-SpeDC fusion proteins with a degraded C-terminal domain in the sulfate-reducing
> Desulfobacteria class of the Thermodesulfobacteriota phylum, for example, Desulfobacter latus (WP_178365596) [...]"

> "Fusion of speD and speDC genes resulted in a bifunctional enzyme; subsequently SpeDC pyruvoyl formation and ODC
> activity was lost in some of the Lentisphaerota-related lineage."

And they use this as the model for the origin of the eukaryotic AdoMetDC, which is itself a fusion keeping a dead
second domain:

> "It is more plausible that the original fusion protein was a SpeD-SpeDC fusion protein with defunct SpeDC domain
> like the V. vadensis fusion protein."

> "[...] the eukaryotic AdoMetDC structure is composed of two class 1b polypeptides fused together, with the
> C-terminal domain lacking a self-cleavage site, and consequently lacking a C-terminal pyruvoyl cofactor (12)."

**How this differs from our case.** Their degraded domain is a pyruvoyl-dependent SpeD-fold ODC; ours is the PLP
fold-III ODC (PF02784/PF00278) fused to the eukaryotic class-2 AdoMetDC (PF01536). Their organisms are bacteria;
ours are red algae and *Percolomonas*. Li 2023 does not mention PLP-ODC fusions or red algae.

## What this changes

- **Novelty claim.** "A fused ODC domain that has lost its catalytic residues but is retained" is not new as a
  phenomenon: Li 2023 shows it biochemically (no self-cleavage, no ODC activity) in bacterial AdoMetDC-ODC fusions,
  and argues the eukaryotic AdoMetDC descends from such a protein. What may still be new is the same pattern in the
  PLP-dependent ODC fused to eukaryotic AdoMetDC, in red algae and *Percolomonas*. Li 2023 should be cited as the
  precedent, and the claim worded as "a second, independent instance in a different decarboxylase family".
- **Li 2023 also gives a model to test against.** Their degraded fusions occur in genomes with no other ODC route
  ("all genomes encoding an anomalous AdoMetDC homolog also lacked canonical ADC (SpeA) or ODC (SpeC)" refers to the
  neofunctionalized SpeD homologs, not to the degraded fusions; the paper does not say whether the degraded-fusion
  genomes have another ODC). p19 found that most of our degraded-fusion species do carry an intact standalone ODC.
- **Bacterial set.** Green 2011 already reported the bacterial ODC-AdoMetDC fusion type, so that set adds organisms,
  not a new fusion type.

Source texts were saved to the session scratchpad only, not to the repository.
