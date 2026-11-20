# p11: AlphaFold Server homodimer jobs (the dimer question)

**Hypothesis.** One dimerisation of a fusion gives both halves their partner at once: AdoMetDC with AdoMetDC, and ODC with ODC. A dead ODC half may then still carry weight as an interface.

Each FASTA holds two identical chains (A, B). Paste both into one AlphaFold Server job as two protein chains, or set one chain with copies = 2.

## Jobs

Each job is two copies of one protein, so the job size is twice the protein length.

| job | protein | aa per chain | job size | role |
|---|---|---|---|---|
| 01 | A0A1Y1HLG4 *Klebsormidium nitens* | 864 | 1,728 | fusion, ODC half intact (13/13) |
| 02 | A0A2V3IZY9 *Gracilariopsis chorda* | 730 | 1,460 | fusion, ODC half degraded; no other ODC with a catalytic Cys found in this species |
| 03 | A0A9C7PQ51 *Galdieria partita* | 762 | 1,524 | fusion, ODC half degraded; the species has standalone ODCs (the contrast) |
| 04 | C1E5E7 *Micromonas commoda* | 1,039 | 2,078 | fusion, the one already reported (Green 2011) |
| 05 | Q8IJ77 *Plasmodium falciparum* | 1,434 | 2,868 | the known bifunctional enzyme |
| 06 | A0AAV8UQI0 *Rhodosorus marinus* | 456 | 912 | control: standalone ODC, 13/13 sites |
| 07 | A0A9D4U7K2 *Adiantum capillus-veneris* | 380 | 760 | control: standalone AdoMetDC |

**Size limit.** As far as I know, AlphaFold Server allows 5,000 tokens per job, and one protein residue is one token. Every job is under that; the largest is 05 at 2,868. Check the Server's current limit when you submit.

**Why job 07 is *Adiantum*.** A0A9D4U7K2 is a standalone AdoMetDC from the dataset (InterPro source, model pLDDT 88, all of Marc's sites kept). *Adiantum* also carries a fusion, A0A9D4ZGS2, so the control and that fusion come from the same species.

## Residue map

Residues were placed by local alignment (BLOSUM62) to human AdoMetDC (P17707) and to human ODC (P11926). For the four fusions, the ODC positions agree with p7.

| job | AdoMetDC: S68 / C82 / S229 / H243 / E247 | AdoMetDC Pfam | ODC Pfam | ODC K69 | ODC C360 |
|---|---|---|---|---|---|
| 01 | S67 / C81 / S232 / H245 / E249 | 7–334 | 440–776 | K463 | C749 |
| 02 | S86 / C100 / S261 / H274 / E278 | 26–318 | 399–706 | **R429** | C680 (3.5 Å off in p7) |
| 03 | S78 / C92 / S238 / H251 / E255 | 19–319 | 434–512 (Pfam match is short) | **R441** | **G710** |
| 04 | S162 / C176 / S352 / H365 / E369 | 105–453 | 634–1003 | K656 | C976 |
| 05 | S73 / C87 / not placed (parasite insertions break the local alignment) | – | – | K880 | C1370 |
| 06 | – | – | whole chain | K115 | C403 |
| 07 | S67 / C81 / S230 / H243 / E247 | 5–330 | – | – | – |

## What to look for

**1. Which domain pairs with which across the two chains.** For each model, check four pairings:

- AdoMetDC(A) with AdoMetDC(B)
- ODC(A) with ODC(B)
- AdoMetDC(A) with ODC(B), and the reverse

The hypothesis predicts the first two together in one model. If you only see AdoMetDC–ODC contacts across chains, or an interface in one domain only, the hypothesis is not supported.

**2. Whether the ODC halves form the canonical ODC active-site interface.**

- In ODC, each active site sits at the dimer interface: Lys69 (the PLP Schiff base) comes from one subunit and Cys360 from the other.
- In a model, measure Lys69(A) to Cys360(B) and Lys69(B) to Cys360(A), using the residue numbers in the map above.
- Don't judge "close enough" against a number from memory. Superpose the model's ODC–ODC pair onto a human ODC dimer crystal structure and measure the same pair there; the canonical arrangement should give a similar distance.
- Do the same with job 06 as the positive control.
- For 02 and 03 the residues at those positions are R/C and R/G, so the question is only whether the positions line up across the interface. The residues themselves are already known to be substituted.

**3. Confidence of the interface itself, not just of the model.**

- Use the inter-chain PAE block restricted to the ODC(A)–ODC(B) and AdoMetDC(A)–AdoMetDC(B) residue ranges. Low PAE there means AlphaFold is confident about the relative placement.
- ipTM is a single number over all chains, so in a two-domain chain it mixes both interfaces. My understanding of AlphaFold's guidance is that above about 0.8 is confident and below about 0.6 has likely failed. Read it alongside the per-domain PAE, not alone.
- AlphaFold Server returns 5 models. Check whether all 5 give the same pairing.

**4. Whether the degraded ODC halves are still interface-competent.** Compare across three groups:

| group | jobs |
|---|---|
| degraded | 02, 03 |
| intact | 01, 04 |
| standalone ODC | 06 |

What to compare:

- whether an ODC–ODC interface forms at all
- its PAE
- its interface pLDDT
- roughly how many residues are in contact

A degraded half that still forms a confident canonical ODC dimer is consistent with a structural (load-bearing) role. One that does not dimerise, or only dimerises with low confidence, is not.

**5. The AdoMetDC side.**

- Compare the fusions' AdoMetDC–AdoMetDC contacts with job 07.
- If 07 comes back as two chains with no confident interface, the plant/algal AdoMetDC may not dimerise on its own. Then a fusion that does pair its AdoMetDC halves would be getting that from the fusion. My recollection, which I have not checked in this session, is that the potato AdoMetDC crystal structure is a monomer, unlike the human (αβ)2 dimer. Verify this before relying on it.

**6. The Plasmodium job (05) is a reference for topology, not a clean positive control.**

- *P. falciparum* AdoMetDC/ODC is described as a heterotetramer, two of each domain.
- The parasite insertions will likely model with low pLDDT. Judge the interface on the core domains only.

## Caveats

- **AlphaFold models proenzymes.** The real AdoMetDC cuts itself at Ser68 into β and α chains and forms the pyruvoyl group. PLP and the pyruvoyl group are absent unless you add ligands, and I am not sure PLP is on the Server's allowed ligand list.
- **A confident homodimer prediction is not evidence the dimer exists.** AlphaFold readily reproduces the dimer seen in homologous crystal structures. This matters especially for 02 and 03, whose ODC halves are still recognisably ODC-fold. The comparison across jobs (degraded vs intact vs control) is more informative than any single model.
- **Analysis tool.** Once the models are back, a small script could measure the four domain–domain contacts, the K69–C360 cross-distances and the per-block PAE from the Server's output files. Not written yet.
