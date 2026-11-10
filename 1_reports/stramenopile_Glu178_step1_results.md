# The putrescine "sensor" residue is missing in every stramenopile AdoMetDC — step-1 results
*Tuy, 5 Nov 2026. 💡 My computational analysis of **your own FOLDTRACE run** (read-only copy; nothing in your project folder was changed). Not Marc's protocol. No wet lab involved.*

## Short version
When I went back to the raw FOLDTRACE data behind the "diatom + brown algae" WOW finding, the pattern turned out to be **bigger and cleaner than we thought**. It isn't just diatoms and brown algae. **Every stramenopile AdoMetDC in the run, 31 of 31, has lost Glu178.** Glu178 is one of the three acidic residues that grip putrescine in the human enzyme. In all other eukaryotes it's kept **509 of 527** times (100% of fungi, ~99% of plants and animals).

That shifts the story from "two algal groups" to **"a whole major branch of the eukaryote tree replaced the putrescine sensor"**. It's a clean, FOLDTRACE-derived, computational result.

---

## What I did (all computational, all on your data)
1. **Copied (read-only)** your WOW-run files from `8_extra_runs/`: pocket calls for all 1,000 candidates and the taxonomy table.
2. Rebuilt the pocket table for the **617 fold-passing proteins** at the human positions L13, F111, D174, T176, E178, E256 and F285. I grouped them by NCBI class/phylum.
3. Counted Glu178 status for **stramenopiles vs every other eukaryote**.
4. **Checked the result with a second, independent method.** I pulled the *Phaeodactylum* protein (B7FPB8) from UniProt and confirmed the 5 residues FOLDTRACE reported (T12, G124, D194, N198, M318). Then I aligned it to human AdoMetDC by **sequence alone** (MAFFT). The sequence alignment matches the structural mapping **at all 7 pocket positions**: human E178 ↔ diatom N198 and D174 ↔ D194, and the local stretch aligns without gaps:
   - human `…PDQTLEILM…`
   - diatom `…TDRVLNIMM…`
   
   Two independent methods give the same answer.
5. Ran a **novelty check** on PubMed (below).

## The result
| | Glu178 kept | Glu178 replaced |
|---|---|---|
| **Stramenopiles** (21 Ochrophyta from 11 classes, 4 Bigyra, 6 Oomycota) | **0** | **31** |
| All other eukaryotes (animals, fungi, plants, green/red algae, haptophytes, other protists) | 509 | 18 |

**What replaced it:**
- **Ochrophyta** (photosynthetic stramenopiles): **Asn** in 20 of 21 and Asp in 1. All 21 also carry **Gly at F111** (100%), and most carry Thr at L13.
- **Bigyra** (*Blastocystis*, *Cafeteria*, *Hondaea*): **Asn** in 4 of 4.
- **Oomycota** (water moulds, *Phytophthora*): **His** in 4 and Gln in 2.

The rest of the pocket is mostly kept. **Glu256 is kept in 31 of 31.** Asp174 is kept in most Ochrophyta and is Glu in oomycetes. So this is a *specific swap of the sensor residue*, not a pocket that has fallen apart.

**Simplest reading:** an Asn at 178 arose **once, near the base of the stramenopiles**. Oomycetes later changed it to His or Gln. Gly111 was added in Ochrophyta. *This is a hypothesis; the tree analysis below would test it.*

**Other groups that lost Glu178 independently (convergence):** Euglenozoa (the trypanosomatids have Ser; *Bodo*, *Neobodo* and *Eutreptiella* have Asn), cryptophytes (*Guillardia* and *Hemiselmis* have Asn, plus a built-in Lys13), a few cestode tapeworms (Gln/Lys, still to be checked), and 5 bacteria (Gln). The same residue has been replaced several separate times.

### Full stramenopile table (FOLDTRACE calls, human numbering)
| Clade | Class | Organism | UniProt | L13 | F111 | D174 | **E178** | E256 | F285 |
|---|---|---|---|---|---|---|---|---|---|
| Ochrophyta | Bacillariophyceae | *Cylindrotheca closterium* | A0AAD2GAL3 | T | G | D | **N** | E | M |
| Ochrophyta | Bacillariophyceae | *Fragilariopsis cylindrus* | A0A1E7FKD8 | T | G | D | **N** | E | M |
| Ochrophyta | Bacillariophyceae | *Nitzschia inconspicua* | A0A9K3LMG5 | T | G | D | **N** | E | M |
| Ochrophyta | Bacillariophyceae | *Phaeodactylum tricornutum* | B7FPB8 | T | G | D | **N** | E | M |
| Ochrophyta | Bolidophyceae | *Triparma laevis* | A0A9W7AQ06 | - | G | D | **N** | E | M |
| Ochrophyta | Bolidophyceae | *Triparma strigata* | A0A9W7EFJ6 | T | G | D | **N** | E | M |
| Ochrophyta | Chrysophyceae | *Spumella elongata* | A0A7S3H9T9 | - | G | Y | **N** | E | F |
| Ochrophyta | Coscinodiscophyceae | *Chaetoceros tenuissimus* | A0AAD3D8K8 | T | G | D | **N** | E | M |
| Ochrophyta | Coscinodiscophyceae | *Cyclotella cryptica* | A0ABD3QT96 | T | G | D | **N** | E | M |
| Ochrophyta | Coscinodiscophyceae | *Proboscia inermis* | A0A7S0C0C1 | L | G | D | **N** | E | M |
| Ochrophyta | Dictyochophyceae | *Octactis speculum* | A0A7S2HLL7 | N | G | T | **N** | E | M |
| Ochrophyta | Dictyochophyceae | *Octactis speculum* | A0A7S2HLM3 | N | G | T | **N** | E | M |
| Ochrophyta | Eustigmatophyceae | *Nannochloropsis salina* | A0A4D9CSB4 | T | G | D | **N** | E | M |
| Ochrophyta | Fragilariophyceae | *Cyclophora tenuis* | A0A7S1CXK8 | - | G | D | **D** | E | M |
| Ochrophyta | Pelagophyceae | *Aureococcus anophagefferens* | F0Y0A9 | T | G | E | **N** | E | F |
| Ochrophyta | Pelagophyceae | *Chrysophaeum taylorii* | A0AAD7XHA9 | T | G | E | **N** | E | F |
| Ochrophyta | Phaeophyceae | *Ectocarpus siliculosus* | D7FY29 | T | G | D | **N** | E | F |
| Ochrophyta | Phaeophyceae | *Ectocarpus sp.* | A0A6H5KT37 | T | G | D | **N** | E | F |
| Ochrophyta | Phaeophyceae | *Ectocarpus sp.* | A0A6H5KVY3 | T | G | D | **N** | E | V |
| Ochrophyta | Pinguiophyceae | *Phaeomonas parva* | A0A7S1XYU7 | T | G | D | **N** | E | F |
| Ochrophyta | Xanthophyceae | *Tribonema minus* | A0A835ZL65 | T | G | D | **N** | E | M |
| Bigyra | Bigyra | *Blastocystis hominis* | D8LUU6 | I | L | S | **N** | E | T |
| Bigyra | Bigyra | *Blastocystis sp.* | A0A196S442 | I | L | G | **N** | E | C |
| Bigyra | Bigyra | *Cafeteria roenbergensis* | A0A5A8CY83 | T | C | D | **N** | E | F |
| Bigyra | Bigyra | *Hondaea fermentalgiana* | A0A2R5GFL0 | I | G | E | **N** | E | F |
| Oomycota | Oomycota | *Aphanomyces invadans* | A0A024TG05 | N | R | E | **H** | E | S |
| Oomycota | Oomycota | *Hyaloperonospora arabidopsidis* | M4BIM5 | N | N | E | **H** | E | Y |
| Oomycota | Oomycota | *Lagenidium giganteum* | A0AAV2Z8E4 | N | S | E | **H** | E | Y |
| Oomycota | Oomycota | *Nothophytophthora sp.* | A0A662X2J7 | - | S | E | **Q** | E | Y |
| Oomycota | Oomycota | *Phytophthora cactorum* | A0A329SCT9 | - | N | E | **Q** | E | Y |
| Oomycota | Oomycota | *Saprolegnia diclina* | T0REK4 | N | R | E | **H** | E | S |

## Why it matters (and the testable prediction)
In the human enzyme, mutating this exact residue (**E178Q**) gives a **"stuck-on" enzyme**: about 4× *more* active without putrescine (Bale et al. 2008). Stramenopiles carry the natural version of that change, so they are **predicted to have putrescine-independent, always-on AdoMetDC**.

That would fit diatom biology. Diatoms build their silica cell walls using long-chain polyamines, so they need a lot of polyamine output (Zhou et al. 2026, *Commun Biol*, links polyamine pools to silicification). **This is a prediction, not a tested fact.**

## Novelty status (first pass, honest)
- PubMed: **no paper reports that stramenopile AdoMetDCs lack the putrescine sensor**, and I found no eukaryote-wide survey of this residue.
- **Lin et al. 2024** (*Biomolecules*, open access) studied *Phaeodactylum* polyamine synthesis and call it "unexplored in diatoms." They overexpressed the diatom enzyme (PtSAMDC) in cells but **did not study its regulation or its pocket**. I searched their full text. They also confirm that *Phaeodactylum* has this canonical enzyme *plus* three separate fusion genes. That supports our distinction from Michael 2011, which describes the fusion machines, not this enzyme.
- **Pösö et al. 1975** found a *putrescine-insensitive* AdoMetDC in the ciliate *Tetrahymena*. That's old biochemical proof that putrescine-insensitive protist enzymes exist. It's a possible **validation point** for us (see next steps), not a pre-emption.
- **Gupta et al. 1987**: *Acanthamoeba* AdoMetDC is putrescine-*activated*, and *Acanthamoeba* **keeps Glu178** in our data. That's consistent with our prediction.
- Still to read in full: **Michael 2011** (on your ILL list).

## Honest limits
- This is a **sequence and structure prediction**. Nobody has yet measured whether stramenopile AdoMetDC actually works without putrescine.
- The sample is limited to the **top 1,000 Foldseek hits** (31 stramenopiles). A wider survey is needed.
- 31 proteins ≠ 31 independent events, because related species share ancestry. The right claim is about **how many times** Glu178 was lost, which the tree analysis will answer.
- The *Phaeodactylum* UniProt entry is labelled a **fragment**. It still covers the whole pocket, but the full gene model should be checked.

## How this becomes a computational paper, with no wet lab
1. **Wider sampling:** every stramenopile AdoMetDC in UniProt reference proteomes, not just the Foldseek top 1,000.
2. **Tree plus ancestral-state reconstruction:** where and how many times was Glu178 lost? (Once at the base of the stramenopiles? Separately in Euglenozoa, cryptophytes and bacteria?)
3. **Two independent lines of evidence:** FOLDTRACE structural calls *and* sequence alignment agree (already done for *Phaeodactylum*; I'll do it for every protein).
4. **Check against biochemistry that's already published:** do enzymes already measured as putrescine-insensitive (*Tetrahymena*, plants) or putrescine-activated (human, *Acanthamoeba*) line up with their residue at 178? If they do, that's validation without doing any new experiment.
5. **Structure figure:** compare the pocket in the human structure with the AlphaFold model of the diatom enzyme. These models already exist, so no GPU is needed.
6. **State the prediction plainly** and propose the experiment as future work.

**Later, if you want wet lab:** the Lin lab (National Taiwan Ocean University) already has *Phaeodactylum* PtSAMDC expression constructs. They'd be a natural collaborator for Marc to contact if this pans out.
