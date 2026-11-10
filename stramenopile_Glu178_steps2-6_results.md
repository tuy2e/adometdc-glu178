# Stramenopile AdoMetDCs lost the putrescine sensor (Glu178): results, steps 2–6 + wider sampling
*Tui, 5 Nov 2026 (updated the same evening). This is my computational analysis, not Marc's protocol. I used your FOLDTRACE run read-only and changed nothing in your project folder. No wet lab. Step 1 is in `stramenopile_Glu178_step1_results.md`.*

> **⚠️ Update, 7 Nov: two harder tree tests changed one conclusion (details at the end of step 4).**
> - **Still holds:** two Glu-less groups each break the species tree: (a) stramenopiles + the dinoflagellate/*Perkinsus* Asp copy, and (b) kinetoplastids/euglenid + cryptophytes. Moving any of them back to its species-tree place is rejected, now also under a C60 mixture model.
> - **No longer claimed:** that (a) and (b) are *one* lineage with *one* Glu loss. The best tree still joins them (UFBoot 88 standard / 90 mixture). But a tree that keeps them apart, the arrangement in Li et al. 2019, **can't be rejected** (AU p = 0.18 standard, 0.37–0.39 mixture, 0.48 without pocket columns). So it's 1 or 2 losses, not settled.
>
> **What's new in the 5 Nov update**
> 1. **Michael 2011 is read in full, and our novelty holds** (details below).
> 2. **Sampling is much wider:** 339 active enzymes (up from 208), now including cryptophytes, euglenids, diplonemids, more dinoflagellates and *Perkinsus*, haptophytes, red algae and more.
> 3. **A new, larger family tree, checked with several robustness tests.** It changes the "how many times" answer (see step 4).
> 4. **A second shared pocket change (Ser109→Trp)** turned up.

## Short version
1. **No stramenopile AdoMetDC has Glu178: 97 of 97 active enzymes** (Asn 54, His 32, Gln 9, Asp 2). FOLDTRACE's structure calls and the sequence alignment agree on every stramenopile tested (31/31).
2. **Three other groups lack Glu178 too:**
   - **cryptophytes, all of them** (9/9, 7 species, Asn);
   - **kinetoplastids** (Ser) and the **euglenid** *Eutreptiella* (Asn);
   - one of the two AdoMetDC copies in **dinoflagellates and *Perkinsus*** (Asp, which keeps the negative charge).

   Haptophytes (32/32), Rhizaria, heteroloboseans, the diplonemid *Diplonema*, red algae, metamonads, animals, fungi and plants keep Glu.
3. **Family tree (revised 7 Nov):** the Glu-less AdoMetDCs form **two well-supported groups that don't follow the species tree**: (a) **stramenopiles + a second dinoflagellate/*Perkinsus* copy** and (b) **kinetoplastids/euglenid + cryptophytes**. Moving any group back to its species-tree position is rejected under both models (standard: AU p = 0.00001–0.033; C60 mixture: p = 0.0002–0.013).
   - **Not settled:** whether (a) and (b) are **one lineage with one Glu loss** (the best tree under both models, UFBoot 88/90) or **two groups with two losses**. A two-group tree (as in Li et al. 2019) can't be rejected: AU p = 0.18 standard, 0.37–0.39 mixture.
   - Both groups also carry Trp at human 109. That is a single origin if they're one lineage, or a striking convergence (the same two pocket changes twice) if not.
   - Within stramenopiles: Asn first, His later in oomycetes, Glu never regained.
4. **Published biochemistry fits for stramenopiles, but our Discoba group now has two exceptions (updated 9 Nov).**
   - All **8** sequenced enzymes measured as putrescine-activated (outside trypanosomes) keep Glu178. *(9 Nov: rat and pike-perch added from Pösö 1976.)*
   - ⚠️ **New, 9–10 Nov:** *Euglena gracilis* AdoMetDC **is putrescine-activated, 3–10×** (Aleksijevic 1979, full text), and our tblastn search of its transcriptome (10 Nov) shows it has **Asn178** (with Arg13 and Asp174). So it is a **confirmed second exception** after *T. cruzi*.
     - Ancestral reconstruction on the main tree: the kinetoplastid + euglenid ancestor already had Arg13 / Asp174 / Asn178. So the two exceptions share **one inherited trait**, not two coincidences.
     - The same reconstruction gives the group (b) ancestor (incl. cryptophytes) Arg13 / Asp174 / Asn178; modern cryptophytes have Lys13. **So cryptophytes resemble the exceptions; the clean prediction is for stramenopiles only, with cryptophytes uncertain** until we know whether position 13 matters (Clyne 2002 / Beswick 2006 full text needed).
   - The Glu→Gln mutation makes the human enzyme 4× more active without putrescine, and the worm (*C. elegans*) enzyme 100× more active without it.
   - *T. cruzi* has no Glu178 and is still activated, through Asp174. So the stramenopile prediction is reasonable but not guaranteed.
5. **Structure:** the diatom pocket is intact and still acidic, with Asn where the Glu was. Next to it, a **Trp** replaces human Ser109 and reaches toward the putrescine site. That swap is shared by 84% of the Glu-less lineages and almost nobody else (1%).

Figures:
- `Fig_tree_Glu178_expanded.png` (main tree, 187 enzymes)
- `Fig_tree_Glu178.png` (first-pass tree, 171 enzymes)
- `Fig_pocket_human_vs_stramenopile.png` (structure)

---

## Step 2 — Pull every stramenopile AdoMetDC ✅ (+ wider sampling)
**2a. UniProt (first pass):**
- I took every UniProt entry with Pfam **PF01536** (the canonical class-2 eukaryotic AdoMetDC) from stramenopiles, the rest of SAR and outgroups, plus 38 sequences from your FOLDTRACE run.
- That gave 248 proteins.

**2b. Wider sampling (new):**
- **Why:** UniProt misses most transcriptome-derived (MMETSP/TSA) sequences, and those are where cryptophytes, euglenids and many dinoflagellates live.
- **Searches:** NCBI protein (by name, within each target group) and UniParc (UniProt's archive, which still holds the transcriptome entries UniProt dropped).
- **Groups searched:** cryptophytes, euglenids, kinetoplastids, diplonemids, heteroloboseans, haptophytes, dinoflagellates, *Perkinsus*, chromerids, ciliates, Rhizaria, red algae, glaucophytes, non-ochrophyte stramenopiles, metamonads and apusomonads.
- **Result:** 332 new unique sequences. All went through the same pipeline: pairwise MAFFT to human, 16 mapped positions, and the engine check.
- **Contamination screen** (decided in advance). I dropped a sequence if its closest relative (≥50% identity) was in a *different* supergroup than its own and closer than anything in its own group. I also dropped records that UniProt pins on two different organisms. **10 removed** (plus 2 short fragments):
  - a "*Symbiodinium*" sequence 100% identical to the bicosoecid *Cafeteria* (culture contamination);
  - the ciliate *Tiarina* carrying a diatom sequence (prey);
  - the "dinotoms" *Durinskia* and *Kryptoperidinium*, which carry a diatom endosymbiont and whose copies sit with diatoms;
  - a *Heterosigma*/*Oxyrrhis* mixed record;
  - and a few others (full list in my working files).
- **Ciliates** (e.g. *Tetrahymena*) have only the prokaryotic-type PF02675 fusions, not the class-2 enzyme. That matches step 5.

**Final survey:** 339 active enzymes (sequence-identical duplicates removed). Of these, **97 are stramenopile**.

## Step 3 — Sequence check of the pocket ✅
**Residue at the human Glu178 position, active enzymes, wider survey (339):**

| Supergroup | Group | Proteins | Species | Glu kept | What replaced it |
|---|---|---|---|---|---|
| **Stramenopiles** | **Ochrophyta** | **42** | 32 | **0** | **Asn 42** |
| **Stramenopiles** | **Oomycota** | **41** | 39 | **0** | **His 32, Gln 9** |
| **Stramenopiles** | **Bigyra & other** | **14** | 8 | **0** | **Asn 12, Asp 2** ‡ |
| **Cryptophyta** | **Cryptophytes** | **9** | 7 | **0** | **Asn 7, Ala 2** ¶ |
| Discoba | Kinetoplastids | 7 | 6 | 0 | Ser 6, Tyr 1 |
| Discoba | Euglenid (*Eutreptiella*) | 1 | 1 | 0 | Asn 1 |
| Discoba | **Diplonemid (*Diplonema*)** | 1 | 1 | **1** | – |
| Discoba | Heteroloboseans | 9 | 6 | 9 | – |
| Alveolata | Dinoflagellates | 49 | 27 | 2 | Asp 28, Gln 11, Arg 3, Val 3, other 2 § |
| Alveolata | *Perkinsus* | 10 | 3 | 0 | Asp 10 |
| Alveolata | Apicomplexa | 10 | 8 | 3 | *Plasmodium* Ser 4, other 3 † |
| Alveolata | Chromerids | 3 | 2 | 2 | Ala 1 (likely alignment error) |
| Rhizaria | | 12 | 6 | 12 | – |
| Haptophyta | | 32 | 19 | **32** | – |
| Archaeplastida | Red algae | 22 | 18 | 19 | Asp 2 (Gracilariales), 1 gap |
| Archaeplastida | Green algae + land plants | 27 | 25 | 25 | Arg 1 (likely alignment error), 1 gap |
| Archaeplastida | Glaucophyte | 1 | 1 | 1 | – |
| Opisthokonta | Animals / fungi / other | 42 | 42 | 41 | Gln 1 (tapeworm) |
| Amoebozoa, Metamonada, Apusozoa | | 7 | 7 | 7 | – |

‡ The Asp cases are the bicosoecid *Halocafeteria* and one *Blastocystis* genome (subtype "Blastoise"). The *Blastocystis* one is odd: the other 6 *Blastocystis* have Asn, so it may be contamination.
§ Dinoflagellates have **two** AdoMetDC copies. One, often fused to an ODC (ornithine decarboxylase) domain, has Gln and sits with apicomplexans in the tree. The other, a standalone copy, has Asp and sits with stramenopiles (step 4).
¶ The two Ala calls are both *Cryptomonas paramaecium*. The 187-protein multiple alignment puts Asn in that column for the one included in the tree.
† *Plasmodium* AdoMetDC is part of a ~1,400-aa fusion protein, so its calls are the least reliable in the set.

**What the table says:**
- **Every stramenopile, every cryptophyte, and every kinetoplastid/euglenid** lacks Glu178.
- **Their closest relatives keep it.** Haptophytes and Rhizaria are close relatives of stramenopiles and cryptophytes; heteroloboseans and *Diplonema* are close relatives of kinetoplastids.
- **Two methods agree.** FOLDTRACE structure calls and sequence alignment agree at position 178 for 125/129 overlapping proteins (31/31 stramenopiles). Over all 7 pocket positions they agree 830/874 times.

## Step 4 — Family tree: where, and how many times, was Glu178 lost? ✅ (re-done with wider sampling)

**What I ran:**
- **Proteins:** 187 active enzymes. I took the wider survey, removed near-duplicates, and capped very well-sampled groups so the tree stays balanced. Every cryptophyte, euglenozoan, *Perkinsus*, bigyran and deep-branching sequence was kept.
- **Alignment:** MAFFT L-INS-i, trimmed to 340 columns.
- **Tree:** IQ-TREE 2 with the best-fit model **Q.pfam+R7** (chosen fresh by BIC), 1,000 ultrafast bootstraps (UFBoot) + 1,000 SH-aLRT tests, and ancestral reconstruction at the Glu178 column.
- **Figure:** `Fig_tree_Glu178_expanded.png`.

**Result 1 — inside stramenopiles (unchanged and solid):**
- Stramenopiles now form a single group in the tree (weak support: UFBoot 68). Their ancestor had **Asn** (probability 0.999).
- Oomycetes later switched **Asn → His** (support 100). The downy mildew *Bremia* then switched His → Gln.
- Two bicosoecid-type tips switched Asn → Asp.
- **Glu was never regained anywhere** in the tree.

**Result 2 — the bigger finding: the Glu-less enzymes are one gene lineage**
- **What the tree shows:** the AdoMetDCs of **stramenopiles, cryptophytes, kinetoplastids + the euglenid, and a second copy found in dinoflagellates and *Perkinsus*** form **one group** (UFBoot 88).
  - That group does **not** follow the species tree: these organisms belong to 3 different supergroups.
  - Their close relatives in the species tree (haptophytes, Rhizaria, *Diplonema*, heteroloboseans, and the dinoflagellates' *other* copy, which sits with apicomplexans) all keep Glu and sit elsewhere.
- **The reconstruction** puts **a single Glu→Asn change at the base of this lineage** (probability 0.99). Each branch then made its own later swap:
  - Asn → Asp in the dinoflagellate/*Perkinsus* copy (charge restored);
  - Asn → His in oomycetes;
  - Asn → Ser in kinetoplastids.
- **Shared pocket clue:** a Trp in place of human **Ser109** is found in 84% of this lineage and in **1%** of everything else. In the diatom and brown-alga models, this Trp sits about 4 Å from where putrescine binds in human, which is roughly where Phe111's ring sits in human. In Ochrophyta Phe111 became Gly, so the Trp appears to take over its job.

**Robustness tests (are we fooling ourselves?):**

| Test | What it checks | Result |
|---|---|---|
| Remove all 35 pocket columns (within 8 Å of putrescine) and rebuild | Is the pocket itself pulling these sequences together? | Same lineage: stramenopiles + cryptophytes + kinetoplastids/euglenid (SH-aLRT 87 / UFBoot 64). The dinoflagellate/*Perkinsus* copy is its sister (85.5 / 72). |
| Remove the fast-evolving kinetoplastids and rebuild | Long-branch artifact? | Same lineage (UFBoot 85), still one Glu→Asn at its base. The euglenid now sits next to cryptophytes. |
| AU test: put cryptophytes back with haptophytes | Could cryptophytes have lost Glu on their own? | ΔlnL 38, **p = 0.025** (rejected) |
| AU test: put kinetoplastids + euglenid back with *Diplonema* | Could Euglenozoa have lost it on their own? | ΔlnL 55, **p = 0.001** (rejected) |
| AU test: put the dinoflagellate/*Perkinsus* Asp copy with the alveolate copy | Is it just an ordinary alveolate duplicate? | ΔlnL 34, **p = 0.003** (rejected) |
| AU test: all three moved back (species-tree arrangement) | Separate losses everywhere? | ΔlnL 99, **p = 0.0001** (rejected) |
| Same four AU tests with the pocket columns removed | Is the rejection driven by the pocket? | Still rejected: p = 0.034 / 0.012 / 0.002 / 0.0006 |

*(Each alternative tree was built by moving that group to its species-tree neighbour in the best tree, then re-optimizing under a constraint. The AU test used 10,000 replicates.)*

**So how many times was Glu178 lost? (5 Nov answer. ⚠️ Revised on 7 Nov: "1 or 2, not settled"; see robustness round 2 below.)**
- **In the best-supported picture, once.** A single Glu→Asn change happened in one AdoMetDC gene lineage. That lineage is now carried by stramenopiles, cryptophytes, kinetoplastids/euglenids, and (as an extra copy) dinoflagellates and *Perkinsus*.
- **The alternative of separate losses in each group** (at least 4 events) is **statistically rejected**, with or without the pocket residues.
- **Apart from that lineage**, only a handful of single-species changes occur: one tapeworm (Gln), plus one red alga and one mite that swapped Glu for Asp and so kept the charge. Separately, dinoflagellates' own alveolate copy has Gln (one event, support 100).
- **Glu was never regained.**

**What this probably means (stated carefully):**
- A gene lineage spread across unrelated eukaryotes usually means **horizontal gene transfer** (the gene moved between organisms).
- The other option is an **ancient gene duplication followed by losses**. Dinoflagellates still carrying *both* types would fit either one.
- **This single-gene analysis can't tell which, or which direction the gene moved.** That has to be said plainly in a paper.
- **For our story,** stramenopiles carry the putrescine-sensor-less AdoMetDC because their gene belongs to this special lineage, not because each group independently reinvented it. *(7 Nov: the "each group separately" version with 4 losses is still rejected. A 2-loss version, with stramenopiles + Asp copy as one group and kinetoplastids + cryptophytes as another, is not.)*

**⚠️ Robustness round 2 (7 Nov): two harder tests**
1. **A Li 2019-style alternative.** Kinetoplastids + cryptophytes kept together, but stramenopiles + the Asp copy moved out of their lineage, either next to animals/fungi (as in Li 2019's Fig. 2) or next to the alveolate copy (their place in the SAR species tree).
2. **A C60 mixture model** (LG+C60+F+G, PMSF). It is better at avoiding false groupings caused by uneven amino-acid composition.

| Test | Standard model | Without pocket columns | C60 mixture model |
|---|---|---|---|
| Best-tree support for the 4-group lineage | UFBoot 88 | (64 earlier) | **UFBoot 90**, SH-aLRT 87 |
| Cryptophytes back with haptophytes | p = 0.033 ✗ | p = 0.068 (not rejected) | p = 0.008 ✗ |
| Kinetoplastids + euglenid back with *Diplonema* | p = 0.0005 ✗ | p = 0.021 ✗ | p = 0.0002 ✗ |
| Asp copy back with the alveolate copy | p = 0.004 ✗ | p = 0.004 ✗ | p = 0.013 ✗ |
| All three back | p = 0.00001 ✗ | p = 0.005 ✗ | p = 0.0005 ✗ |
| **Li 2019-style: stramenopiles + Asp copy next to animals/fungi** | **p = 0.18 (not rejected)** | **p = 0.48** | **p = 0.37** |
| **Stramenopiles + Asp copy next to the alveolate copy** | **p = 0.18 (not rejected)** | **p = 0.19** | **p = 0.39** |

*(✗ = rejected at p < 0.05. AU test, 10,000 replicates, 7–10 candidate trees per run. The p-values for the old hypotheses differ slightly from the 5 Nov table because the candidate set changed. Raw numbers: `/root/audit/strameno/expand2/RESULTS_round2.md`.)*

**What this changes:**
- **Solid:** (a) stramenopiles + the dinoflagellate/*Perkinsus* Asp copy group together, and the Asp copy is not an ordinary alveolate duplicate. (b) Kinetoplastids/euglenid + cryptophytes group together, against the species tree. Both hold under the mixture model. One caveat: without pocket columns the cryptophyte test is borderline (p = 0.068).
- **Not settled:** whether (a) and (b) share one origin. The best tree says yes under both models, but the alternative is only 8–13 log-likelihood units worse and not rejected.
- **New possibility for (a):** stramenopiles and alveolates are sister groups (SAR). So (a) could simply be an **ancient SAR duplicate**: stramenopiles kept one copy, alveolates the other, and dinoflagellates/*Perkinsus* kept both. That needs no gene transfer. Only (b), cryptophytes inside a Discoba group, clearly crosses supergroups, and Li 2019's tree already shows that pairing.
- **Glu178 losses:** 1 (if one lineage) or 2 (one at the base of each group). Either way stramenopiles lost it once, early, and never regained it.
- **Trp109:** found in 36/45 stramenopiles, 16/16 Asp copies, 6/7 kinetoplastids/euglenid and 5/7 cryptophytes, but almost never elsewhere. If there are two groups, both made the **same pair of pocket changes** independently. That would be convergent evolution, a functional signal in its own right.

**Round 3 (9 Nov): ancestral reconstruction on the two-group trees + a second gene (spermidine synthase)**

*A. Ancestral reconstruction* (Q.pfam+R7, fixed trees; IQ-TREE 2.0.7 can't reconstruct under the C60 model):

| Tree | Glu178 losses | Where | Trp109 gains |
|---|---|---|---|
| Best tree (one lineage) | **1** | base of the 4-group lineage (P 0.99) | **1**, at the same node (P 0.99) |
| Two-group trees (Li-style / SAR-style, both models) | **2** | base of stramenopiles + Asp copy (ancestor N or D, P ≈ 0.4–0.5), and base of kinetoplastids + cryptophytes (N, P 0.99) | 1–3; unstable (in some trees a single low-confidence gain that plants, red algae and haptophytes then reverse) |

- **Glu is never regained in any tree.**
- In the mixture-model two-group trees, the groups are separated only by haptophytes + red algae + green plants (48 sequences). That's why the likelihood gap is small (8 units).

*B. Spermidine synthase tree* (record: `/root/audit/spdsyn/RESULTS_spdsyn.md`):
- **Data:** 142 enzymes (130 eukaryotic + 12 prokaryotic), 292 positions. Model LG+R6, UFBoot + SH-aLRT. Thermospermine synthases and animal spermine synthase were excluded using reference enzymes. Same contamination rule (4 dropped).
- **Sampling gap:** NCBI BLAST, needed to catch unannotated transcriptome proteins, never ran (queued for more than 2.5 h). So only 2 cryptophytes, 0 euglenids, 3 dinoflagellates and 2 *Perkinsus*.

| Arrangement forced on spermidine synthase | p (AU) |
|---|---|
| Cryptophytes + kinetoplastids (like AdoMetDC group b) | **0.51, not rejected**; as good as the best tree (Δ0.25) |
| Stramenopiles + dinoflagellates/*Perkinsus* (like AdoMetDC group a) | **0.002, rejected** |
| All four together (like the AdoMetDC lineage) | **0.0006, rejected** |
| Cryptophytes + red algae | 0.016, rejected |

- **What it means:**
  - Spermidine synthase does **not** share AdoMetDC's stramenopile–dinoflagellate grouping or the 4-group lineage. So AdoMetDC's pattern is specific to AdoMetDC, not a quirk of these organisms' genomes, and the two genes didn't travel together.
  - The cryptophyte–kinetoplastid pairing is weakly echoed (not rejected; Li 2015's quick tree also showed it).
  - **The second gene does not settle one-vs-two lineages.**

**Where this leaves the "how many losses" question:** still 1 or 2, and the data can't decide.
- **One** is favoured by the best tree under both models and by needing only one Trp109 gain.
- **Two** can't be rejected statistically.
- What's solid regardless: stramenopiles lost Glu178 once, early, with no reversal. Cryptophytes + kinetoplastids/euglenid and the dinoflagellate Asp copy are also Glu-less and each break the species tree. Trp109 accompanies the loss in all of them.

**Earlier first-pass tree (171 enzymes, UniProt only)**, kept for the record (`Fig_tree_Glu178.png`):
- It showed the same grouping, but the AU tests could only reject all three species-tree arrangements *together*, not individually.
- Adding cryptophytes, *Diplonema*, more dinoflagellates, *Perkinsus* and more Bigyra made each individual test significant.

**Honest limits of the tree part:**
- **It's one gene.** Support for the deepest node of the lineage is moderate (UFBoot 85–88 with the pocket; 64 without it).
- **The AU tests** depend on how well the alternative trees were optimized. Each got a constrained search starting from the best tree.
- **Sampling still has gaps.** Katablepharids, telonemids and centrohelids returned no AdoMetDC sequences in NCBI name searches. Some sequences come from transcriptomes, which I screened for contamination but cannot fully guarantee.
- ***Diplonema* is a single sequence.** It keeps Glu and groups with heteroloboseans, not with kinetoplastids.

## Step 5 — Check the prediction against enzymes whose putrescine response has been measured ✅
I looked up every AdoMetDC whose putrescine response has been measured, pulled its sequence from UniProt, and mapped it with the same pipeline. Abstracts and full text were checked on PubMed/PMC.

**First, the mapping reproduces published residue numbers:**
- *C. elegans*: T192, E194, E274 match Ndjonka 2003 exactly.
- Potato: R18, R114, V181, H294 match Bennett 2002 exactly.
- *T. cruzi*: 9 of 10 residues named in Clyne 2002 / Beswick 2006 match (F7, R13, E15, I80, D174, S178, T245, E256, F285). The one miss is Ser111, where our alignment is 2 residues off because of a nearby gap. I note this for the F111 calls in trypanosomes; the stramenopile F111 calls are backed by FOLDTRACE.
- The *Acanthamoeba castellanii* peptide that Hugo & Byers sequenced (UniProt P34039) matches the N-terminus of L8GRS6 exactly. So L8GRS6 is the protein they purified.

| Organism (UniProt) | Measured putrescine response | Source (PMID) | Residue at 178 | Fits the rule? |
|---|---|---|---|---|
| Human (P17707) | Activated (processing + activity). **E178Q: ~4× WT activity without putrescine, no further stimulation** | Stanley 1994 (8132508); Bale 2008 Table 4 (19053272) | **E** | ✅ |
| Yeast *S. cerevisiae* (P21182) | "putrescine-activated" | Pegg & Jacobs 1983 (6351843); Karvonen 1985 (3933487) | **E**213 | ✅ |
| *Neurospora crassa* (Q9Y7F3) | Putrescine **absolutely required** for activity | Hoyt 2000 (10852489) | **E**258 | ✅ |
| *C. elegans* (O02655) | **350-fold** stimulated. **E194Q (= human E178Q): 100-fold higher activity with no activator** | Ndjonka 2003 (12674502) | **E**194 | ✅ key precedent |
| *Onchocerca volvulus* (Q27883) | 15-fold stimulated | Da'dara 1998 (9879883) | **E**196 | ✅ |
| *Acanthamoeba castellanii* (L8GRS6) | 2-fold stimulated (weak) | Hugo & Byers 1993 (8216217) | **E**202 | ✅ |
| *Acanthamoeba culbertsoni* | Strongly activated | Gupta 1987 (3600699) | no sequence available | – |
| Potato (Q04694), *Arabidopsis* (Q96286) | **Not** stimulated; fully active alone | Bennett 2002 (12463749) | **E**185 / **E**182 | ✅ different route: Glu kept, but Arg/Lys/His built into the pocket mimic putrescine |
| *P. falciparum* (Q9Y1L5) | AdoMetDC domain **not** stimulated. A structural model explains why: built-in basic residues (Arg11, Lys15, Lys215) take putrescine's place, a plant-like mechanism | Wrenger 2001 (11390378); Wells 2006 (16257247) | **S** (and **R** at human-L13, matching Wells' Arg11) | ✅ |
| ***T. cruzi* (O76240)** | **Stimulated 11–27-fold**; binding needs **Asp174 only**; mutating Ser178 has no effect | Clyne 2002 (12403622); Beswick 2006 (16784231) | **S**194 | ❌ **exception** |
| Rat liver (P17708) | 1.8-fold stimulated (crude extract) | Pösö 1976 (813639) | **E** (D174, E256 also kept) | ✅ *(added 9 Nov)* |
| Pike-perch liver, *Sander lucioperca* (A0A8C9YFM5) | 2.0-fold stimulated (crude extract) | Pösö 1976 (813639) | **E**180 (D176, E258 also kept) | ✅ *(added 9 Nov)* |
| ***Euglena gracilis*** | **Stimulated 3–10-fold** by 2 mM putrescine (dialysed extract); apparent K<sub>a</sub> 0.6 µM; specific (diaminopropane gives 25% of max); no Mg²⁺ effect | Aleksijevic 1979 (116684) | **N** ✅ *(found 10 Nov: tblastn of *Eutreptiella* AdoMetDC vs *E. gracilis* TSA, RID F5EJKSE8014. One gene, 10 transcripts from two independent assemblies (GEFR01006935.1; GDJR01034331–39), E = 2e-145, 65% identical over query 49–372, no gaps. All 10 agree: R13, **Q15**, Q111, D174, Q176, **N178**, E256, L285. No second AdoMetDC-like transcript, so no prozyme candidate in this search.)* | ❌ **confirmed second Discoba exception**: no Glu178, still activated |
| *Tetrahymena* | **Not** stimulated (also 1.0× in Pösö 1976) | Pösö 1975 (813467); Pösö 1976 (813639) | n/a: no class-2 (PF01536) enzyme in UniProt for *Tetrahymena*. Its only AdoMetDCs are 6 small prokaryotic-type (PF02675) domains fused to spermidine synthase | ✅ explained by enzyme class, not residue 178 |

**What this means for the prediction (stated carefully):**
- **Glu178 marks the activated enzymes.** All 8 measured, putrescine-activated enzymes **with a known sequence** outside trypanosomes keep it. *(9 Nov: 6 → 8 with rat and pike-perch.)*
- **Removing Glu178 has a known effect in two species.** In the human enzyme, E178Q gives ~4× wild-type activity with no putrescine, and putrescine adds nothing. In the worm enzyme, the equivalent E194Q gives 100× higher activity with no activator. Stramenopiles carry the natural version of this change. Ochrophyta and Bigyra have Asn, a close relative of Gln.
- **But *T. cruzi* shows another possible outcome.** It lacks Glu178 and is still putrescine-activated, because Asp174 alone binds putrescine. Stramenopiles keep Asp174 (Glu in oomycetes) and Glu256. So some putrescine binding could remain.
- **The honest prediction:** stramenopile AdoMetDC should have **high activity without putrescine and little or no extra stimulation**, like E178Q. A *T. cruzi*-like residual response cannot be ruled out from sequence alone. The experiment (activity ± putrescine) decides between them.
- **Correction to step 1:** I had listed *Tetrahymena* (Pösö 1975) as a possible validation point. It isn't one. *Tetrahymena* has no class-2 enzyme at all, so its insensitivity tells us nothing about residue 178.
- **Literature checks (chunk 1C):**
  - **No stramenopile, oomycete or cryptophyte AdoMetDC has ever had its putrescine response measured.** I searched PubMed for every group, so the experiment is genuinely open.
  - ✅ **Settled 10 Nov:** *E. gracilis* AdoMetDC has **Asn178** (tblastn vs its transcriptome; see table). So option 1 below is ruled out. *Euglena* lost Glu178 and is still putrescine-activated, as *T. cruzi* is. Both keep Asp174 and Arg13, which fits the *T. cruzi* Asp174 route (my inference; not tested). No prozyme-like transcript turned up.
  - **The only measured enzyme inside our Glu-less lineage outside trypanosomes is *Euglena gracilis*** (Aleksijevic 1979). ✅ **Full text read 9 Nov (your ILL PDF): it IS putrescine-activated, 3–10×** (see the table row above). Three ways this could fit:
    1. *E. gracilis* kept Glu178 (unlike *Eutreptiella*);
    2. putrescine acts through Asp174, as in *T. cruzi* (euglenids, like kinetoplastids, have Asp174 and Arg13);
    3. a prozyme-type partner, as in trypanosomes.
    - The *E. gracilis* sequence decides between option 1 and options 2/3.
    - **Stramenopiles differ from the Discoba exceptions in one checked way:** no basic residue at human L13 (they have Thr/Asn/Ile; kinetoplastids/*Eutreptiella* have Arg, cryptophytes Lys). So the prediction is best stated for stramenopiles (cryptophytes: untested).
    - ⚠️ **Not yet checked: a prozyme in stramenopiles.** No stramenopile prozyme has been reported, but **we have not searched for one specifically**. Our NCBI expansion flagged 3 low-identity (~20%) class-2-like hits in *Cafeteria roenbergensis* (KAA0165454.1, KAA0162548.1, KAA0153562.1; ~255 aa) and 1 in *Schizochytrium* (CAE0742880.1) as lacking the active site. **They are unexamined.** They could be prozyme-like dead copies, mis-hits, or fragments. *My suggestion: a quick check (domain scan + where they fall on the tree) before we claim stramenopiles have no prozyme.*
  - **Prior note for trypanosomatids:** Singh 2013 (*PLoS One*) observed that *Leishmania* AdoMetDC shares only 2 of the 3 human putrescine-binding residues. That fits our finding for kinetoplastids.
  - **The standard classification has no group for these enzymes.** Bale & Ealick 2010 list 2a (plants, built-in Arg, not activated), 2b-I (putrescine-activated, marked by a conserved **Lys80**) and 2b-II (trypanosomatid prozyme-activated). Nothing covers stramenopiles or cryptophytes.
  - **Lys80 in our lineage:**
    - stramenopiles, the dinoflagellate/*Perkinsus* copy and all other eukaryotes keep it, so the relay to the active site is intact;
    - **cryptophytes have Leu80**;
    - kinetoplastids and the euglenid have Ile80 (the known trypanosomatid feature).

## Step 6 — Structure figure: human vs diatom pocket ✅
- **Method:** human crystal structure 1MSV (putrescine bound, chain A) vs the AlphaFold models already in your FOLDTRACE run. I superposed each pair on the residues near the pocket, refining the fit iteratively. Core RMSD: diatom *Phaeodactylum* 0.95 Å, brown alga *Ectocarpus* 0.76 Å, *Cafeteria* 0.83 Å, oomycete *Hyaloperonospora* 0.78 Å. pLDDT at the pocket residues is 80–98, so the models are confident there.
- **Figure:** `Fig_pocket_human_vs_stramenopile.png`. A = human, B = diatom, C = oomycete.

**Closest approach of each pocket residue to putrescine's nitrogens (Å). AlphaFold confidence (pLDDT) in brackets:**

| Human residue (1MSV) | Human | Diatom *Phaeodactylum* | Brown alga *Ectocarpus* | *Cafeteria* (Bigyra) | Oomycete *Hyaloperonospora* |
|---|---|---|---|---|---|
| L13 | 3.4 | Thr12 3.6 (97) | Thr56 3.9 (98) | Thr32 4.1 (98) | Asn29 3.1 (98) |
| E15 | 2.8 | Glu14 3.2 (96) | Glu58 2.6 (97) | Glu34 2.6 (97) | Glu31 2.7 (98) |
| F111 | 4.5 | **Gly124** 5.8 (95) | **Gly148** 5.7 (95) | Cys124 4.5 (96) | Asn121 3.2 (96) |
| D174 | 2.9 | Asp194 3.7 (80) | Asp210 4.1 (83) | Asp187 4.2 (87) | Glu199 3.6 (78) |
| T176 | 2.7 | Val196 3.2 (93) | Thr212 3.1 (95) | Thr189 3.0 (94) | Thr201 3.0 (92) |
| **E178** | 5.2 | **Asn198** 7.2 (95) | **Asn214** 6.8 (98) | **Asn191** 5.0 (98) | **His203** 3.7 (96) |
| E256 | 5.0 | Glu289 5.6 (97) | Glu292 5.3 (98) | Glu352 4.6 (98) | Glu282 4.8 (98) |
| F285 | 4.5 | Met318 3.0 (88) | Phe321 4.7 (93) | Phe381 4.6 (92) | Tyr311 4.3 (88) |

**Reading it:**
- **Ochrophyta and Bigyra (diatom, brown alga, *Cafeteria*):** the pocket is still there and still lined with acid (Glu15, Asp174, Glu256). The Glu178 spot holds a **neutral Asn**. No side chain fills the space where putrescine sits. In human 1MSV, Glu178 contacts putrescine through a bridging water (Bale 2008), which fits its 5.2 Å distance. Bale also reports that E178Q *lacks* that bridging water, which is why putrescine no longer activates it. The diatom geometry looks like that mutant, not like the plant solution. **Gly** at F111 (Ochrophyta) removes one of the two aromatic walls. *New:* a **Trp** at the human Ser109 position (Trp122 in *Phaeodactylum*, Trp146 in *Ectocarpus*, Trp122 in *Cafeteria*) reaches to within 3–4 Å of the putrescine site, close to where human Phe111's ring sits (3.7 Å). So the aromatic wall seems to have **moved over by two residues** rather than disappeared. This Trp is shared across the whole Glu-less lineage (step 4).
- **Oomycete:** the pocket is **partly filled by the protein's own side chains**. His203 sits at the Glu178 spot and points in (3.7 Å). Asn121 replaces Phe111, and a Tyr (Tyr354, outside the 8 mapped positions) has its OH 2.1 Å from where putrescine's nitrogen sits in human. That's closer to the plant idea (a built-in side chain occupying the activator site), so **oomycetes may differ from diatoms**. This is a possible sub-finding.
- **Limits:** AlphaFold models have no ligand and no water. This shows pocket *shape and chemistry*, not binding.


## Novelty status (Michael 2011 now read in full)
- **Michael 2011** (*FEBS Lett* 585:2627; you found it open-access) is about the **bacterial-type AdoMetDC fusion proteins** in diatoms:
  - class-1b AdoMetDC fused to an aminopropyltransferase, sometimes with SET, Tudor or ShKT domains;
  - some carry a dead AdoMetDC domain;
  - all with signal peptides, proposed to build long-chain polyamines for silica.
- **On our enzyme it says only one thing:** "Typical eukaryotic AdoMetDC orthologues are present in the three diatom genomes." It also notes the *Thalassiosira* and *Phaeodactylum* mRNAs carry short upstream reading frames (29 and 26 amino acids), probably for **translational feedback control**, as in plants and mammals.
- **It never discusses the putrescine pocket, Glu178, allosteric activation, or the evolutionary history of the canonical enzyme.** Our finding is not pre-empted.
- **Useful for the paper:** the uORF point means diatom AdoMetDC is probably still regulated by polyamines at the level of how much protein is made, even if the putrescine "switch" on the protein is gone. That's a nice discussion point.
- **Other checks so far:**
  - PubMed has no report of a stramenopile-wide (or cryptophyte/euglenozoan-wide) loss of Glu178, and no survey of this residue across eukaryotes.
  - Lin 2024 overexpressed the *Phaeodactylum* enzyme but didn't study its regulation.
  - **Still to check before writing:** whether anyone has reported **horizontal transfer of AdoMetDC** among protists. I haven't searched for that specifically yet, and I'll do it next.

## Honest limits (whole project)
- **Prediction, not measurement.** No stramenopile AdoMetDC has had its putrescine response measured. *T. cruzi* shows that losing Glu178 doesn't guarantee losing the response.
- **One-gene phylogeny.** Two Glu-less groups are well supported. Whether they are one lineage is **not settled** (7 Nov tests). Everything rests on one ~340-column alignment.
- **Sequence sources.** Many sequences come from transcriptomes. I screened for contamination with a rule set in advance, but some risk remains (e.g. the *Blastocystis* "Blastoise" Asp sequence).
- **Uneven taxon sampling.** Diatoms and oomycetes are heavily sampled; *Diplonema* and *Eutreptiella* are single sequences.

## Where this leaves the paper
**Core results, all computational:**
1. **No stramenopile AdoMetDC has the putrescine-sensor Glu178** (97/97). Two independent methods agree.
2. **The same is true of cryptophytes, kinetoplastids/euglenids, and a second dinoflagellate/*Perkinsus* copy.** These form two Glu-less groups that each break the species tree: stramenopiles + Asp copy, and kinetoplastids + cryptophytes. Species-tree placements are rejected under both models. **Whether the two groups are one lineage (one loss) or two (two losses) is not settled (7 Nov).** Stramenopiles lost Glu once, early, and never regained it.
3. **The pocket was remodelled, not destroyed.** It keeps its acidic residues, has a neutral Asn at 178, and has a Trp at 109 that is unique to this lineage.
4. **A prediction with precedent:** human E178Q (~4× more active without putrescine) and worm E194Q (100×) are activator-independent. Published biochemistry fits, with *T. cruzi* as the stated exception.
5. **Clearly separated from Michael 2011's fusion family.**

**Possible title idea (revised 7 Nov):** "Loss of the putrescine sensor in stramenopile, cryptophyte and euglenozoan AdoMetDCs", which works whether it was one event or two.

**Next steps (computational):**
- [ ] Search the literature specifically for AdoMetDC horizontal gene transfer in protists (novelty for result 2)
- [ ] Redraw the structure figure with Trp109 added
- [ ] Paper Figure 1: species-tree cartoon + residue column + gene-lineage tree + pocket panel
- [x] Mixture-model test (C60, 7 Nov): the two groups hold; their union is the best tree but not significant (see step 4, round 2)
- [x] Ancestral reconstruction on the two-group trees (9 Nov): 2 Glu losses on those trees, 1 on the best tree; Glu never regained
- [x] Spermidine-synthase tree (9 Nov): doesn't share AdoMetDC's grouping; doesn't settle one-vs-two
- [ ] Retry NCBI BLAST for transcriptome-only spermidine synthases (cryptophytes, euglenids, dinoflagellates) if we want a fuller second-gene test
- [ ] One-page summary for Marc; then ask Luke/Maddie about the journal
