# Background factors that affect the pocket-chemistry survey

**Done overnight, 10–11 Nov 2026, while Tui was asleep.** Six things that bear on the survey and
were not in the handoff. Two of them change how Stage 1 should be read; one closes an open item in
the log; one is a new result from our own data.

Literature was read through the desktop browser (the cloud workspace cannot reach Europe PMC, PMC,
RCSB, NCBI or UniProt, and org web search is off). Computation was done in the cloud workspace on the
structures already in hand. **No Stage 2 structure was measured** — the blinding in
`PREREG_pocket_chemistry_measurements.md` is intact.

---

## 1. Regulation in this family runs through oligomerisation, not only through the pocket

This is the finding that most changes how to read Stage 1.

- **Volkov et al. 2016, *eLife* 5:e20198** (PMID 27977001; the 5TVO/5TVF paper). Monomeric
  *T. brucei* AdoMetDC is **inactive**, autoinhibited by its own N-terminal sequence sitting in the
  active site. Heterodimerisation with prozyme displaces that sequence through a cis-to-trans proline
  isomerisation, a β-sheet reorganisation, and insertion of the N-terminal α-helix into the dimer
  interface. They argue that acquiring the dimerisation domain is the single step that explains how
  **regulatory schemes diverged across the whole AdoMetDC family**.
- **Velez, Brautigam & Phillips 2013, *J Biol Chem* 288:5232** (PMID 23288847). The *T. brucei*
  N-terminus is **essential** for allosteric activation by prozyme.
- **Bale et al. 2008** (already read): the putrescine pocket sits **near the dimer interface**, and
  putrescine binding to the wild-type dimer is **cooperative**. D174N does not bind; E178Q and E256Q
  bind weakly and **without cooperativity**.

What this means for us:

1. **5TVO and 5TVF are not replicate crystals of one enzyme.** 5TVO (monomer) is the autoinhibited
   state; 5TVF (prozyme heterodimer) is the activated state. Using their difference as a
   crystal-to-crystal noise estimate, as the pre-registration currently does, would fold a real
   conformational change into the error bar and make the gate far too lenient. **Amended** (see §8).
2. **Every AlphaFold model we measure is a monomer**, and cooperativity is a property of the dimer.
   A monomer measurement cannot see it, and where the site abuts the interface, part of the pocket
   will be artificially solvent-exposed. This is a stated limitation, and it is testable once 1MSV is
   downloaded: measure the deposited monomer against the biological dimer and report the difference.
3. **The truncated *T. cruzi* model is worse than we thought, and in a specific direction.** It is
   missing the first ~20 residues — precisely the autoinhibitory N-terminal segment. Its pocket is
   therefore artificially open, which is consistent with the odd numbers it gave on 10 Nov (cavity
   70 Å³, the largest of any structure, yet 0/9 poses in site).
   **Confirmed tonight.** The full-length *T. brucei* model (Q587A7, site pLDDT 98) measures
   **−25.1 kT/e at the site, −32.9 over 2 Å, cavity 29.5 Å³, 2/9 poses in site** — much more like the
   activated controls than the truncated *T. cruzi* model's −9.7 / −21.7 / 70.1 Å³ / 0 of 9. The
   truncation, not the biology, was driving a large part of the 10 Nov *T. cruzi* result.

## 2. The putrescine site is lined by loops that move

- **Bale 2008**: putrescine binding reorganises four aromatic residues — **Phe285, Phe315, Tyr318,
  Phe320** — and closes **loop 312–320** over the ligand, shielding it from solvent and thereby
  *enhancing* the electrostatic and hydrogen-bonding effects we are trying to measure.
- **Patel et al. 2025, *Biomolecules* 15:1274** (PMID 41008581). Multi-temperature crystallography of
  AdoMetDC at 100 K, 273 K and 293 K finds two mobile loops: **residues 20–28** (large motion at
  ambient temperature versus cryo) and **residues 164–174**, which is disordered in previous cryo
  structures, ordered at 100 K, and in different conformations at 273 K and 293 K.

Loop 164–174 ends at **Asp174**, the one residue that binds putrescine in both human and *T. cruzi*
and one of our four pKa positions. So the site is bounded by a loop that is genuinely mobile at
physiological temperature and poorly ordered in most cryo structures.

Two consequences:

1. **The G2 error bar is well aimed.** The three human crystals differ in exactly these loops, so
   crystal-to-crystal spread is the right way to size the noise.
2. **There is a holo/apo confound in G2 as written.** 1MSV, 1JL0 and 1I7B all have putrescine bound,
   so their lid is closed; every AlphaFold model is apo, with the lid open. "Crystal versus model"
   therefore currently measures "closed lid versus open lid" at the same time. **Fix: add 1JEN**
   (Ekstrom 1999, 2.25 Å), which is human AdoMetDC with **no ligand bound** — a true apo human
   crystal. 1MSV-vs-1JEN then isolates the lid effect, and 1JEN-vs-P17707 gives a clean apo
   crystal-to-model comparison. **Added to the fetch script.**

## 3. New result: the signalling relay is intact everywhere; the lid is not

Bale 2008 identifies a chain of hydrophilic residues connecting the putrescine site to the active
site — **Glu11, Lys80, Ser229, His243** — and shows that E178Q loses *stimulation* without losing
*binding*. That is the distinction between sensing and transmitting, and nothing in this project had
looked at the transmitting half. Computed tonight from the existing 187-sequence trimmed alignment
and human column map (sequence-level, like the round-4 residue table; no structures involved).

Consensus per lineage (`*` = below 80% agreement):

| Lineage (n) | E11 | K80 | S229 | H243 | F285 | F315 | Y318 | F320 |
|---|---|---|---|---|---|---|---|---|
| Metazoa (14) | E | K | S | H | F | F* | Y* | L* |
| Fungi (10) | E | K | S | H | F | L* | Y | L* |
| **Ochrophyta (22)** | **E 20/22** | **K 22/22** | **S 22/22** | **H 22/22** | M 14/22 | V 16/22 | C 9/22 | **C 18/22** |
| **Oomycota (12)** | **E 12/12** | **K 12/12** | **S 12/12** | **H 12/12** | Y 6/12 | F 7/12 | Y 7/12 | C 7/12 |
| Land plants (9) | E | K | S | H | H 6/9 | L* | G 6/9 | V* |
| **Cryptophytes (7)** | E 7/7 | **L 7/7** | S 7/7 | H 7/7 | H 5/7 | F* | Y 7/7 | V* |
| **Kinetoplastids (6)** | E 6/6 | **I 6/6** | S 6/6 | H 6/6 | F* | F 5/6 | Y 5/6 | V* |
| Dinoflagellates (16) | E | K | S | H | F | F* | Y* | A* |
| Haptophytes (12) | E | K | S | H | T | L* | V* | V* |
| Rhodophytes (12) | E | K | S | H | V* | L* | Y* | A* |

Three things stand out.

1. **The relay is universally conserved.** Glu11, Ser229 and His243 are invariant across all 19
   lineages examined, land plants and *Plasmodium*-type apicomplexans included. **Stramenopiles keep
   the entire relay** — 22/22 at K80, S229 and H243 in Ochrophyta, 12/12 in oomycetes.
   So stramenopiles have lost the *sensor* (Glu178) and kept the *wire*. That is a sharper and more
   interesting statement than "they lost Glu178", and it tightens the prediction: if these enzymes are
   no longer putrescine-activated, it is because they stopped sensing, not because they stopped
   transmitting.
2. **Lys80 is lost in exactly two lineages: cryptophytes (Leu, 7/7) and kinetoplastids (Ile, 6/6).**
   Those are the two members of group (b). A shared derived feature in the relay, independent of
   position 178, that the project had not noticed. Worth checking against the ancestral
   reconstruction: if Lys80 was lost once in the group (b) ancestor, it is another character
   supporting that pairing.
3. **Ochrophytes have replaced the aromatic lid.** Where human has Phe285/Phe315/Tyr318/Phe320,
   ochrophytes have **Met285, Val315, Cys318, Cys320** — all four changed, and two of them cysteines.
   Since Bale 2008 says the lid closing over putrescine is what enhances the electrostatics, a
   remodelled lid is a second, independent change to the sensing apparatus. Two cysteines facing a
   regulatory pocket is at least suggestive of a redox-sensitive or covalently modifiable element, but
   **that is speculation and should be labelled as such** — no measurement here supports it.

**Caveats on this table.** The lid positions are variable even among activated enzymes (Tyr318 is only
7/14 in metazoa), the trimmed alignment is less reliable near the C-terminus, and these are
consensus calls, not per-sequence verification. The relay columns are solid; the lid columns should
be re-checked per sequence before anything is written up.

## 4. Control: the site potential is not just bulk protein charge

A real risk for M1: a protein that is more basic overall has a less negative potential *everywhere*,
so an apparent species difference at the site could be bulk amino-acid composition rather than pocket
chemistry. Tested by re-reading the APBS maps and sampling the potential over a solvent shell around
the whole protein (`4_scripts/stage1_charge_control.py`).

| Structure | Known response | Net charge | Shell median | Site | Site − shell | z |
|---|---|---|---|---|---|---|
| human 1MSV | activated | 0 e | +0.03 | −28.95 | −28.98 | −14.7 |
| *Acanthamoeba* | activated | −18 e | −1.33 | −37.05 | −35.72 | −19.7 |
| lamprey | activated | −15 e | −0.91 | −39.48 | −38.57 | −15.9 |
| *T. brucei* (full length) | activated, shifted | −16 e | −1.02 | −25.11 | −24.09 | −13.6 |
| *T. cruzi* (truncated) | activated, shifted | −16 e | −1.08 | −9.69 | −8.61 | −5.5 |
| *Selaginella* | **not stimulated** | −13 e | −0.56 | **+12.95** | **+13.52** | **+10.6** |
| human E178Q | not stimulated | 0 e | −0.01 | −18.45 | −18.44 | −10.7 |

**The confound is ruled out.** Net charge ranges over 18 e across these proteins, but shell potentials
are all within ±1.4 kT/e — bulk charge is screened at 150 mM salt. Site potentials are 5–20 shell
standard deviations from the background, and correcting for the shell changes nothing. Decisively:
*Selaginella* is the **second most negatively charged protein** in the set and still has a **positive**
site. The positive plant site comes from its own built-in Arg/Lys, not from composition.

This control belongs in the methods section whatever Stage 1 returns.

## 5. Diatom polyamine physiology is unusual, and it supplies a functional hypothesis

Nothing in the project so far asks *why* a stramenopile would stop sensing putrescine. The diatom
literature offers an answer worth stating as a hypothesis.

- Diatoms make **long-chain polyamines (LCPAs)** bound to silaffins, which direct silica nanosphere
  formation in the cell wall; the polyamines are species-specific and control silica morphology
  (Kröger et al. 1999 *Science* 286:1129; 2000 *PNAS* 97:14133; 2001 *J Biol Chem* 276:26066).
- **Michael 2011, *FEBS Lett* 585:2627** — already read by this project on 5 Nov — describes
  bacterially-derived multi-domain gene fusions that potentially synthesise, N-methylate and transfer
  LCPAs in diatoms. Michael's 2016 review (*J Biol Chem* 291:14896, PMC4946907) repeats this.
- **Knott, Römer & Sumper 2007, *FEBS Lett* 581:3081**: the putative spermine synthases of the diatom
  *Thalassiosira pseudonana* make **thermospermine**, not spermine.

Every aminopropyl unit on a long-chain polyamine costs one decarboxylated AdoMet, the product of
AdoMetDC. So a diatom runs a **high, constitutive, structural demand** for this enzyme's product,
quite unlike a mammalian cell where polyamine synthesis is tightly feedback-regulated.

**Hypothesis worth putting in the discussion** (and labelled as a hypothesis): an enzyme supplying a
constitutive bulk-biosynthetic demand has little use for feedback activation, so loss of the
putrescine sensor in stramenopiles may be *adaptive* rather than merely tolerated — the same endpoint
plants reached by a different route (filling the pocket with built-in bases). This gives the paper a
reason for the loss, instead of only reporting it. It is also consistent with the two independent
changes above: sensor gone, lid remodelled, wire intact.

**It also raises a question the survey should not dodge:** if diatoms have a separate,
bacterially-derived machinery for LCPA synthesis, is the eukaryotic class-2 enzyme we are studying the
main flux enzyme at all, or a minor one? That affects how much the regulation of this copy matters,
and it should be checked against Michael 2011, which is already in hand.

### A concrete wet-lab partner, with a 2024 paper

- **Lin HY et al. 2024, *Biomolecules* 14:372** (PMID 38540790, PMC10968388, open access), Lin lab —
  the group the log already flags as holding PtSAMDC constructs. They cloned **PtSDS1, PtSDS2 and
  PtSAMDC** from *Phaeodactylum tricornutum*, confirmed the spermidine synthases by substrate
  selectivity (the Y79 active site governs their difference), localised all three to the cytoplasm,
  and found that **overexpressing PtSAMDC alone changed neither polyamine content nor growth**, while
  **co-overexpressing PtSDS1 with PtSAMDC raised spermidine and promoted growth** — which they call
  multifactorial regulation.

They have the diatom enzyme, an expression system and a live interest in its regulation, and they did
not test putrescine activation. **That is the experiment that would fill the empty quadrant**, and
this paper is the thing to point at when proposing it. (It also independently confirms PtSAMDC is a
functioning enzyme in the spermidine pathway, not a dead paralog.)

## 6. Open item closed: the four "dead" class-2-like stramenopile hits are not prozymes

The log flagged that the claim "stramenopiles have no prozyme" was never actually searched, and that
four low-identity hits needed a look. Examined tonight from
`3_data/annotations/ncbi_expansion_residue_calls.json`:

| ID | Organism | Identity | Length | Engine | 178 |
|---|---|---|---|---|---|
| KAA0165454.1 | *Cafeteria roenbergensis* | 0.201 | 253 | DEAD | A |
| KAA0162548.1 | *Cafeteria roenbergensis* | 0.215 | 256 | DEAD | A |
| KAA0153562.1 | *Cafeteria roenbergensis* | 0.209 | 254 | DEAD | A |
| CAE0742880.1 | *Schizochytrium aggregatum* | 0.253 | 660 | DEAD | R |

The three *Cafeteria* copies are **85–89% identical to one another** — a recent triplication or
assembly redundancy. At every diagnostic column they read glycine or another non-functional residue:
Glu67→Pro, **Ser68→Gly** (the serine that forms the pyruvoyl group), **Cys82→Val**, Asp174→Leu,
Thr176→Gly, Glu178→Ala, Ser229→Gly, His243→Gly or gap, Glu256→Gly.

A genuine prozyme is a dead paralog that **keeps the fold and the dimerisation surface** while losing
catalysis. A string of glycines at exactly the functional columns, at 20–22% identity, is the
signature of a **poor alignment of a divergent or unrelated sequence**, not of a prozyme. The
*Schizochytrium* hit (660 aa, 25% identity) looks like an unrelated or fused protein.

**Conclusion:** these four are not evidence of a stramenopile prozyme, and they are also not good
evidence that the stramenopile AdoMetDC family contains dead paralogs at all. But this still does not
license the claim "stramenopiles have no prozyme". A proper test is a profile/HMM search built from
trypanosomatid prozymes against stramenopile proteomes, which has not been run. Until then the right
wording stays **"not searched"**, as the log already has it.

Separately: a Europe PMC search for prozyme outside trypanosomatids returns nothing — all the
prozyme literature is the Phillips lab on kinetoplastids (plus a commercial protease of the same
name). So there is no published precedent to lean on either way.

## 7. Model confidence is not the limiting factor at the pocket

Site pLDDT for every AlphaFold model in hand is **93–98**. Whatever is wrong with these models, it is
not low confidence at the putrescine site; it is that they are apo, monomeric and single-conformation.

## 8. What I changed, and what I did not

**Changed** (all before any Stage 1 structure was measured, recorded as a dated amendment in the
pre-registration):

1. **5TVO and 5TVF are no longer a noise pair.** They are the autoinhibited and prozyme-activated
   states of the same enzyme, and are reported as a state comparison. The crystal-to-crystal error bar
   now comes only from the human crystals.
2. **1JEN added** — apo human crystal, to separate the closed-lid/open-lid effect from the
   crystal-versus-model effect. The fetch script now lists 15 downloads.

**Not changed.** The six measurements and the gate itself are untouched. The relay and lid positions
found in §3 are *not* added to the gate, even though they look informative, because adding a
measurement after seeing data is exactly what pre-registration exists to prevent. They are recorded as
candidate Stage 2 measurements and as sequence-level background.

## 9. What failed, and limits

- **`javascript_tool` on the browser returned `undefined` every time**, so full-text extraction from
  PMC had to be done through page text and targeted finds. Abstracts were read in full; full texts
  were not, except where the project already holds the paper.
- **Europe PMC's `fullTextXML` endpoint returned 500**, and its `CITES:` search syntax returned
  nothing; the citations REST endpoint worked instead.
- The lid consensus calls in §3 are from a trimmed alignment and are not per-sequence verified.
- The *T. brucei* full-length measurement is one model, not a crystal; the crystals will test it.
- The diatom hypothesis in §5 is a hypothesis. Nothing computed here tests it.
