# Which "WOW" finding is the best publication bet? — honest ranked assessment
*Tuy, 4 Nov 2026. You asked me to look past the prozyme angle and see which of the other FOLDTRACE findings could lead to a paper. I re-read your `WOW_findings.md`, did first-pass novelty checks on PubMed, and verified the identity of the key proteins on UniProt. Below is a straight ranking — including a meta-warning you should know.*

## Meta-warning first (this matters for all of them)
AdoMetDC is a **heavily worked enzyme by a small number of labs** — Phillips and **Anthony Michael** (both at UT Southwestern, the same place), and Ealick (Cornell). They've already published comparative/evolutionary work on it. So *any* "surprising AdoMetDC pattern" carries real pre-emption risk, and every candidate below has to be checked against their papers specifically. That's not a reason to give up — it's a reason to novelty-check hard before investing, which is exactly what sank the over-claim last time.

---

## #1 (best bet) — Diatom + brown-algae canonical AdoMetDC is "stuck-on" (a second, independent route to putrescine-independence)
**The finding (your WOW part 2b):** 6 of 8 diatoms and all 3 brown algae that passed the fold gate share the *same* pocket rewiring — **Asn at 178, Gly at Phe111, Thr at Leu13** — and *all 15 proteins* with that exact combination are **Ochrophyta** (the stramenopile group = diatoms + brown algae). Plants reach putrescine-independence a *different* way (built-in Arg/Lys); Ochrophyta did it by the **Glu178→Asn "stuck-on" swap**. Two independent evolutionary solutions to the same problem.

**Why it's the strongest candidate:**
- **Crisp, testable prediction with an existing precedent.** Human **E178Q is already known to be "stuck-on"** — ~4× more active *without* putrescine (Bale et al. 2008). So the prediction "diatom AdoMetDC is constitutively active / putrescine-independent" is sharp and has a direct human-mutant analog. That's a clean experiment.
- **Genetically tractable model organism.** The diatom protein is **Phaeodactylum tricornutum B7FPB8** (322 aa, gene PHATRDRAFT_3362) — *Phaeodactylum* is a workhorse model alga with transformation tools, so the lab test is realistic and Marc can find a collaborator.
- **Compelling functional story.** Diatoms pour polyamines (long-chain polyamines) into their **silica cell walls** — they need huge aminopropyl/dcSAM flux. A constitutively-on AdoMetDC fits that demand. The "why" writes itself.
- **FOLDTRACE-native.** It's literally a regulatory-site comparative map across a lineage — exactly what the method is for.

**Novelty status (checked, honest):** In all of PubMed, **only ONE** paper links AdoMetDC to diatoms/brown algae/stramenopiles: **Michael 2011, FEBS Lett (PMID 21827754)**. Crucially, that paper is about a **different protein** — giant *multidomain fusion machines* (AdoMetDC fused to aminopropyltransferase ± methyltransferase) that build long-chain polyamines for silica, and it even notes some carry a *dead* AdoMetDC domain. **Your WOW proteins are the separate, canonical standalone AdoMetDC** (B7FPB8 322 aa; Ectocarpus D7FY29 388 aa, gene "SAMDC"; Fragilariopsis A0A1E7FKD8 369 aa) — I confirmed this on UniProt. So the **regulation of the housekeeping diatom enzyme (the stuck-on pocket) appears unexplored.** ⚠️ But Michael 2011 is paywalled and sits right next door to this topic — you must read its full text (ILL it) and frame your novelty carefully against it before committing.

**Confirmation path:** (1) re-pull the exact Ochrophyta pocket list from your run and confirm the convergence is clean + monophyletic (one origin in the stramenopile lineage); (2) model/overlay a diatom pocket to show it's structurally "filled/on"; (3) **lab:** express *Phaeodactylum* AdoMetDC, measure activity ± putrescine → does it behave like human E178Q? That's the paper.

---

## #2 — The Bodo/kinetoplastid prozyme (your current angle) — NOT a dead end
Honestly it got **stronger** with the interface mapping (the dead copy keeps M150/Y154, the active enzyme keeps the P37 switch). Its weakness is novelty: the "regulatory pseudoenzyme predates parasitism in free-living *Bodo*" theme is **partly pre-empted** (the PRMT prozyme was already shown in Bodo). Best home: **co-authorship on Marc's FOLDTRACE paper** as a case study, or a short note. Keep it alive; it's your safest co-authorship, but #1 is a cleaner shot at a *first/independent* contribution.

## #3 — Eukaryotic-type AdoMetDC hiding in specific bacterial lineages (your WOW part 7)
54 bacteria (notably **Bdellovibrionota** [predatory bacteria] and **Shewanella**) carry a full **eukaryote-sized, eukaryote-type (class-2) AdoMetDC with an intact engine** — unusual for bacteria, which normally use the small class-1 enzyme. Possible horizontal gene transfer / unappreciated distribution. **Novelty looks open** (my search found nothing on this). **But** it's a distribution/phylogenetics story with **no crisp functional test** and a murkier "so what" — harder for an undergrad to turn into a lab-backed paper. Good as a *bioinformatics note* or a section of a bigger paper, not a first solo project.

## #4 / lower — one-off curiosities (keep as color, not projects)
- **Palpitomonas** independently lost the AdoMetDC engine (non-kinetoplastid) — a neat convergence, but a single predicted structure; risky to build on.
- **Shared rare trio** (His13+Glu66+Asp256) in an unrelated bacterium + archaeon → possible HGT; too niche to stand alone.
- **Invariant Gly83/Pro246** near the engine → descriptive, supporting-figure material.

---

## Bottom line / recommendation
Pursue **#1 (the diatom "stuck-on" convergence)** as your lead independent project, and keep **#2 (Bodo prozyme)** running as the co-authorship track. #1 wins because it has the crispest testable prediction (human E178Q precedent), a tractable model organism (*Phaeodactylum*), a clear functional rationale (silica polyamines), and — on first check — cleaner novelty, as long as you read and carefully distinguish Michael 2011.

**Immediate next steps I can do now (no lab needed):**
1. Re-pull the exact list of Ochrophyta "stuck-on" proteins from your run and build a clean cross-lineage pocket table + a tree, to show the convergence is real and single-origin.
2. Map the pocket residues the way I did for the prozyme, and model a diatom AdoMetDC pocket to show the "on" conformation.
3. Draft the one-paragraph novelty framing vs Michael 2011 so you can sanity-check it with Marc.

**For your ILL pile:** Michael AJ 2011, *FEBS Lett* 585:2627 (PMID 21827754) — the one adjacent diatom-AdoMetDC paper; you need its full text to lock the novelty of #1.
