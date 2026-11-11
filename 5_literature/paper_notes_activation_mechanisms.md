# Paper notes: how AdoMetDC is switched on in trypanosomes, plants and *Plasmodium* (10 Nov 2026)

These are notes from a full read of 5 PDFs (Clyne 2002, Beswick 2006, Bennett 2002, Wells 2006, Williams 2011). The read was done in a separate chat, on Tui's copies of the PDFs, and checked here against our alignment. The figure and table numbers are as reported in that read; reopen the PDFs before quoting any number in the paper.

## Clyne et al. 2002, *Biochemistry* 41:13207 (*T. cruzi*)

- Wild type: Kd(putrescine) 150 µM (human 5–6 µM); stimulated 11-fold (Table 2).
- **D174V:** no binding (Fig 5F) and no activation (0.8-fold). **Asp174 is the only binding residue shared with human.**
- E15A: binds tighter (Kd 8 µM) but activation drops to 2-fold. E256A: binds normally, activation 4-fold. So these carry the signal, not the ligand.
- **I80K/S178E** (rebuilding the human Lys80/Glu178 pair): no effect (12-fold vs 11-fold). **Positions 80 and 178 do nothing in *T. cruzi*.**
- **R13L:** the proenzyme doesn't process, and putrescine can't rescue it (Fig 2c). **Arg13 is a processing determinant, not an activation residue.** Arg13, Gln176 and Met113 probably block the human-type site, so putrescine binds a shifted, more surface-exposed site.

## Beswick et al. 2006, *Biochemistry* 45:7797 (*T. cruzi*)

- Wild type: 29-fold activation, Kact 0.42 mM, which matches trypanosome putrescine levels (Table 1).
- D174V: Kact > 45 mM, no activation. S111R: 3-fold. **F285H: 360-fold**, so one change can raise responsiveness enormously.
- Active-site F7A and T245A keep a normal Kact but drop activation to 4-fold (Table 2), so activation also depends on coupling residues outside the pocket.
- The site is promiscuous: cadaverine, diaminopropane and CGP 48664A also activate it.
- **Assay caveat:** contaminants in commercial AdoMet masked activation, and purified human enzyme carries about 0.5 mol putrescine per mol. Published "not activated" results can be artefacts.

## Bennett et al. 2002, *Biochemistry* 41:14509 (potato, PDB 1MHM)

- There is no room for putrescine. **Arg18 (human 13), Arg114 (human 111) and His294 (human 285)** fill the cavity, and **Asp174 becomes Val181**. Glu178 is kept as Glu185.
- The enzyme is constitutively as active as putrescine-activated human enzyme, and it is a monomer.
- Bennett's warning: the Arg roles "were not predicted by examining the sequence alone, nor ... by homology modeling".

## Wells et al. 2006, *J Mol Graph Model* 24:307; Williams et al. 2011, *Mol Biochem Parasitol* 180:17 (*Plasmodium*)

- *Plasmodium* **keeps** the human relay residues Glu11, Lys80, **Glu178 (Glu219)** and Glu256, but **Asp174 becomes Lys215**. The 11-RVVIK-15 motif (Arg11 = human 13, Lys15 = human 17) plus Lys215 fill the putrescine site.
- Mutagenesis: R11L loses about 90% of activity and doesn't process; K15A and K215A each lose about 50% (Fig 8).
- Williams: monofunctional *Pf*AdoMetDC is not stimulated by putrescine and doesn't need it to process. It is slow (kcat 5.3 min⁻¹, Table 2). The proposed new subclass is 2b-III.

## Checked against our 200-enzyme alignment (human numbering)

| Group (n) | 13 | 15 | 80 | 111 | 113 | 174 | 176 | 178 | 256 | 285 |
|---|---|---|---|---|---|---|---|---|---|---|
| human | L | E | K | F | S | D | T | E | E | F |
| **Ochrophyta (22)** | T19 | E20 | **K22** | **G22** | S19 | **D16 E4** | V13 T8 | N22 | E22 | M14 F7 |
| Bigyra (11) | I8 | E7 Q3 | K11 | L4 G4 | S8 | E5 S3 D2 | F4 T3 N3 | N9 D2 | E11 | F9 |
| **Oomycota (12)** | N12 | E12 | K12 | **R5** S4 N2 | S11 | E11 D1 | T12 | **H11** Q1 | E12 | Y6 S5 |
| Dino/*Perkinsus* Asp copy (15) | R15 | E15 | K15 | C11 | S15 | D13 E2 | T5 Q5 | D15 | E14 | F10 |
| Cryptophytes (7) | K6 | E7 | L7 | Q7 | S7 | D5 | M4 S2 | N7 | E7 | **H5** |
| Kinetoplastids (6) = *T. cruzi* | R | E | I | S | M | D | Q | S | E | F |
| Euglenids (5) | R5 | E3 Q2 | I4 | Q4 | S5 | D5 | Q5 | N4 K1 | E5 | L4 |
| Diplonemids (8) | K7 | E8 | K8 | S8 | S8 | D8 | T8 | E8 | E8 | F7 |
| Land plants (9) | R8 | E7 | K9 | **K6 R3** | T8 | **V5** | T8 | E8 | E9 | **H6** |
| Apicomplexans + chromerids (5) | R4 | mixed | K5 | L4 | mixed | **not D** | mixed | E5 | E5 | mixed |

- **The papers' mappings check out in our data.**
  - Plants: R13, K/R111, V174, H285, Glu178 kept (Bennett).
  - Apicomplexans: R13, no Asp174, Glu178 kept (Wells).
  - Kinetoplastids: R13, I80, S111, M113, Q176, S178 (Clyne).
- **Stramenopiles keep Asp174.** In Ochrophyta 20 of 22 are acidic there (16 Asp, 4 Glu); oomycetes are Glu. They also keep Glu15, Glu256 and Lys80.
- **The ochrophyte pocket is human-like at the positions that make *T. cruzi* different.** It has no Arg13, no Gln176, no Met113 and no Ile80. Where ochrophytes differ from human is 178 (Asn), 111 (Gly, so the Phe111 stacking residue is lost) and 109 (Trp, which sits about where Phe111 would).
- **Oomycetes look plant-like in part:** Arg111 in 5 of 12 (potato Arg114 is at the same position), His178, and Tyr/Ser285.
- **Cryptophytes have His285 in 5 of 7** (potato's His294 is at the same position) and Leu80.

## What this changes

1. **Glu178 alone is not the predictor, and the paper must say so.** *T. cruzi* is activated without it; *Plasmodium* has it and isn't activated. Our prediction has to rest on the whole pocket, not one residue.
2. **The stramenopile prediction weakens from "likely" to "plausible".**
   - In favour: the ochrophyte pocket matches human everywhere *T. cruzi* has rewired, so the closest precedent is the human E178Q mutant (more active without putrescine, not stimulated), not *T. cruzi*.
   - Against: ochrophytes keep Asp174, the residue *T. cruzi* binds with, so putrescine binding may survive. Gly111/Trp109 is a change no studied enzyme has. And single substitutions can swing responsiveness widely (F285H, 360-fold).
   - **Only a measurement settles it.**
3. **Oomycetes should be treated separately.** Arg111 + His178 look like plant-type built-in stand-ins for putrescine. If so, the oomycete enzyme would be constitutively active by a different route.
4. **Position 13 is a processing determinant** (*T. cruzi* R13, potato R18, *Pf* R11), not an activation residue. The earlier reason for calling cryptophytes "uncertain" (ancestral Arg13) is replaced. They remain uncertain because of His285, Met176 and Leu80.
5. ***Rhabdomonas* Lys178 is probably not the plant mechanism.** Plants keep Glu178 and add bases at 13, 111 and 285. Treat it as an open curiosity.
6. **Say which "activation"**: stimulation of processing, of catalysis, or both. *T. cruzi* has catalysis stimulated but not processing; *Plasmodium* has neither.
7. **Literature "not activated" negatives** can be assay artefacts (Beswick). Use them with care.
8. **Docking (to-do 3B/3C) now matters more.** Putrescine docked into the diatom and oomycete models, alongside 1MSV (human) and 1MHM (potato), is the cheapest computational test of point 2. **Done 10 Nov** (`3_data/test_results/RESULTS_docking.md`; the plant control was *Selaginella* rather than 1MHM). Result: most Glu-less pockets don't take putrescine in the human position; they are shifted (like *T. cruzi*) or filled (oomycetes, like the plant). That supports point 3 but can't settle point 2.
