# Putrescine docking into AdoMetDC pockets (10 Nov 2026, to-do 3B/3C)

**Short version**

- **The method passes its controls.**
  - Putrescine redocked into the human crystal (1MSV) lands in the crystal position (RMSD 0.7 Å).
  - It also docks into the human-type site of AlphaFold models of two putrescine-activated enzymes, *Acanthamoeba* and lamprey.
  - It does **not** dock into the site of a plant, whose pocket is filled by built-in residues (about 5 Å³ of space). That matches the potato crystal (Bennett 2002).
  - In *T. cruzi* it docks **6 Å away** from the human site, which fits Clyne 2002's "shifted site".
- **None of the Glu-less test enzymes binds like the activated controls.**
  - Diatom, brown alga, cryptophyte and the *Perkinsus* Asp copy: 0 poses in the human site. Their best poses sit 4–6 Å away, as in *T. cruzi*.
  - Two of three oomycetes have nearly filled pockets (2–10 Å³), like the plant.
  - *Cafeteria* and the third oomycete dock in the site only occasionally (2 of 9).
- **What this means:** the Glu-less pockets differ from the activated human-type pocket in shape, as models. But docking **can't tell** "not activated" from "activated through a shifted site". *T. cruzi* shows the second can happen.
  - The oomycete result supports their separate, plant-like hypothesis.
  - The stramenopile prediction stays **plausible, not proven**. Only a measurement will settle it.
- **Not usable:** the charge-aware AutoDock4 scoring **failed its control** (it put human putrescine 9 Å from the crystal site), so it is not used. Docking also can't test the human E178Q effect: E178Q docks exactly like wild type, because Glu178 changes signalling, not binding.

## Method

- **Structures:** human 1MSV chain A, the S68A proenzyme with putrescine bound (crystal), plus an in-silico E178Q copy (Glu178 renamed Gln). Then 13 AlphaFold models from the FOLDTRACE run (single-chain proenzymes, as 1MSV is), read from the copy in Downloads and left unchanged.
- **Superposition:** onto 1MSV, using Cα atoms of human residues within 18 Å of putrescine, mapped by pairwise MAFFT and refined twice to atoms under 2 Å. That gave 103–145 core Cα, RMSD 0.6–0.9 Å, mean site pLDDT 93–97.
- **Receptors:** Open Babel, hydrogens at pH 7.4, Gasteiger charges.
- **Ligand:** putrescine as the dication `[NH3+]CCCC[NH3+]` (RDKit 3D).
- **Docking:** AutoDock Vina 1.2.5, exhaustiveness 32, 9 modes.
  - Wide box: 22 Å centred on the crystal putrescine, 3 seeds, 27 poses per structure.
  - Focused box: 12 Å, 1 seed, 9 poses.
- **"In site":** the pose's heavy-atom centre is within 2.5 Å of the crystal putrescine centre.
- **Cavity:** free volume (probe 3 Å) within 4 Å of the crystal putrescine centre.
- **AutoDock4 scoring:** charge-aware, run through Vina's AD4 mode with AutoGrid4 maps. It failed the human control, so it is reported but not used.

## Results

| Structure | 174 | 178 | Wide box: poses in site | Focused box: best Vina score (distance from site) | Focused box: poses in site | Free volume at site (Å³) |
|---|---|---|---|---|---|---|
| **Human 1MSV crystal** | D174 | E178 | 17/27 | −4.42 (0.1 Å) | 7/9 | 33.9 |
| Human E178Q (in silico) | D174 | Q178 | 17/27 | −4.44 (0.1 Å) | 7/9 | 33.9 |
| **Lamprey** (activated-type, model) | D | E | 13/27 | −4.04 (1.4 Å) | 6/9 | 25.5 |
| ***Acanthamoeba*** (activated, model) | D | E | 22/27 | −4.04 (0.5 Å) | 8/9 | 42.0 |
| ***T. cruzi*** (activated via a shifted site) | D | S | 0/27 | −3.65 (5.9 Å) | 0/9 | 70.2* |
| **Plant *Selaginella*** (built-in, not stimulated) | G | E | 0/27 | −2.82 (4.9 Å) | 0/9 | 4.6 |
| Diatom *Phaeodactylum* | D | N | 0/27 | −3.11 (6.1 Å) | 0/9 | 50.5 |
| Brown alga *Ectocarpus* | D | N | 0/27 | −4.51 (4.6 Å) | 0/9 | 32.8 |
| *Cafeteria* (Bigyra) | D | N | 5/27 | −3.85 (2.2 Å) | 2/9 | 32.6 |
| Oomycete *Hyaloperonospora* | E | H | 0/27 | −3.58 (6.2 Å) | 0/9 | 9.6 |
| Oomycete *Saprolegnia* (Arg111, Lys109) | E | H | 0/27 | −3.19 (4.5 Å) | 0/9 | 1.8 |
| Oomycete *Aphanomyces* (Arg111, Lys109) | E | H | 8/27 | −3.98 (2.4 Å) | 2/9 | 30.5 |
| Cryptophyte *Guillardia* | D | N | 0/27 | −3.87 (4.2 Å) | 0/9 | 14.9 |
| *Perkinsus* Asp copy | D | D | 0/27 | −4.08 (4.8 Å) | 0/9 | 14.6 |
| Haptophyte *Emiliania* (keeps Glu178) | E | E | 0/27 | −3.19 (5.0 Å) | 0/9 | 1.0 |

\* The *T. cruzi* model lacks its first ~20 residues (no positions 13/15), so its site looks emptier than it is.

- **Contacts of the in-site poses:** human Glu15, Asp174, Thr176 (and Glu256 in some runs), as in the crystal. *Acanthamoeba*: the same four. *Aphanomyces*: Asn13, Glu15, **Arg111**, Ser113, Thr176. In that one, Arg111 sits where potato Arg114 does.
- **The brown alga's best pose** (−4.51, the best score of any structure) sits 4.6 Å from the human site, against Glu15, Asp174 and Thr176 equivalents. That is a shifted site, as in *T. cruzi*, still using Asp174.

## How to read it

1. **Activated, Glu178-keeping enzymes dock in the human site; the Glu-less test enzymes mostly don't.** That fits their pockets being remodelled, not just missing one residue.
2. **Shifted site or filled pocket.** Diatom, brown alga, cryptophyte and the *Perkinsus* copy look like *T. cruzi* (a shifted site). Two of three oomycetes look like the plant (a filled pocket).
   - *T. cruzi* is still activated through its shifted site, so "shifted" doesn't mean "not activated".
   - The plant is constitutively active, so the oomycete result fits a plant-like, built-in mechanism.
3. **The haptophyte keeps Glu178 but its model pocket is filled** (Arg109, Lys13). Either haptophytes have their own built-in variant, or this apo model packs side chains into the cavity. Either way, it is a reminder that a model with nothing bound can close a pocket. *(Haptophyte responses have never been measured.)*
4. **Docking doesn't test signalling.** Human E178Q docks exactly like wild type. The human mutant's loss of stimulation is about transmitting the signal, not binding, and no static docking can show that.

## Limits

- The models are AlphaFold predictions with nothing bound; side chains are not relaxed and only one conformation is used.
- Vina has no explicit electrostatics. The charge-aware alternative failed its control here.
- Score differences under about 1 kcal/mol aren't meaningful.
- Proenzyme models were used; processing could reshape the site.
- The *T. cruzi* model is truncated.

## Files

`8_full_workspace_backup/workspace_part5_docking.zip` holds the aligned receptors, all poses, logs, scripts and summary JSON.
