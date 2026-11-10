# 💡 Claude's suggestions (NOT Marc's protocol, NOT in the project folder)

These results live only in this chat. Nothing here was saved to your computer.

## Suggestion 1: Are the "extra AdoMetDC copies" in animals real, and are any dead?

**Short answer: no. None of them is a real second copy, and none is dead.**

To check, I downloaded each neighbor protein, aligned it to human and read the engine positions. This was a sequence check, not FOLDTRACE. Each case turned out to be one of these:

| Species | What the "second copy" really is |
|---|---|
| *Papilio machaon* (swallowtail butterfly), *Orchesella cincta* (springtail), *Dendrothele bispora* (a fungus, Agaricales), *Dictyocaulus viviparus* (bovine lungworm), *Apostichopus japonicus* (Japanese sea cucumber) | One gene that the genome annotation **split into two pieces**. One piece covers the start of the protein and the other covers the rest. |
| *Xenopus tropicalis* (tropical clawed frog), *Paramuricea clavata* (a soft coral, Octocorallia) | **The same protein listed twice** (100% identical) |
| Drosophila pseudoobscura | **Versions of the same gene** (97% identical) |
| *Thamnophis sirtalis* (a snake) | Versions of the same gene at the same spot |

So the only **real** side-by-side pairs, with one working and one dead copy, are in the kinetoplastids. This also corrects **Finding C in the README I saved earlier**: those animal cases should not be called second copies.

## Suggestion 2: A better tree for the kinetoplastids (`kinetoplastid_prozyme_tree.png`)

**Proteins:** 61 sequences in total.
- 54 unique kinetoplastid AdoMetDC and prozyme sequences from UniProt (Leishmania, Trypanosoma, Leptomonas, Strigomonas, Angomonas, Novymonas, Porcisia and Bodo).
- My hits from Neobodo, Palpitomonas, Eutreptiella and Crithidia.
- Human, fission yeast and a plant protein for comparison.

**Method:**
- Aligned with MAFFT and trimmed with ClipKIT.
- Built the tree with IQ-TREE (best model JTT+R5) and 1,000 bootstrap replicates.
- Read the engine residues from the alignment.

**What it shows:**

1. **The family splits cleanly into two groups.** One group has the working engine (support 96). The other is the prozymes, with a broken engine (support 100). Of the 27 trypanosomatid species in this set, 24 have one copy in each group. *Angomonas deanei* has only the prozyme in the databases, and *Leishmania tarentolae* and *Crithidia acanthocephali* have only the working copy. The other copy may simply be missing from the databases.
2. **The prozyme branch includes my WOW proteins.**
   - Bodo's dead copy sits at the base of the trypanosomatid prozymes (99).
   - Neobodo sits just below it (100).
   - So both are prozymes, not separate accidents.
3. **They all share the same "dead" signature.** 27 of the 28 prozymes, Bodo and Neobodo included, have **Arg where human has His243**; the exception is *Angomonas*, which has Lys. 19 of the 28 also have **Ser at Glu8**. A change shared like this points to one loss event in a common ancestor.
4. **Palpitomonas is different.** It falls outside both groups and has **Ser at 243, not Arg**. It lost its engine on its own, separately from the prozymes.
5. **When the duplication happened.** Eutreptiella, a euglenid and a more distant relative, has only a working copy, and it sits outside the kinetoplastid split. So the gene duplicated **after the euglenids branched off but before Neobodo and Bodo branched off from the trypanosomatids**. In other words, the prozyme existed in the common ancestor of Neobodo, Bodo and the trypanosomatids. (Their groups were checked in NCBI Taxonomy: all are in Kinetoplastea, subclass Metakinetoplastina.)
6. **One loose end.** Angomonas deanei only has a prozyme in UniProt. Its working copy is probably missing from the database.

**Caveats:**
- This is a single-gene tree.
- The prozyme branches are long, meaning those sequences changed fast.
- Engine positions here come from a sequence alignment, which can be off by a residue or two in poorly aligned regions. FOLDTRACE uses 3D structure instead.

## Files (in the chat only)

- `kinetoplastid_prozyme_tree.png`: the tree picture
- `kinetoplastid_tree.nwk`: the tree file, with support values
- `kinetoplastid_engine_residues.tsv`: engine residues for every sequence
