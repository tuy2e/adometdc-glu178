# p22 ColabFold heterodimers: fused ODC half + same-species standalone ODC

Each `*_COLAB.txt` is `FUSION:STANDALONE` on one line (chain A = fusion, chain B = standalone). Settings as in `data/p16_colabfold/README.md` (alphafold2_multimer_v3, unpaired_paired, 4 seeds, save_all). Fusion sequences are from Marc's models (read only), standalones from the UniProt records cached in `data/p17_standalone/`.

| job | chain A (fusion) | aa | chain B (standalone) | aa | total residues | file |
|---|---|---|---|---|---|---|
| a | A0A9C7PQ51 | 762 | A0A9C7PSP6 | 540 | **1302** | `a_Galdieria_partita_A0A9C7PQ51_fusion+A0A9C7PSP6_standalone_COLAB.txt` |
| b | A0A2V3IZY9 | 730 | A0A2V3IYX3 | 433 | **1163** | `b_Gracilariopsis_chorda_A0A2V3IZY9_fusion+A0A2V3IYX3_standalone_COLAB.txt` |

## (a) Galdieria partita: A0A9C7PQ51 + A0A9C7PSP6

- A0A9C7PQ51: ODC half degraded (p9 gated: 1 kept / 7 lost / 5 unresolved; lost K69R H197S E274N P275V G276D R277S C360G).
- A0A9C7PSP6: ODC-family standalone, called intact in p19 by the >= 11/13 rule. It is **not 13/13**: by structure 11 kept / 1 lost (Y331S) / 1 unresolved (p10); Y331 is outside p21b's core list.
- The question: does the fused ODC half pair with the standalone ODC at all (inter-chain PAE between the fusion's ODC half and chain B), and in the ODC dimer orientation?

## (b) Gracilariopsis chorda: A0A2V3IZY9 + A0A2V3IYX3

**The test: whether K69 of the standalone and C360 of the fusion come together the way they do in human ODC 1D7K.** The two halves are broken at opposite positions, and in a real ODC dimer those two residues come from different chains (each active site takes Lys69 from one subunit and Cys360 from the other; in 1D7K, K69 of one chain to C360 of the other is 8.7 / 8.5 A CA-CA, p13). So one of the two active sites of this heterodimer could be complete and the other doubly broken.

| residue (human ODC numbering) | fusion A0A2V3IZY9 | standalone A0A2V3IYX3 |
|---|---|---|
| K69 | **Arg429** (lost; p9 structure, pLDDT 86) | **Lys97** (kept; p17 sequence alignment) |
| C360 | **Cys680** (kept, but 3.5 A off in p9, close to the 4 A limit) | **Arg380** (lost; p17 sequence, p10 structure) |

What to measure: CA-CA distance Lys97(B) to Cys680(A), and Arg429(A) to Arg380(B), against 8.5-8.7 A in 1D7K, plus the inter-chain PAE between the two ODC parts. The standalone positions come from a sequence alignment (p17), not a structural one; check them on the predicted model before measuring.

Caveat: the fusion's ODC half has other losses too (H197E, E274V, P275V, R277A, p9), so a complete K69 + C360 pair would not by itself make a working active site.
