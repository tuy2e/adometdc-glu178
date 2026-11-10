# Tree robustness round 2 (7 Nov 2026)
Files: /root/audit/strameno/expand2/
- Alternatives built: H_li (stramenopiles + Asp copy grouped with opisthokonts, as in Li 2019 Fig. 2; KE+CR kept together), H_sar (stramenopiles + Asp copy grouped with the alveolate Gln/apicomplexan clade; KE+CR together). Constraint trees list only constrained taxa, plus green algae/plants as an anchor group. Start trees were made by regrafting in the ML tree; constrained search used -nstop 20.
- Mixture model: LG+C60+F+G, PMSF (site frequencies from ML tree expA), UFBoot 1000 + SH-aLRT 1000 -> pmsfA.treefile, logL -85574.974.
  - Lineage ST+DP+KE+CR: SH 87.4 / UFB 90; ST+DP 62.3/87; KE+CR 99/100; DP 99.8/100; ST 78.5/90; CR 88/99; KE 98.6/100.
  - (Q.pfam+R7 ML: lineage 93.5/88; ST+DP 7.2/56; KE+CR 99.4/100.)
AU tests (10,000 RELL, seed 1)
| Tree | Q.pfam+R7 (au7) | no pocket, Q.pfam (au7NP) | C60 PMSF (au10C60) |
| ML Q.pfam (expA) | 0.936 | 0.791 | 0.431 |
| H_crypto | 0.0328 | 0.0684 | 0.00783 |
| H_kin | 0.000474 | 0.0208 | 0.000165 |
| H_dino | 0.0042 | 0.00418 | 0.0131 |
| H_all | 1.19e-05 | 0.00492 | 0.000451 |
| H_li (Q.pfam-optimized) | 0.183 | 0.482 | 0.0661 |
| H_sar (Q.pfam-optimized) | 0.182 | 0.191 | 0.112 |
| ML C60 (pmsfA) | - | - | 0.759 |
| P_li (C60-optimized) | - | - | 0.373 |
| P_sar (C60-optimized) | - | - | 0.391 |
ΔlnL: Q.pfam H_li 13.4, H_sar 11.3; C60 P_li 8.4, P_sar 8.2; no-pocket H_li 2.7, H_sar 6.4.
Residues by group (exp_trim columns for human 109 / 178):
ST(45): 109 W36 K5 F2 Y2; 178 N31 H11 D2 Q1 | DP(16): W16; D16 | KE(7): W6 Y1; S6 N1 | CR(7): W5 Y2; N7 | rest(112): 109 mostly R/D/Y/S/F; 178 E98 Q10 D2 gap2

## Ancestral reconstruction on the two-group trees (Q.pfam+R7, fixed trees, -asr) — 8 Nov
Position 178 (human Glu178):
- ML tree (expA): one E->N at the 4-group lineage base (E 0.999 -> N 0.991), then N->D in the dino/Perkinsus copy, N->H in oomycetes, N->S in kinetoplastids. 0 reversals to E.
- H_li tree: TWO losses: E->D(0.46) at base of stramenopiles+Asp copy (61), then D->N(0.97) at stramenopile base; and E->N(0.99) at base of kinetoplastids+euglenid+cryptophytes (14). 0 reversals.
- H_sar tree: TWO losses: E(0.78)->N(0.44) at stramenopiles+Asp copy base; E(0.78)->N(0.99) at kineto+crypto base. 0 reversals.
Position 109 (human Ser109; Trp in the Glu-less groups):
- ML tree: one H(0.78)->W(0.99) at the lineage base; later W->Y/K in a few subgroups and single tips.
- H_li: THREE separate gains of W: stramenopiles+Asp copy (Y0.70->W0.91), cryptophytes (Y0.64->W0.86), kinetoplastids (Y0.62->W0.99). Ancestral states uncertain.
- H_sar: ambiguous (one deep low-probability W 0.60 on a 123-tip node with a reversal to Y in 48 tips).
Interpretation: the two-group alternative needs parallel losses of Glu178 (2) and parallel gains of Trp109 (2–3); the one-lineage tree needs one of each. Parsimony favours one lineage, but these are linked pocket residues, so convergence is plausible and AU does not reject the two-group trees.

## Ancestral reconstruction on the two-group trees (9 Nov)
Q.pfam+R7, fixed topology (-te), IQ-TREE 2.0.7 (-asr not supported with PMSF/C60 in this version, so C60-optimized topologies were reconstructed under Q.pfam+R7).
Glu178 (column 178):
- ML tree expA: one E->N at the base of the 4-group lineage (from earlier run; P 0.99).
- H_li (Q.pfam-opt): E->D(0.46) at base of ST+DP (61 tips), D->N(0.97) at stramenopiles; E->N(0.99) at base of KE+CR (14 tips) => 2 losses.
- H_sar: E->N(0.44) at ST+DP base, N->D(0.97) for DP; E->N(0.99) at KE+CR => 2 losses.
- P_li / P_sar (C60-opt): E(0.81-0.83)->N(0.42-0.44) at ST+DP; E->N(0.99) at KE+CR => 2 losses.
- All trees: Glu never regained; dino alveolate copy E->Q once (9 tips).
Trp109 (column 109):
- ML expA: one H(0.78)->W(0.99) at the lineage base (UFB 88), with later W->Y/K reversals in a few tips/subclades.
- H_li: three gains (ST+DP 0.91; cryptophytes 0.86; kinetoplastids 0.99).
- H_sar, P_li, P_sar: a single low-confidence W gain (0.60-0.65) placed at a 123-tip node that also includes plants, red algae and haptophytes (which then revert to Y) => unstable; not interpretable.
Smallest split containing all 4 groups: ML 75 (exactly the 4 groups); P_li/P_sar/H_sar 123 (extra: haptophytes 12, red algae 12, green 20, 2 ODC-fusion dinos, 1 metamonad); H_li 133.
