# BLAST searches A–H (NCBI web BLAST, 10–11 Nov 2026)

NCBI deletes web BLAST results about a day after a search runs. The permanent record is the reports in this folder, made with **Download All → Text** and named `<RID>-Alignment.txt` (or similar, starting with the RID). Use the table below to see which search each RID is.

Accession numbers don't expire, so any hit's full sequence can be fetched again from NCBI by its accession when a tree needs it.

Settings, hit counts and checks below were read from each live results page (Search Summary) on 10–11 Nov.

| Search | RID | Job title | Program / database | Taxon limit | Query | Max targets / E | Sequences hit | Status |
|---|---|---|---|---|---|---|---|---|
| A | F5G9RPR9016 | A Cryptophyte SpdSyn | blastp / nr | Cryptophyceae (3027) | *Guillardia* SpdSyn (L1IR35), 305 aa | 500 / 1e-20 | 4 | ✓ |
| A2 | F5J4PJTH014 | A2 Cryptophyte SpdSyn TSA | tblastn / TSA (26 projects) | Cryptophyceae (3027) | *Guillardia* SpdSyn (L1IR35) | 500 / 1e-20 | 64 (19 species) | ✓ Includes 5 *Roombia* hits. *Roombia* is a **katablepharid** (NCBI files katablepharids under Cryptophyceae); label separately. |
| B | F5J6K5RN014 | B Euglenid Diplonemid SpdSyn TSA | tblastn / TSA (17 projects) | Euglenida (3035) + Diplonemea (191814) | *T. brucei* SpdSyn (Q38EH6), 298 aa | 500 / 1e-20 | 55 | ✓ |
| C | F5HFU3NK014 | C Dino Perkinsus SpdSyn | blastp / TSA proteins (tsa_nr) | Dinophyceae (2864) + Perkinsozoa (2497438) | human SpdSyn (P19623), 302 aa | 500 / 1e-20 | 100 (37 species) | ✓ Not cut off: the last hit is at 7e-22, at the cutoff. Includes the dinotoms *Kryptoperidinium*/*Durinskia* (diatom-derived copies), so apply the contamination screen. |
| D | F5HK30WR014 | D Rhizaria Haptophyta SpdSyn TSA | tblastn / TSA | Rhizaria (543769) + Haptophyta (2830) | human SpdSyn (P19623) | **100 / 0.05 (defaults)** | 100 | ✗ **Cut off at the 100 cap** (100th hit at 2e-86, almost all haptophytes). Replaced by D2 + D3. |
| D2 | F5J81YRJ014 | D2 Rhizaria SpdSyn TSA | tblastn / TSA (44 projects) | Rhizaria (543769) | human SpdSyn (P19623) | 500 / 1e-20 | 142 (41 species) | ✓ |
| D3 | F6EP33G2016 | D3 Haptophyta SpdSyn TSA | tblastn / TSA | Haptophyta (2830) | human SpdSyn (P19623) | 500 / 1e-20 | 319 (29 species) | ✓ Probably includes other aminopropyltransferases (e.g. thermospermine synthases); filter with the family tree before use. |
| E | F5HNB4Z6014 | E Katablepharid Centrohelid AdoMetDC TSA | tblastn / TSA | Kathablepharidacea (339960) + Centroplasthelida (193537) | human AdoMetDC (P17707), 334 aa | 100 / 0.05 (defaults) | 3, all centrohelids (*Raphidocystis* ×2, *Choanocystis* ×1) | ✓ The cap wasn't reached and 0.05 is looser than 1e-20, so the defaults lose nothing. **No katablepharid hit.** |
| E2 | F5J9HKT2014 | E2 Katablepharid AdoMetDC TSA | tblastn / TSA (13,410 sequences searched) | Kathablepharidacea (339960) | human AdoMetDC (P17707) | 500 / 1e-20 | 0 | ✓ (no hits) |
| E3 | F6EYV227014 | E3 Katablepharid AdoMetDC Guillardia query | tblastn / TSA (13,410 sequences) | Kathablepharidacea (339960) | *Guillardia* AdoMetDC (L1JMA4), 406 aa | 500 / 1e-5 | 0 | ✓ (no hits). An earlier attempt, F6ESCJH3014 (same but E 1e-20, no title, 0 hits), isn't kept. |
| F | F5HSK9YB016 | F Euglenid Diplonemid AdoMetDC TSA | tblastn / TSA (17 projects) | Euglenida (3035) + Diplonemea (191814) | human AdoMetDC (P17707) | 500 / 1e-20 | 35 | ✓ *E. gracilis*, *E. longa*, *Rhabdomonas costata*, *Eutreptiella* + 8 diplonemid species |
| G | F5EJKSE8014 | G Euglena AdoMetDC | tblastn / TSA (*E. gracilis*: GDJR01, GEFR01, HBDM01) | *Euglena gracilis* (3039) | *Eutreptiella* AdoMetDC (UPI0016DBF433), 375 aa | 100 / 0.05 | 10 (one gene, 2 assemblies) | ✓ Cap not reached |
| H | F5FJXSNR014 | H Euglena SpdSyn | tblastn / TSA (*E. gracilis*) | *Euglena gracilis* (3039) | *T. brucei* SpdSyn (Q38EH6) | 100 / 0.05 | 7 | ✓ Cap not reached. The hits are 1 SpdSyn, 1 thermospermine synthase-like and 1 divergent PABS gene (see `../test_results/RESULTS_euglena.md`). |

## What the searches settle

- **Katablepharids:** no AdoMetDC was detected in the available katablepharid transcriptome (E, E2, E3), using both the human query at a loose cutoff and the closer *Guillardia* query. The data is one small, fragmentary *Roombia* set, so this means "not detected", not "lost".
- **Centrohelids** do have AdoMetDC (search E, 3 transcripts). Their Glu178 status hasn't been checked yet.
- **Sampling gaps from the 9 Nov SpdSyn tree are now filled:** cryptophytes (A2), euglenids and diplonemids (B), dinoflagellates and *Perkinsus* (C), Rhizaria (D2) and haptophytes (D3).
- **New euglenid AdoMetDCs** (F): *Euglena longa* and *Rhabdomonas costata*, plus diplonemids. These can be checked for Glu178 and added to the AdoMetDC tree.
- *E. gracilis* (G, H) has already been added to both trees: see `../test_results/RESULTS_euglena.md`.
