# Supplementary Figures 1-10, provenance

These image files were regenerated from the deposited result tables
(`04_补充数据_supplementary_data/ST*.csv`) by
`src/r15/10_make_supplementary_figures.py`.

| Figure | Content | Source table(s) | File |
|---|---|---|---|
| Supplementary Fig. 1 | Phenome-wide MR specificity scan | ST80 | `Supplementary_Fig_1_phenomewide.png` |
| Supplementary Fig. 2 | HRS longitudinal validation and early-warning signals | ST83, ST86 | `Supplementary_Fig_2_hrs.png` |
| Supplementary Fig. 3 | OneK1K independent cell-type MR (FDR-significant pairs) | ST88 | `Supplementary_Fig_3_onek1k.png` |
| Supplementary Fig. 4 | Blood-cell-trait decomposition MR | ST89 | `Supplementary_Fig_4_bloodtrait.png` |
| Supplementary Fig. 5 | Single-cell ageing atlas (FDR-significant age slopes) | ST90 | `Supplementary_Fig_5_sc_atlas.png` |
| Supplementary Fig. 6 | IL18 / TNFSF14 cell of origin | ST93 | `Supplementary_Fig_6_target_cells.png` |
| Supplementary Fig. 7 | Disease-signature transfer (GSE193816, descriptive) | GSE193816 airway allergen-challenge h5ad | `Supplementary_Fig_7_signature_transfer.png` |
| Supplementary Fig. 8 | Onset-stratified asthma MR | ST94 | `Supplementary_Fig_8_onset.png` |
| Supplementary Fig. 9 | Molecular mediation MR | ST99 | `Supplementary_Fig_9_mediation.png` |
| Supplementary Fig. 10 | GTEx tissue-context eQTL | ST100 | `Supplementary_Fig_10_gtex.png` |

Supplementary Fig. 7 was regenerated from the allergen-challenge airway dataset
(GSE193816; `src/r7/03_disease_signature_transfer.py`); the reproduced cell-level
correlations match the SI (ageing-vs-Th2 r = 0.178; ageing-vs-alarmin r = 0.007).
It remains descriptive (cell-level pseudoreplication), per Supplementary Note 5(vii).
The single-cell atlas figures (5, 6) derive from Tabula Sapiens 2.0 h5ad files
(Blood, Lung, large/small intestine, Immune).
