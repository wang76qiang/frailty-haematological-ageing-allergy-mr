# Frailty-associated haematological ageing, not epigenetic age, shares causal architecture with allergic susceptibility

**Reproducibility package** for the manuscript submitted to *Nature Aging*.

Author: Qiang Wang, Chinese PLA Center for Disease Control and Prevention, Beijing, China.

## Study summary

This study adjudicates which dimension of biological ageing causally shapes allergic
susceptibility, using a pre-registered decision tree over five orthogonal genetic
instruments of ageing (four epigenetic clocks and a frailty index) tested against
FinnGen GWAS for asthma, allergic rhinitis and atopic dermatitis by Mendelian
randomization (MR), with population validation in 14,878 NHANES participants,
cell-type-resolved MR, methylome convergence, plasma proteome MR with colocalization,
and cross-database pharmacovigilance. The evidence adjudicates a frailty-specific
causal dimension (path C): a functional, haematologically visible ageing axis, not
epigenetic or chronological age, shares causal architecture with allergic disease.

## Contents

| Path | Description |
|------|-------------|
| `scripts/` | R5 analysis scripts (32 Python + 5 shell), including the five main-figure renderers `fig1_path_c.py`–`fig5_translation.py` and `nature_style.py` |
| `tables/` | All 56 derived result tables (ST01–ST56) as machine-readable CSV plus `supplementary_data_index.csv` |
| `README.md` | This file |
| `CITATION.cff` | Citation metadata for this deposit |
| `LICENSE` | MIT license |
| `requirements.txt` | Python dependencies |

## Analysis layers

- **MR of ageing instruments**: `01_orthogonal_aging_mr.py`, `02_bidirectional_mr.py`,
  `03_heterogeneity_governance.py`, `03_heterogeneity_presso_raps.py`, `04_winners_curse.py`,
  `04_astle_sensitivity.py`
- **5q31/IL4–IL13 locus**: `06_il4_locus_coloc.py`, `06_il4_locus_adjudication.py`
- **NHANES population validation**: `07_nhanes_survey_weighted.py`, `08_age_interaction_pooled.py`
- **Plasma proteome MR**: `09_pqtl_drug_target_mr.py`
- **Epigenome-wide convergence**: `10_ewas_convergence.py`
- **Cell-type-resolved MR**: `12_celltype_mr_matrix.py`, `12_onek1k_replication_*.py`
- **Pharmacovigilance**: `13_pharmacovigilance_three_source.py`, `13_three_source_vigilance.py`
- **Figure rendering**: `fig1_path_c.py`–`fig5_translation.py`, `nature_style.py`

## Reproducing

Requires Python 3.13 with the packages in `requirements.txt`. Data inputs are publicly
available GWAS/eQTL summary statistics and NHANES public-use files (accessions listed in
the manuscript's Data availability section); derived tables in `tables/` are version-frozen
with this submission.

## Citation

Please cite the manuscript and this Zenodo deposit; see `CITATION.cff`.

## License

MIT (see `LICENSE`).
