# Frailty and haematological type-2 genetic architecture in allergic disease

**Reproducibility package** for the manuscript submitted to *Genome Medicine*:

> Dimension-resolved Mendelian randomization implicates a functional frailty
> dimension and a rhinitis-anchored blood-cell composite in allergic disease
> without a consistent epigenetic-clock signal

Author: Qiang Wang, Chinese PLA Center for Disease Control and Prevention, Beijing, China.

## Study summary

We compared five ageing instruments (four epigenetic clocks and a frailty index)
plus a composite blood-cell-trait index against FinnGen R12 asthma, allergic
rhinitis and atopic dermatitis by Mendelian randomization under an a priori (by
file provenance, unregistered) adjudication tree, with independent replication,
MHC-free LD score regression, multivariable MR, NHANES and HRS validation,
cell-type and cis-pQTL MR, colocalization and pharmacovigilance. The epigenetic
clocks show no consistent, concordant signal; a functional (frailty) dimension
and, most replicably for allergic rhinitis, an eosinophil-weighted blood-cell
composite are implicated. The association is not HLA-driven and has not been
shown to be causal.

## Contents

| Path | Description |
|------|-------------|
| `scripts/` | Analysis scripts, including the R5 core (orthogonal ageing MR, reverse MR, LDSC, colocalization, NHANES, proteome MR, pharmacovigilance) and the editorial-revision scripts `r16_01`-`r16_04` |
| `tables/` | Derived result tables ST01-ST117 (ST116b companion) and ST124-ST136 (ST135b companion) as machine-readable CSV/TSV, plus `supplementary_data_index.csv` |
| `README.md` | This file |
| `LICENSE` | MIT license |

## Key editorial-revision scripts

- `r16_01_eosinophil_free_mr.py` - genetic-level eosinophil decomposition (ST133)
- `r16_01b_harmonisation_diagnostic.py` - palindromic-handling diagnostic
- `r16_02_rg_contrast.py` - formal frailty-versus-clock rg contrasts (ST134)
- `r16_03_steiger_directionality.py` - per-SNP Steiger directionality audit (ST135/ST135b)
- `r16_04_eosinophil_free_panukbb.py` - independent-cohort replication of the decomposition (ST136)

## Reproducing

Requires Python 3.13. Data inputs are publicly available GWAS/eQTL summary
statistics and NHANES public-use files (accessions listed in the manuscript's
Data availability section); derived tables in `tables/` are version-frozen with
this deposit.

## License

MIT (see `LICENSE`).
