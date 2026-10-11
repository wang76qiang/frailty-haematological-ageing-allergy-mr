# Main tables

## Table 1 | Study populations and data resources

| Resource | n | Ancestry | Role | Source |
|---|---|---|---|---|
| NHANES cycle D (2005-06) | 7,854 (IgE) | US, admixed | population validation (objective IgE) | NHANES |
| NHANES cycle J (2017-18) | 7,021 (1,092 events) | US, admixed | population validation (allergic outcomes) | NHANES |
| NHANES D+J (pooled) | 14,878 | US, admixed | set-point model / target-trial emulation | NHANES |
| HRS (RAND, 1992-2022) | 45,234 (723,744 person-waves) | US | longitudinal validation of frailty construct | RAND HRS |
| FinnGen R12 | - | Finnish (EUR) | discovery GWAS (3 outcomes) | FinnGen |
| Han/Y UKB asthma | - | EUR | replication | GCST010042 |
| TAGC asthma | - | EUR | replication | GCST005212 |
| Neale UKB rhinitis/asthma | - | EUR | replication | Neale lab |
| Budu-Aggrey atopic dermatitis | - | EUR | replication | GCST90244787 |
| Pan-UK Biobank | - | EUR | independent replication | Pan-UKB |
| Biobank Japan | - | East Asian | cross-ancestry replication | BBJ |
| eQTLGen | - | EUR | whole-blood cis-eQTL | eQTLGen |
| DICE | - | EUR | immune cell-type cis-eQTL | DICE |
| OneK1K | 982 donors | EUR | single-cell cis-eQTL (independent) | OneK1K |
| Tabula Sapiens 2.0 | 592,317 immune; 65,847 lung cells | - | single-cell atlas | TS 2.0 |
| pQTL platforms | - | EUR | proteome MR (SCALLOP/Luminex/HELIC/AGES/INTERVAL) | Zhao 2023 |
| FAERS / Canada Vigilance | - | - | pharmacovigilance | FAERS; CV |

*Note: the pooled NHANES analytic set (n = 14,878; as used in the pooled age-interaction model, ST16/ST17/ST17b) is not the arithmetic sum of the cycle-specific outcome subsets (D 7,854 for IgE; J 7,021 for allergic outcomes, which sum to 14,875); the model-specific complete-case sets are smaller (set-point drift model n = 9,360; no-CRP sensitivity n = 9,510; target-trial emulation n = 4,403). Sources: ST14, ST16, ST17, ST146, ST153.*

## Table 2 | Dimension adjudication: MR effect and genetic correlation

| Instrument | Outcome | MR OR (95% CI) | rg (SE) |
|---|---|---|---|
| Frailty | ALLERG_ASTHMA | 3.44 (2.91-4.06) | 0.458 (0.038) |
| Frailty | ALLERG_RHINITIS | 1.69 (1.46-1.96) | 0.353 (0.046) |
| Frailty | L12_ATOPIC | 1.46 (1.31-1.62) | 0.250 (0.033) |
| GrimAge | ALLERG_ASTHMA | 1.12 (1.05-1.19) | 0.165 (0.070) |
| GrimAge | ALLERG_RHINITIS | 1.00 (0.94-1.06) | -0.004 (0.085) |
| GrimAge | L12_ATOPIC | 0.97 (0.93-1.02) | 0.091 (0.073) |
| Hannum | ALLERG_ASTHMA | 1.02 (1.00-1.03) | 0.113 (0.065) |
| Hannum | ALLERG_RHINITIS | 0.97 (0.96-0.99) | 0.036 (0.084) |
| Hannum | L12_ATOPIC | 1.01 (1.00-1.02) | 0.116 (0.071) |
| IEAA | ALLERG_ASTHMA | 0.99 (0.98-0.99) | 0.083 (0.060) |
| IEAA | ALLERG_RHINITIS | 1.01 (1.00-1.01) | -0.032 (0.073) |
| IEAA | L12_ATOPIC | 1.00 (0.99-1.00) | 0.161 (0.066) |
| PhenoAge | ALLERG_ASTHMA | 0.94 (0.94-0.95) | 0.082 (0.078) |
| PhenoAge | ALLERG_RHINITIS | 0.99 (0.98-1.00) | -0.045 (0.093) |
| PhenoAge | L12_ATOPIC | 0.99 (0.99-1.00) | 0.274 (0.102) |

## Table 3 | Set-point drift model parameters

| Cycle | Parameter | Estimate (95% CI) |
|---|---|---|
| J | drift_beta_per_decade | 0.089 (0.070-0.108) |
| J | S0_at_age50 | 0.052 (-0.015-0.119) |
| D | drift_beta_per_decade | 0.056 (0.038-0.074) |
| D | S0_at_age50 | 0.058 (0.004-0.113) |
| pooled | drift_beta_per_decade | 0.074 (0.060-0.088) |
| pooled | S0_at_age50 | 0.057 (0.013-0.100) |
| J | hinge best knot (y) | 60 |
| D | hinge best knot (y) | 55 |
| pooled | hinge best knot (y) | 55 |
| pooled | pre_knot_slope_per_decade | 0.039 (0.014-0.064) |
| pooled | post_knot_increment_per_decade | 0.098 (0.041-0.155) |
| pooled | post_knot_slope_per_decade | 0.137 (0.098-0.176) |

## Table 4 | Translational targets (IL18 / TNFSF14)

| Target | pQTL coloc PP.H4 | Platforms replicated | Pooled OR (95% CI) | Safety-scan direction | Gene-burden | Druggability |
|---|---|---|---|---|---|---|
| IL18 | 0.968 | SCALLOP, Luminex, HELIC (3) | 1.19 (1.15-1.23) | concordant (no adverse signal) | null | tadekinig alfa (IL-18BP) |
| TNFSF14 | 0.976 | AGES, INTERVAL SomaScan (2) | 1.13 (1.05-1.22) | inverse septicaemia/fracture | null | caution (single-instrument) |
| TNF / TRAIL | 0.010 / 0.282 | - | - | - | - | rejected (linkage) |
| IL33 (comparator) | - | - | - | protective | protective across 3 diseases | alarmin counterpoint |

*Note: the colocalization and cross-platform replication evidence for IL18 and TNFSF14 is anchored to atopic dermatitis (European ancestry), the only outcome arm with a cis-pQTL signal passing FDR < 0.1; replication uses the same three cis instruments across platforms (assay-level, not instrument-level, independence). The discovery meta-analysis is the SCALLOP Consortium study, so the SCALLOP and INTERVAL replication layers share its infrastructure and are non-independent, whereas HELIC, AGES-Reykjavik and the Ahola-Olli Luminex cohort are separate resources (see Supplementary Information).*

Source files: ST01, ST60, ST153, ST155, ST156, ST45, ST106, ST128, ST74, ST101.
