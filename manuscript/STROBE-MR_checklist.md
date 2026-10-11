# STROBE-MR checklist

**Manuscript:** One law, two dimensions of ageing: a drifting haematological type-2 set-point couples functional ageing to allergic disease, while molecular clocks remain largely decoupled from the primary allergic endpoints (v3.2, R38)

This checklist maps each reporting requirement to its location in the manuscript and to the supporting additional files (Supplementary Tables, ST; Supplementary Figures) and to the eight main composite figures of the current version. Figure/panel pointers follow the current scheme: Figure 1–8 (further detailed in `02_主图_figures/lancet/MAIN_FIGURES_blueprint.md`) and Supplementary Figures 1–16 (renumbered to first-citation order). The manuscript follows the STROBE-MR statement (Skrivankova et al., JAMA 2021).

| **Item** | **Location** |
| --- | --- |
| Pre-specified decision tree and thresholds | Methods 2.1; Fig. 1a |
| Instrument selection and F-statistics | Methods 2.2; Table 2; ST2 |
| Data sources (FinnGen, Pan-UKBB, Astle, eQTLGen, DICE, pQTL consortia, NHANES, FAERS/CV) | Methods; Data availability |
| Primary estimator (IVW) and sensitivity estimators (Egger, weighted median, PRESSO, RAPS) | Methods 2.2; Fig. 1a; ST3-ST5 |
| Pleiotropy assessment (Egger intercepts; multivariable MR adjusted for BMI/smoking/education; joint frailty + index MVMR; immune-locus-excluded instruments) | Methods 2.2; Fig. 1d; ST61-ST65 |
| Reverse MR | Methods 2.2; ST6-ST7 |
| Winner's-curse and selection sensitivity | Methods 2.2; ST8-ST9 |
| Positive/negative controls (pQTL) | Fig. 7c; ST42, ST45 |
| Colocalization adjudication | Fig. 5a,b; ST22, ST29-ST31 |
| Independent replication and meta-analysis | Methods 2.2; Fig. 1b; ST57-ST58 |
| LD score regression (h2, rg; MHC-free verification, instrument consistency and simulation calibration) | Methods 2.3; Fig. 1e, Fig. 1g; ST59-ST60, ST76-ST79 |
| Multi-signal colocalization and SMR/HEIDI at 5q31 | Methods 2.6; Fig. 5b, Fig. 5h, Fig. 5i; ST66-ST69 |
| NHANES index validity | Methods 2.4; Fig. 3i, Fig. 3j; ST70-ST73 |
| On-target safety scan | Methods 2.7; Fig. 7d; ST74-ST75 |
| Phenome-wide MR specificity scan and sensitivity battery | Methods 2.2; Fig. 1k; ST80-ST82; Supplementary Fig. 2 |
| HRS longitudinal validation (incident disease, mortality, early-warning signals) | Methods 2.5; Fig. 4h, Fig. 4i; ST83-ST86; Supplementary Fig. 7 |
| OneK1K remediation and independent cell-type MR | Methods 2.6; Fig. 5e; ST87-ST88; Supplementary Fig. 8 |
| Blood-cell-trait decomposition MR | Methods 2.2; Fig. 6l; ST89; Supplementary Fig. 5 |
| Single-cell ageing atlas, index anchoring and target cell localization | Methods 2.6; Fig. 6a, Fig. 6e, Fig. 6f; ST90-ST93; Supplementary Figs. 10, 13 |
| Disease-signature transfer (descriptive) | Methods 2.6; Supplementary Fig. 9 |
| Onset-stratified asthma MR and LDSC | Methods 2.2; Fig. 1h; ST94-ST98; Supplementary Fig. 3 |
| Molecular mediation MR (IL18/TNFSF14 downstream scan) | Methods 2.7; Fig. 7e; ST99; Supplementary Fig. 15 |
| Tissue-context eQTL analysis (GTEx) | Methods 2.7; Fig. 7k; ST100; Supplementary Fig. 16 |
| Human gene-burden lookup (natural knockout) | Methods 2.7; Fig. 7f; ST101 |
| Coloc prior sensitivity (IL4 locus and proteome) | Methods 2.6; Fig. 5k; ST102 |
| LDSC external benchmark versus official LDSR | Methods 2.3; ST103 |
| Canonical Bulik-Sullivan LDSC re-computation (munged summary statistics, published Phase 3 HapMap3 weights, MHC-free; MHC-included sensitivity) and canonical frailty-versus-clock contrasts | Methods 2.3; Fig. 1g; ST137-ST138 |
| Donor-level single-cell atlas expression-age analysis (regulatory-T-cell five-test family, translational four-test family, ambient-RNA falsification, lung localization) | Methods 2.6; Fig. 6a, Fig. 6g; ST139-ST145; Supplementary Figs. 11, 12, 14 |
| Target trial emulation (stabilized IPTW, survey-weighted, PSU-clustered, E-values; drug-subgroup emulations descriptive) | Methods 2.4; Fig. 4g, Fig. 4l; ST146, ST152; Supplementary Fig. 6 |
| Formal set-point drift model (OU homeostat; drift, variance-resilience, drift shape, gain invariance) | Methods 2.4; Results; Fig. 3a, Fig. 3c, Fig. 3e; ST153-ST156, ST154b; Supplementary Fig. 4 |
| Gene-level target triangulation (MR, colocalization, PRISM, ChEMBL, FAERS) and LINCS drug-connectivity screening (ageing and alarmin signatures, two 11-drug families) | Methods 2.8; Fig. 8c, Fig. 8d, Fig. 8e; ST147-ST151 |
| HRS companion analyses (component-exclusion Cox; mean-matched early-warning) | Methods 2.5; ST104-ST105 |
| Cross-cohort/cross-platform pQTL replication (SCALLOP, AGES, INTERVAL SomaScan, Luminex, HELIC) | Methods 2.7; Fig. 7a, Fig. 7b; ST106, ST127, ST128 |
| Incident-allergy cohort pipeline (synthetic self-validation; awaits UKB/AoU extracts) | Methods; Discussion (Limitations); ST112 |
| TNFSF14 multi-platform and leave-one-out sensitivity | Discussion; Fig. 8l; ST113 |
| IL4 independent single-cell eQTL checks (OneK1K; independent-resource search) | Results 3.8; Discussion; ST114, ST115 |
| Full Tabula Sapiens 2.0 Immune atlas re-analysis (donor-level, complexity-adjusted) | Results 3.9; ST116, ST116b |
| Epigenetic-clock validity check (GSE40279 methylation) | SI table legend; ST117 |
| Index weight sensitivity (equal-weight and eosinophil-free indices) | Results 3.4, 3.5; Discussion (Limitations); Fig. 4f; ST124 |
| Independent ageing exposures (telomere length, longevity negative controls) | Results 3.1; ST125 |
| Independent Python-3 LD score regression replication | Results 3.1; Fig. 1l; ST126 |
| gnomAD loss-of-function constraint and tractability (proxy for gene-burden lookup) | Supplementary Note 9; ST129 |
| Non-HLA (immune-locus-excluded) frailty re-adjudication | Results 3.1; Fig. 1c; ST131 |
| Type-2 disease-spectrum MR (Pan-UK Biobank asthma, rhinitis, urticaria; non-type-2 negative controls) | Results 3.3; Fig. 2h; ST130 |
| External healthspan dimension re-adjudication (Timmers 2019, MHC-excluded) | Results 3.1; ST132 |
| Genetic-level eosinophil decomposition of the composite index (eosinophil-free / eosinophil-only instruments) | Results 3.5; Discussion (Limitations); ST133 |
| Formal frailty-versus-clock genetic-correlation contrasts (12 contrasts, two LDSC implementations) | Results 3.1; Fig. 1f; ST134 |
| Steiger SNP-level directionality audit (forward MR; per-SNP companion ST135b) | Methods 2.2; Results 3.1; Fig. 1i; ST135, ST135b |
| Eosinophil decomposition: independent-cohort replication in Pan-UK Biobank | Results 3.5; Discussion (Limitations); ST136 |
| East Asian (Biobank Japan) cross-ancestry replication | Methods 2.2; Results 3.3; Fig. 2i; ST107-ST111 |
| Limitations (AD replication underpowered; compartment mismatch; multi-signal coloc approximations; underpowered mortality; pharmacovigilance confounding) | Discussion |
| Knowledge-layer enrichment (17-gene type-2/alarmin panel membership and Fisher enrichment in the 500 curated human Aging Atlas genes; IL18 absence) | Methods 2.8; Results 3.12; Fig. 8f; ST162 (KB01) |
| Full-instrument Steiger directionality audit including the magnitude-carrying composite index (37% correct-direction for asthma; 19 wrong vs 3 correct among significant SNPs) | Methods 2.2; Results 3.1; Limitation (vii); ST135, ST135b |
| Objective-IgE vs eosinophil component decomposition (eosinophil-free index null for IgE; eosinophil-only index carries it) | Results 3.4, 3.5; Fig. 4e; ST124; Supplementary Note 11 |
| Incident-allergy cohort pipeline is SYNTHETIC only (simulated validation; no empirical cohort) | Methods; Discussion (Limitations); ST112 |
