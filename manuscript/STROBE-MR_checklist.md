# STROBE-MR checklist

**Manuscript:** One law, two dimensions of ageing: a drifting haematological type-2 set-point couples functional ageing to allergic disease, while molecular clocks remain largely decoupled from the primary allergic endpoints (v3.2, R35/R36)

This checklist maps each reporting requirement to its location in the manuscript and to the supporting Additional files (Supplementary Tables, ST) and Supplementary Figures, and to the eight main composite figures of the current version (Figure 1–8; the authoritative panel map is given in `02_主图_figures/lancet/MAIN_FIGURES_blueprint.md`). Earlier drafts used a different figure/panel scheme, so a small number of legacy "Fig." panel pointers in this checklist are retained and should be read against the current Figure 1–8 panels. The manuscript follows the STROBE-MR statement (Skrivankova et al., JAMA 2021).

| **Item** | **Location** |
| --- | --- |
| Pre-specified decision tree and thresholds | Methods; Fig. 1a |
| Instrument selection and F-statistics | Methods; Figs. 1b, 2d |
| Data sources (FinnGen, Pan-UKBB, Astle, eQTLGen, DICE, pQTL consortia, NHANES, FAERS/CV) | Methods; Data availability |
| Primary estimator (IVW) and sensitivity estimators (Egger, weighted median, PRESSO, RAPS) | Methods; Figs. 1c, 2e |
| Pleiotropy assessment (Egger intercepts; multivariable MR adjusted for BMI/smoking/education, ST61; joint frailty + index MVMR, ST64-ST65; immune-locus-excluded instruments, ST62-ST63) | Fig. 1g; Methods |
| Reverse MR | Fig. 2e,f |
| Winner's-curse and selection sensitivity | Fig. 2g,h |
| Positive/negative controls (pQTL) | Fig. 6d,e |
| Colocalization adjudication | Figs. 5b,c |
| Independent replication and meta-analysis | Methods; ST57-ST58 |
| LD score regression (h2, rg; MHC-free verification, instrument consistency and simulation calibration) | Methods; ST59-ST60, ST76-ST79 |
| Multi-signal colocalization and SMR/HEIDI at 5q31 | Methods; ST66-ST69 |
| NHANES index validity | Methods; ST70-ST73 |
| On-target safety scan | Methods; ST74-ST75 |
| Phenome-wide MR specificity scan and sensitivity battery | Methods; ST80-ST82; Supplementary Fig. 1 |
| HRS longitudinal validation (incident disease, mortality, early-warning signals) | Methods; ST83-ST86; Supplementary Fig. 2 |
| OneK1K remediation and independent cell-type MR | Methods; ST87-ST88; Supplementary Fig. 3 |
| Blood-cell-trait decomposition MR | Methods; ST89; Supplementary Fig. 4 |
| Single-cell ageing atlas, index anchoring and target cell localization | Methods; ST90-ST93; Supplementary Figs. 5-6 |
| Disease-signature transfer (descriptive) | Methods; Supplementary Fig. 7 |
| Onset-stratified asthma MR and LDSC | Methods; ST94-ST98; Supplementary Fig. 8 |
| Molecular mediation MR (IL18/TNFSF14 downstream scan) | Methods; ST99; Supplementary Fig. 9 |
| Tissue-context eQTL analysis (GTEx) | Methods; ST100; Supplementary Fig. 10 |
| Human gene-burden lookup (natural knockout) | Methods; ST101 |
| Coloc prior sensitivity (IL4 locus and proteome) | Methods; ST102 |
| LDSC external benchmark versus official LDSR | Methods; ST103 |
| Canonical Bulik-Sullivan LDSC re-computation (munged summary statistics, published Phase 3 HapMap3 weights, MHC-free; MHC-included sensitivity) and canonical frailty-versus-clock contrasts | Methods; ST137-ST138 |
| Donor-level single-cell atlas expression-age analysis (regulatory-T-cell five-test family, translational four-test family, ambient-RNA falsification, lung localization) | Methods; ST139-ST145; Supplementary Figs. 11-13 |
| Target trial emulation (stabilized IPTW, survey-weighted, PSU-clustered, E-values; drug-subgroup emulations descriptive) | Methods; ST146, ST152; Supplementary Fig. 14 |
| Formal set-point drift model (OU homeostat; drift, variance-resilience, drift shape, gain invariance) | Methods (Set-point drift model); Results; ST153-ST156, ST154b; Supplementary Fig. 15 |
| Gene-level target triangulation (MR, colocalization, PRISM, ChEMBL, FAERS) and LINCS drug-connectivity screening (ageing and alarmin signatures, two 11-drug families) | Methods; ST147-ST151 |
| HRS companion analyses (component-exclusion Cox; mean-matched early-warning) | Methods; ST104-ST105 |
| Cross-cohort/cross-platform pQTL replication (SCALLOP, AGES, INTERVAL SomaScan, Luminex, HELIC) | Methods; ST106, ST127, ST128 |
| Incident-allergy cohort pipeline (synthetic self-validation; awaits UKB/AoU extracts) | Methods; Discussion (Limitations); ST112 |
| TNFSF14 multi-platform and leave-one-out sensitivity | Discussion; ST113 |
| IL4 independent single-cell eQTL checks (OneK1K; independent-resource search) | Results; Discussion; ST114, ST115 |
| Full Tabula Sapiens 2.0 Immune atlas re-analysis (donor-level, complexity-adjusted) | Results; ST116, ST116b |
| Epigenetic-clock validity check (GSE40279 methylation) | SI table legend; ST117 |
| Index weight sensitivity (equal-weight and eosinophil-free indices) | Results; Discussion (Limitations); ST124 |
| Independent ageing exposures (telomere length, longevity negative controls) | Results; ST125 |
| Independent Python-3 LD score regression replication | Results; Fig. 3b; ST126 |
| gnomAD loss-of-function constraint and tractability (proxy for gene-burden lookup) | Supplementary Note 9; ST129 |
| Non-HLA (immune-locus-excluded) frailty re-adjudication | Results; ST131 |
| Type-2 disease-spectrum MR (Pan-UK Biobank asthma, rhinitis, urticaria; non-type-2 negative controls) | Results; ST130 |
| External healthspan dimension re-adjudication (Timmers 2019, MHC-excluded) | Results; ST132 |
| Genetic-level eosinophil decomposition of the composite index (eosinophil-free / eosinophil-only instruments) | Results; Discussion (Limitations); ST133 |
| Formal frailty-versus-clock genetic-correlation contrasts (12 contrasts, two LDSC implementations) | Results; ST134 |
| Steiger SNP-level directionality audit (forward MR; per-SNP companion ST135b) | Methods; Results; ST135, ST135b |
| Eosinophil decomposition: independent-cohort replication in Pan-UK Biobank | Results; Discussion (Limitations); ST136 |
| East Asian (Biobank Japan) cross-ancestry replication | Methods; ST107-ST111 |
| Limitations (AD replication underpowered; compartment mismatch; multi-signal coloc approximations; underpowered mortality; pharmacovigilance confounding) | Discussion |
| Knowledge-layer enrichment (17-gene type-2/alarmin panel membership and Fisher enrichment in the 500 curated human Aging Atlas genes; IL18 absence) | Methods 2.8; Results 3.12; KB01 |
| Full-instrument Steiger directionality audit including the magnitude-carrying composite index (37% correct-direction for asthma; 19 wrong vs 3 correct among significant SNPs) | Methods 2.2; Results 3.1; Limitation (vii); ST135, ST135b |
| Objective-IgE vs eosinophil component decomposition (eosinophil-free index null for IgE; eosinophil-only index carries it) | Results 3.4, 3.5; ST124; Supplementary Note 11 |
| Incident-allergy cohort pipeline is SYNTHETIC only (simulated validation; no empirical cohort) | Methods; Discussion (Limitations); ST112 |
