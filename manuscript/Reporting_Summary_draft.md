# Reporting Summary (eBioMedicine), draft for author confirmation (v3.2, R35)

## Statistics

| Item | Response |
|---|---|
| **Data exclusions** | NHANES: participants missing the composite index, outcome or covariates were excluded (complete-case); the J cycle was primary for allergic outcomes (n = 7,021) and the D cycle for IgE biomarkers (n = 7,854); target-trial emulation restricted to adults >=18 with complete data (n = 4,403). HRS: 45,234 participants with deficit-accumulation index data across 16 waves. No outcome-based exclusions. |
| **Replication** | Genetic findings were replicated in independent GWAS (Han UK Biobank, TAGC, EAGLE/Budu-Aggrey, Neale lab), Pan-UK Biobank endpoints, Biobank Japan (cross-ancestry), an independent Python-3 LDSC implementation and the canonical Bulik-Sullivan implementation; pQTL findings replicated across SCALLOP (Olink), Luminex and HELIC (Olink) cohorts using the same three cis instruments (assay-level, not instrument-level, independence). Because the discovery meta-analysis is the SCALLOP Consortium study, the SCALLOP and INTERVAL replication layers share its infrastructure (non-independent) whereas HELIC, AGES-Reykjavik and the Ahola-Olli Luminex cohort are separate resources. Attempted replications that failed (EAGLE atopic dermatitis; Biobank Japan atopic dermatitis; OneK1K cell-type corroboration) are reported as boundaries, not excluded. |
| **Sample-size determination** | No statistical methods were used to predetermine sample size; all analyses used the maximal available sample within each publicly available dataset. Power statements are given where tests were underpowered (mortality, clocks, positive controls, drug-subgroup emulations). |
| **Randomization** | Not applicable (observational genetic-epidemiological study). The NHANES target-trial emulation used stabilized inverse-probability weighting, not randomization. |
| **Blinding** | Not applicable (summary-statistics-based analyses). Adjudication thresholds were fixed before discovery by file provenance (declared in Methods). |
| **Data collection** | All data were collected previously and are publicly available (FinnGen R12, Pan-UK Biobank, GWAS Catalog, Neale lab, eQTLGen, DICE, OneK1K, Tabula Sapiens 2.0, GTEx v8, pQTL consortia, NHANES, RAND HRS, FAERS, Canada Vigilance). |

## Data availability

All data are publicly available; source URLs are listed in the Data Availability statement of the manuscript. Derived tables are provided as Supplementary Tables 1–162 (with companions): the six new tables fill ST118–ST123 (ONSET01, COUP01, DRIFT01, DRIFT02, GER01, GER02) and the ATP-axis and knowledge-base tables are ST157–ST162, and are archived at Zenodo (https://doi.org/10.5281/zenodo.22152039).

## Code availability

Analysis code is available at https://github.com/wang76qiang/frailty-haematological-ageing-allergy-mr and archived at Zenodo (https://doi.org/10.5281/zenodo.22152039, concept DOI). The analysis plan is registered on the Open Science Framework (https://osf.io/f4v79; registrations qm2hy and f9d3w, post hoc relative to the R5/R6 analyses; internally pre-declared families for post hoc layers are documented in the archived analysis package).

## Ethics

The study used publicly available, de-identified data only; no new individual-level data were collected and no additional institutional review board approval was required (see Ethics declaration in the manuscript).

---
*Author action: confirm or edit each field in the eBioMedicine submission system; this draft mirrors manuscript content as of v3.2/R35 (2026-10-09).*
