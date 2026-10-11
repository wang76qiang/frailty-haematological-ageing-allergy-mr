# Frailty / haematological type-2 set-point and allergic disease — reproducibility package

Reproducibility package for the manuscript submitted to **eBioMedicine** (v3.2, R35/R36):

> **One law, two dimensions of ageing: a drifting haematological type-2 set-point
> couples functional ageing to allergic disease, while molecular clocks remain
> largely decoupled from the primary allergic endpoints**

Author: Qiang Wang, Chinese PLA Center for Disease Control and Prevention, Beijing, China.

## Study summary
We compared five ageing instruments (four epigenetic clocks and a frailty index) plus a
composite haematological type-2 index against FinnGen R12 asthma, allergic rhinitis and
atopic dermatitis by Mendelian randomization under a pre-specified adjudication tree, with
independent replication, MHC-free LD score regression (three implementations including the
canonical Bulik-Sullivan tool), multivariable and reverse MR, an 85-test outcome-wide
coupling map, NHANES and HRS validation, a formal set-point drift model, cell-type and
cis-pQTL MR, colocalization, single-cell atlases, and pharmacological/drug-target layers.
A functional (frailty) dimension — not the epigenetic clocks — carries the genetic overlap
with allergic disease; the associated blood-cell axis behaves as an age-invariant type-2
set-point whose level drifts with age; and IL18 is nominated as a platform-robust,
barrier-resident translational node. The association is not shown to be causal.

## Contents
- `tables/` — all derived result tables ST01–ST162 (with companions), `supplementary_data_index.csv`,
  the machine-readable `supplementary_tables_csv.zip`, and the re-laid-out compendium
  (`Supplementary_Tables_1-162.docx` / `.pdf`). The six formerly unnumbered supplementary
  tables fill ST118–ST123 (ONSET01, COUP01, DRIFT01, DRIFT02, GER01, GER02); ATP-axis tables
  and the knowledge-base table are ST157–ST162.
- `manuscript/` — `manuscript_main.md` / `.docx` (v3.2), `main_tables.*`, cover letter,
  Reporting Summary and STROBE-MR checklist.
- `supplementary/` — `supplementary_information.md` / `.docx` and the supplementary figures
  (including `Supplementary_Fig_16_participant_flow`).
- `scripts/` — analysis scripts (R5 core plus the R9–R24 revision layers).

## Reproducing
Requires Python 3.13+. Data inputs are publicly available GWAS/eQTL/pQTL summary statistics
and NHANES public-use files (accessions listed in the manuscript's Data availability section);
derived tables in `tables/` are version-frozen with this deposit.

## Note on versions
This repository previously carried an earlier ("dimension-resolved MR", Genome Medicine-oriented)
package; that content is retained in the git history. The current `main` branch is the
eBioMedicine v3.2 package. Zenodo concept DOI: 10.5281/zenodo.22152039.

License: MIT (see `LICENSE`).
