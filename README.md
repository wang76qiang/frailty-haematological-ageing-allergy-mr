# Frailty and allergic disease share genetic architecture, partly through overlapping item content

Reproducibility package for the **Clinical & Experimental Allergy** submission (Route B, v8.6; 21 October 2026).

## Current version

The current submission archive is `CEA_RouteB8.6_submission_package.zip` in this repository. The Zenodo release is:

- Version DOI: https://doi.org/10.5281/zenodo.23533619
- Concept DOI: https://doi.org/10.5281/zenodo.22152039

The manuscript asks whether part of the genetic correlation between frailty and allergic disease is definitional, because the 49-item UK Biobank frailty index includes asthma (item 28) and hay fever/allergic rhinitis/eczema (item 33) as two of its 49 items.

## Key reproducible results

- FI item audit: 2/49 outcome items; 9/49 pain items.
- FinnGen frailty–asthma rg: 0.467 (SE 0.039); rhinitis 0.361; atopic dermatitis 0.246.
- European-ancestry TAGC replication: random-effects rg 0.378 (SE 0.048; P = 1.9 × 10−15); fixed-effects 0.362 (SE 0.046); intercept −0.003 (SE 0.006); N = 127,669.
- Non-HLA MR: FinnGen asthma 3.85 (2.69–5.52); FinnGen rhinitis 2.67 (1.92–3.71); Pan-UKB asthma 2.38 (1.86–3.04).
- NHANES blood-count index: OR 1.31 (1.21–1.41) per SD; equal-weight 1.31 (1.20–1.42); eosinophil-free 1.18 (1.08–1.28); eosinophil-only 1.24 (1.16–1.33).
- Log-IgE component betas: full index 0.442; without eosinophils −0.127; eosinophils only 0.301.

## Package contents

- `manuscript.md`, `manuscript.docx`, `manuscript.pdf`
- `cover_letter.md`, `cover_letter.docx`
- five main figures (`.png`, `.pdf`, `.tiff`) and a 300-dpi graphical abstract
- `Table_1-3.md`, `Table_1-3.docx`
- `Supplementary_information.md`, `Supplementary_information.docx`
- current `STROBE-MR_checklist.docx`
- `source_data/` for Figures 2f, 4 and 5a
- `README_submission_manifest.md`

## Historical files

The `manuscript/`, `figures/`, `tables/` and `supplementary/` directories contain the earlier eBioMedicine-oriented reproducibility package (v3.x). They are retained for provenance. The v8.6 CEA package is the zip at the repository root and the Zenodo version linked above.

## Analysis plan

OSF: https://osf.io/f4v79
