# Canonical LDSC re-computation notes (R18)

Official bulik/ldsc v1.0.1 run under Python 3.13 with mechanical interpreter-compatibility
patches only (see scripts/r17_00_patch_py3.py and scripts/r17_01_canonical_ldsc.py).
Protocol: munge_sumstats.py on published GWAS (McCartney clocks GCST90014288-92; Atkins
frailty GCST90020053; FinnGen R12 outcomes; w_hm3.snplist), then ldsc.py --rg on published
1000 Genomes Phase 3 HapMap3 weights (MHC excluded; .l2.M derived as reference SNP counts;
--not-M-5-50). Outputs: ST137 (rg, MHC-free + MHC-included sensitivity) and ST138
(formal frailty-versus-clock contrasts, both SE conventions). Raw GWAS downloads and the
patched ldsc tree are re-creatable from the scripts; per-run logs are archived under
docs/canonical_ldsc_logs/ in this commit.

Headline: frailty rg reproduced (asthma 0.467, rhinitis 0.361, atopic dermatitis 0.246);
median absolute difference from in-house engine 0.025.
