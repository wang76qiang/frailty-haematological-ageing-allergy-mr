#!/usr/bin/env python3
"""
R5-03 Heterogeneity governance for the index -> allergy MR.

(a) Rebuild the harmonised 93-SNP index instruments vs the 3 FinnGen allergy
    outcomes and run the iterative PRESSO-style outlier screen
    (r1_utils.mr_presso_outliers, alpha = 0.05, Bonferroni over the SNP set).
    -> results/r5/tables/r5_presso_outliers.csv

(b) MR-RAPS: Huber-robust IVW (IRLS, k = 1.345, weights 1/sy^2,
    overdispersion phi = max(1, robust residual chi2 / df)) benchmarked against
    IVW fixed/random, MR-Egger and weighted median.
    -> results/r5/tables/r5_raps_vs_ivw.csv
"""

import os
import sys
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r5x_utils as C

OUTCOMES = list(C.FINNGEN_OUTCOMES)
N_EXP = 400000
ALPHA = 0.05

log_lines = []


def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)


def load_index_instruments():
    path = os.path.join(C.ROOT, "results", "r4", "tables", "r4_index_instrument_diagnostics.csv")
    df = pd.read_csv(path)
    return df[["rsid", "ref", "alt", "beta_I", "se_I"]].copy()


def main():
    log("=" * 72)
    log("R5-03 PRESSO outlier screen + MR-RAPS benchmark")
    log("=" * 72)

    inst = load_index_instruments()
    log(f"Index instruments (R4 table): {len(inst)} SNPs")
    sample_sizes = C.r1_utils.load_finngen_sample_sizes()

    presso_rows = []
    presso_summary = []
    raps_rows = []

    for out in OUTCOMES:
        try:
            har = C.harmonise_index_with_outcome_robust(
                inst, C.FINNGEN_OUTCOMES[out], sample_sizes=sample_sizes, n_exp=N_EXP)
            if har is None or har.empty:
                log(f"{out}: harmonisation empty")
                continue
            har = har.reset_index(drop=True)
            log(f"{out}: harmonised n={len(har)}")
            bx = har["beta"].values.astype(float)
            by = har["beta_outcome"].values.astype(float)
            sy = har["se_outcome"].values.astype(float)
            n = len(har)

            # ---- (a) PRESSO ----
            keep, iters = C.r1_utils.mr_presso_outliers(bx, by, sy, alpha=ALPHA)
            b_full, _, _ = C.r1_utils.mr_ivw(bx, by, sy, random=False)
            if keep.sum() >= 2:
                b_corr, se_corr, p_corr = C.r1_utils.mr_ivw(bx[keep], by[keep], sy[keep])
            else:
                b_corr, se_corr, p_corr = b_full, np.nan, np.nan
            resid = by - bx * b_corr
            z = resid / sy
            p_raw = 2 * (1 - C.stats.norm.cdf(np.abs(z)))
            p_bonf = np.clip(p_raw * n, 0, 1)
            for i in range(n):
                presso_rows.append({
                    "outcome": out, "snp": har["snp"].iloc[i],
                    "beta": bx[i], "se": har["se"].iloc[i],
                    "beta_outcome": by[i], "se_outcome": sy[i],
                    "resid_z": z[i], "p_raw": p_raw[i], "p_bonferroni": p_bonf[i],
                    "is_outlier": bool(not keep[i]), "presso_iterations": iters,
                })
            presso_summary.append({
                "outcome": out, "n_snps": n, "n_outliers": int((~keep).sum()),
                "ivw_beta_full": b_full, "ivw_beta_corrected": b_corr,
                "ivw_se_corrected": se_corr, "ivw_p_corrected": p_corr,
                "or_corrected": np.exp(b_corr), "iterations": iters,
            })
            log(f"  PRESSO: outliers={int((~keep).sum())}/{n}, iters={iters}, "
                f"IVW beta {b_full:.4f} -> {b_corr:.4f} (p={p_corr:.3g})")

            # ---- (b) MR-RAPS vs others ----
            b, se, p = C.r1_utils.mr_ivw(bx, by, sy, random=False)
            raps_rows.append({"outcome": out, "method": "IVW_fixed", "n": n, "beta": b,
                              "se": se, "p": p, "or": np.exp(b),
                              "or_lower": np.exp(b - 1.96 * se), "or_upper": np.exp(b + 1.96 * se),
                              "phi": 1.0})
            b, se, p = C.r1_utils.mr_ivw(bx, by, sy, random=True)
            raps_rows.append({"outcome": out, "method": "IVW_random", "n": n, "beta": b,
                              "se": se, "p": p, "or": np.exp(b),
                              "or_lower": np.exp(b - 1.96 * se), "or_upper": np.exp(b + 1.96 * se),
                              "phi": np.nan})
            try:
                slope, se_s, p_s, icpt, p_i = C.r1_utils.mr_egger(bx, by, sy)
                raps_rows.append({"outcome": out, "method": "MR_Egger", "n": n, "beta": slope,
                                  "se": se_s, "p": p_s, "or": np.exp(slope),
                                  "or_lower": np.exp(slope - 1.96 * se_s),
                                  "or_upper": np.exp(slope + 1.96 * se_s), "phi": np.nan,
                                  "egger_intercept": icpt, "egger_intercept_p": p_i})
            except Exception:
                log(f"  Egger failed:\n{traceback.format_exc()}")
            try:
                b, se, p = C.r1_utils.weighted_median(bx, by, sy)
                raps_rows.append({"outcome": out, "method": "weighted_median", "n": n,
                                  "beta": b, "se": se, "p": p, "or": np.exp(b),
                                  "or_lower": np.exp(b - 1.96 * se),
                                  "or_upper": np.exp(b + 1.96 * se), "phi": np.nan})
            except Exception:
                log(f"  weighted median failed:\n{traceback.format_exc()}")
            try:
                b, se, p, phi, it = C.mr_raps_huber(bx, by, sy, k=1.345)
                raps_rows.append({"outcome": out, "method": "MR_RAPS_huber", "n": n,
                                  "beta": b, "se": se, "p": p, "or": np.exp(b),
                                  "or_lower": np.exp(b - 1.96 * se),
                                  "or_upper": np.exp(b + 1.96 * se), "phi": phi,
                                  "irls_iterations": it})
                log(f"  RAPS: beta={b:.4f} se={se:.4f} p={p:.3g} phi={phi:.3f} iters={it}")
            except Exception:
                log(f"  RAPS failed:\n{traceback.format_exc()}")
        except Exception:
            log(f"{out}: FAILED\n{traceback.format_exc()}")

    if presso_rows:
        df = pd.DataFrame(presso_rows)
        path = os.path.join(C.TABLES_DIR, "r5_presso_outliers.csv")
        df.to_csv(path, index=False)
        log(f"Saved {path} ({len(df)} rows)")
    if presso_summary:
        df = pd.DataFrame(presso_summary)
        path = os.path.join(C.TABLES_DIR, "r5_presso_summary.csv")
        df.to_csv(path, index=False)
        log(f"Saved {path}")
    if raps_rows:
        df = pd.DataFrame(raps_rows)
        path = os.path.join(C.TABLES_DIR, "r5_raps_vs_ivw.csv")
        df.to_csv(path, index=False)
        log(f"Saved {path} ({len(df)} rows)")

    with open(os.path.join(C.LOGS_DIR, "r5_03_run.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    log("R5-03 complete.")


if __name__ == "__main__":
    main()
