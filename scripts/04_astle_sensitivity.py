#!/usr/bin/env python3
"""
R5-04 Astle instrument-selection sensitivity.

Use the Astle 2016 eosinophil-count and WBC GWAS (UKBB + UK BiLEVE meta,
GRCh37, rsid-based) directly as exposures: select p < 5e-8 SNPs, LD-clump
(r2 < 0.001, 10 Mb, 1KG EUR), and run MR against the 3 FinnGen allergy
outcomes (rsid matching only, since Astle is GRCh37 and FinnGen is GRCh38).
Effect directions are compared with the R4 Pan-UKBB-based index IVW.

Caveat (logged): Astle includes UKBB samples that overlap Pan-UKBB; we do not
re-estimate Pan-UKBB without that overlap because the Pan-UKBB files are
GRCh38 without rsids and cannot be lifted over.
"""

import os
import sys
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r5x_utils as C

OUTCOMES = list(C.FINNGEN_OUTCOMES)
TRAITS = {"eosinophil_astle": "eosinophil count (Astle 2016)",
          "wbc_astle": "WBC count (Astle 2016)"}
P_THRESH = 5e-8

log_lines = []


def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)


def main():
    log("=" * 72)
    log("R5-04 Astle instrument-selection sensitivity")
    log("=" * 72)
    log("NOTE: Astle (UKBB+UK BiLEVE meta) partially overlaps Pan-UKBB samples; "
        "Pan-UKBB cannot be re-estimated without overlap (GRCh38, no rsids -> "
        "no liftover). Overlap inflates two-sample MR Type-I error only mildly; "
        "this analysis is used for direction/magnitude concordance.")

    # forward R4 index IVW for direction comparison
    try:
        r4 = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables", "r4_mr_diagnostics.csv"))
        r4 = r4.set_index("outcome")
    except Exception:
        log(f"R4 diagnostics read FAILED\n{traceback.format_exc()}")
        r4 = pd.DataFrame()

    rows = []
    out_path = os.path.join(C.TABLES_DIR, "r5_astle_selection_sensitivity.csv")

    def save_partial():
        if not rows:
            return
        res = pd.DataFrame(rows)
        cols = ["trait", "outcome", "method", "n", "n_instruments", "p_threshold",
                "beta", "se", "p", "or", "or_lower", "or_upper",
                "egger_intercept", "egger_intercept_p",
                "r4_index_beta", "r4_index_p", "direction_concordant_with_r4"]
        cols = [c for c in cols if c in res.columns]
        res[cols].to_csv(out_path, index=False)

    for trait, desc in TRAITS.items():
        try:
            sig = C.astle_significant(trait, p_thresh=P_THRESH)
            log(f"{trait} ({desc}): {len(sig)} SNPs at p<{P_THRESH}")
            if sig.empty:
                continue
            cl = C.clump_instruments(sig, r2=0.001, kb=10000)
            log(f"  clumped instruments: {len(cl)}")
            if cl.empty or len(cl) < 3:
                log(f"  too few clumped SNPs for {trait}; skipping MR")
                continue
            exp = cl[["snp", "ea", "oa", "beta", "se", "p"]].copy()
            # load FinnGen outcomes once for the union of clumped rsids
            out_frames = {}
            for out in OUTCOMES:
                df = C.retry_call(
                    lambda o=out: C.r1_utils.load_finngen(
                        C.FINNGEN_OUTCOMES[o], rsids=set(cl["snp"]), chunksize=200000),
                    label=f"finngen_{out}")
                out_frames[out] = df
            for out in OUTCOMES:
                try:
                    df = out_frames[out]
                    if df is None or df.empty:
                        log(f"  {trait} x {out}: no FinnGen rsid matches")
                        continue
                    har = C.r1_utils.harmonise_pair(exp, df)
                    if har.empty or len(har) < 3:
                        log(f"  {trait} x {out}: <3 harmonised SNPs (n={len(har)})")
                        continue
                    r4_beta = np.log(r4.loc[out, "ivw_fixed_or"]) if (not r4.empty and out in r4.index) else np.nan
                    r4_p = r4.loc[out, "ivw_fixed_p"] if (not r4.empty and out in r4.index) else np.nan
                    for m in C.mr_battery(har):
                        concordant = (np.sign(m["beta"]) == np.sign(r4_beta)) if (
                            np.isfinite(m["beta"]) and np.isfinite(r4_beta)) else np.nan
                        m.update({"trait": trait, "outcome": out,
                                  "n_instruments": len(cl), "p_threshold": P_THRESH,
                                  "r4_index_beta": r4_beta, "r4_index_p": r4_p,
                                  "direction_concordant_with_r4": concordant})
                        rows.append(m)
                    log(f"  {trait} x {out}: harmonised n={len(har)}")
                except Exception:
                    log(f"  {trait} x {out}: FAILED\n{traceback.format_exc()}")
            save_partial()
            log(f"  {trait}: partial results saved ({len(rows)} rows)")
        except Exception:
            log(f"{trait}: FAILED\n{traceback.format_exc()}")
            save_partial()

    if rows:
        save_partial()
        log(f"Saved {out_path} ({len(rows)} rows)")
        res = pd.DataFrame(rows)
        ivw = res[res["method"] == "IVW_fixed"]
        for _, r in ivw.iterrows():
            log(f"  IVW {r['trait']} -> {r['outcome']}: beta={r['beta']:.4f} "
                f"(p={r['p']:.3g}), R4 index beta={r['r4_index_beta']:.4f}, "
                f"concordant={r['direction_concordant_with_r4']}")
    else:
        log("No results produced.")

    with open(os.path.join(C.LOGS_DIR, "r5_04_run.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    log("R5-04 complete.")


if __name__ == "__main__":
    main()
