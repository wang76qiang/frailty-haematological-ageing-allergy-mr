#!/usr/bin/env python3
"""
R5-12 completion: fill the OneK1K replication table for the 5 cells whose
streaming failed with OOM during the concurrent first pass (cd4et, cd4nc,
cd8et, monoc, nk).  Memory-frugal: 100k-row chunks, one cell at a time.

Reads  results/r5/tables/r5_celltype_mr_matrix.csv  (significant pairs)
Writes results/r5/tables/r5_celltype_replication.csv (all 19 pairs)
       results/r5/logs/r5_12b_fill.json
"""
import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "r1"))
from utils import load_onek1k  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_TABLES = os.path.join(BASE, "results", "r5", "tables")
OUT_LOGS = os.path.join(BASE, "results", "r5", "logs")
CACHE = os.path.join(OUT_LOGS, "cache")
MATRIX = os.path.join(OUT_TABLES, "r5_celltype_mr_matrix.csv")
REP_OUT = os.path.join(OUT_TABLES, "r5_celltype_replication.csv")

DICE2ONEK = {
    "TH1": "cd4et", "TH2": "cd4et", "TH17": "cd4et", "TFH": "cd4et",
    "THSTAR": "cd4et", "CD4_STIM": "cd4et", "TREG_MEM": "cd4et",
    "CD4_NAIVE": "cd4nc", "TREG_NAIVE": "cd4nc",
    "CD8_NAIVE": "cd8et", "CD8_STIM": "cd8et",
    "MONOCYTES": "monoc", "M2": "monoc",
    "NK": "nk", "B_CELL_NAIVE": "bin",
}

LOG = []


def log(m):
    print(m, flush=True)
    LOG.append(m)


def main():
    mr = pd.read_csv(MATRIX)
    sig = mr[mr["significant"] == True]  # noqa: E712
    pairs = sig[["cell", "gene", "snp_mr", "ea", "oa", "wald_beta", "outcome"]]
    pairs = pairs.rename(columns={"snp_mr": "snp"}).drop_duplicates()
    log(f"significant pairs to replicate: {len(pairs)}")

    rows = []
    for ocell in sorted({DICE2ONEK[c] for c in pairs["cell"] if c in DICE2ONEK}):
        cache = os.path.join(CACHE, f"onek1k_{ocell}.csv")
        cell_pairs = pairs[pairs["cell"].map(DICE2ONEK) == ocell]
        genes = sorted(cell_pairs["gene"].unique())
        df = pd.DataFrame()
        try:
            if os.path.exists(cache):
                df = pd.read_csv(cache, dtype={"chrom": str})
                df = df[df["gene"].isin(genes)]
                log(f"  {ocell}: {len(df)} rows from cache")
            else:
                df = load_onek1k(ocell, genes, chunksize=100_000)
                if not df.empty:
                    tmp = cache + ".tmp"
                    df.to_csv(tmp, index=False)
                    os.replace(tmp, cache)
                log(f"  {ocell}: streamed {len(df)} rows (genes={genes})")
        except Exception as e:
            log(f"  WARN {ocell} load failed: {e}")
        for _, pr in cell_pairs.iterrows():
            hit = df[(df["gene"] == pr["gene"]) & (df["snp"] == pr["snp"])] \
                if not df.empty else pd.DataFrame()
            if hit.empty:
                rows.append({**pr.to_dict(), "onek1k_cell": ocell,
                             "onek1k_beta": np.nan, "onek1k_p": np.nan,
                             "onek1k_beta_aligned": np.nan,
                             "allele_aligned": None, "same_direction": None,
                             "note": "SNP not found in OneK1K" if not df.empty
                                     else "cell load failed/empty"})
                continue
            h = hit.iloc[0]
            ea = str(pr.get("ea", "")).upper()
            a1, a2 = str(h["ea"]).upper(), str(h["oa"]).upper()
            if ea and ea == a2:
                ok_beta, aligned = h["beta"], True
            elif ea and ea == a1:
                ok_beta, aligned = -h["beta"], True
            else:
                ok_beta, aligned = h["beta"], False
            same = bool(np.sign(ok_beta) == np.sign(pr["wald_beta"])) \
                if aligned else None
            rows.append({**pr.to_dict(), "onek1k_cell": ocell,
                         "onek1k_beta": h["beta"], "onek1k_p": h["p"],
                         "onek1k_beta_aligned": ok_beta,
                         "allele_aligned": aligned, "same_direction": same,
                         "note": "" if aligned else "allele mismatch"})

    rep = pd.DataFrame(rows)
    rep.to_csv(REP_OUT, index=False)
    valid = rep[rep["same_direction"].notna()]
    n_same = int(valid["same_direction"].sum()) if len(valid) else 0
    p_sign = stats.binomtest(n_same, len(valid), 0.5).pvalue if len(valid) else np.nan
    log(f"replication rows={len(rep)}; testable={len(valid)}; "
        f"same-direction={n_same}/{len(valid)} (binomial p={p_sign})")
    with open(os.path.join(OUT_LOGS, "r5_12b_fill.json"), "w") as f:
        json.dump({"log": LOG, "n_pairs": len(rep), "n_testable": len(valid),
                   "n_same_direction": n_same, "binom_p": p_sign},
                  f, indent=2, default=str)
    log("DONE 12b")


if __name__ == "__main__":
    main()
