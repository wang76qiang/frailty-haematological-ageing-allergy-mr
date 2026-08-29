#!/usr/bin/env python3
"""Temporary helper: pre-cache OneK1K exact-SNP rows for 12_celltype_mr_matrix.
Streams only the not-yet-cached cell files. Safe to re-run. Deleted after use."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import importlib.util
spec = importlib.util.spec_from_file_location(
    "ctmr", os.path.join(HERE, "12_celltype_mr_matrix.py"))
ctmr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctmr)

GENES = ctmr.ALL_GENES
# rsids needed = all lead SNPs passing F>10 (recompute quickly from DICE)
import pandas as pd
leads = ctmr.collect_lead_snps()
rsids = leads.loc[leads["instrument"] == "pass", "snp"].unique()
cells = ["bin", "cd4et", "cd4nc", "cd8et", "monoc", "nk"]
for cell in cells:
    cache = os.path.join(ctmr.CACHE_DIR, f"onek1k_full_{cell}.csv")
    if os.path.exists(cache):
        print(f"skip {cell} (cached)", flush=True)
        continue
    print(f"streaming {cell} ...", flush=True)
    df = ctmr.load_onek1k_full(cell, GENES, rsids=rsids)
    print(f"  {cell}: {len(df)} rows", flush=True)
print("DRIVER DONE", flush=True)
