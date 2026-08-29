"""Importable worker for R5-09 pQTL file scanning (Windows spawn-safe)."""

import gc
import json
import os
import time

import pandas as pd

PQTL_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
             "beta", "standard_error", "effect_allele_frequency", "p_value",
             "variant_id", "rsid", "n"]


def _scan_once(path, chrom, tss, region_window, build_window, chunksize):
    region_rows = []
    n5mb = 0
    best = None
    reader = pd.read_csv(path, sep="\t", usecols=PQTL_COLS, chunksize=chunksize,
                         dtype={"chromosome": str})
    for chunk in reader:
        chunk["chromosome"] = chunk["chromosome"].str.replace("chr", "", regex=False)
        sub = chunk[chunk["chromosome"] == str(chrom)]
        if not sub.empty:
            pos = sub["base_pair_location"]
            n5mb += int(((pos >= tss - build_window) & (pos <= tss + build_window)).sum())
            keep = sub[(pos >= tss - region_window) & (pos <= tss + region_window)]
            if not keep.empty:
                region_rows.append(keep)
        idx = chunk["p_value"].idxmin()
        row = chunk.loc[idx]
        if best is None or row["p_value"] < best[0]:
            best = (row["p_value"], {
                "chromosome": row["chromosome"],
                "base_pair_location": int(row["base_pair_location"]),
                "effect_allele": row["effect_allele"],
                "other_allele": row["other_allele"],
                "beta": row["beta"], "standard_error": row["standard_error"],
                "effect_allele_frequency": row["effect_allele_frequency"],
                "p_value": row["p_value"], "variant_id": row["variant_id"],
                "rsid": row["rsid"], "n": row["n"]})
        del chunk
    region = (pd.concat(region_rows, ignore_index=True)
              if region_rows else pd.DataFrame(columns=PQTL_COLS))
    return region, n5mb, best


def process_pqtl_file(args):
    """Scan one pQTL file; cache TSS+/-1.1Mb SNPs + JSON summary.

    args = (protein, path, chrom, tss, cache_dir, region_window, build_window)
    Resumable: returns early if cache + JSON already exist.
    Never raises; errors returned in dict.
    """
    protein, path, chrom, tss, cache_dir, region_window, build_window = args
    out = {"protein": protein, "ok": False}
    os.makedirs(cache_dir, exist_ok=True)
    region_csv = os.path.join(cache_dir, f"{protein}_region.csv")
    gmin_csv = os.path.join(cache_dir, f"{protein}_globalmin.csv")
    json_path = os.path.join(cache_dir, f"{protein}_scan.json")
    if os.path.exists(json_path) and os.path.exists(region_csv):
        try:
            with open(json_path) as f:
                cached = json.load(f)
            cached["protein"] = protein
            cached["cached"] = True
            return cached
        except Exception:  # noqa: BLE001
            pass
    for attempt in range(3):
        try:
            region, n5mb, best = _scan_once(path, chrom, tss, region_window,
                                            build_window, chunksize=250_000)
            region.to_csv(region_csv, index=False)
            pd.DataFrame([best[1]]).to_csv(gmin_csv, index=False)
            out.update({"ok": True, "n5mb": n5mb, "n_region": int(len(region)),
                        "global_min_p": float(best[0]), "global_min": best[1]})
            with open(json_path, "w") as f:
                json.dump({k: v for k, v in out.items() if k != "global_min"}, f)
            break
        except MemoryError:
            gc.collect()
            time.sleep(5 * (attempt + 1))
            out["error"] = "MemoryError"
        except Exception as e:  # noqa: BLE001
            import traceback
            out["error"] = repr(e)[:160] + " @ " + traceback.format_exc().splitlines()[-3][:120]
            gc.collect()
            time.sleep(3 * (attempt + 1))
    return out
