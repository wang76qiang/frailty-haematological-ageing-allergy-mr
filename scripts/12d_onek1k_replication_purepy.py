#!/usr/bin/env python3
"""
R5-12 (completion): authoritative OneK1K replication for the 33 FDR-significant
DICE cell x gene MR pairs.

Why this exists: every pandas-based attempt to stream the 1.5-2.3 GB OneK1K
`*_eqtl_table.tsv.gz` files died with out-of-memory under the current RAM
pressure (even 6.87 MiB allocations fail).  This script streams the gz files
with the plain `gzip` + line-split path (a few MB of RAM, no pandas chunk
allocation) and keeps only rows whose GENE is one of the replication targets.

Replication definition (unchanged from the first-pass script):
  * exact same SNP present in the OneK1K full eQTL table for that gene;
  * OneK1K beta = SPEARMANS_RHO (reported w.r.t. allele A2); allele-aligned to
    the DICE effect allele (flip sign if DICE EA == OneK1K A1);
  * same_direction = sign(aligned OneK1K beta) == sign(DICE eQTL beta)
    (pair-level) and == sign(wald_beta) (row-level, per outcome);
  * binomial sign test over unique testable pairs.

Reads : results/r5/tables/r5_celltype_mr_matrix.csv
Writes: results/r5/tables/r5_celltype_replication.csv   (33 rows)
        results/r5/logs/r5_12d_replication.json
        results/r5/logs/cache/onek1k_full_<cell>.csv    (gene-filtered caches)
"""
import csv
import gzip
import json
import os

import numpy as np
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_TABLES = os.path.join(BASE, "results", "r5", "tables")
OUT_LOGS = os.path.join(BASE, "results", "r5", "logs")
CACHE = os.path.join(OUT_LOGS, "cache")
ONEK = os.path.join(BASE, "data", "real", "onek1k")
MATRIX = os.path.join(OUT_TABLES, "r5_celltype_mr_matrix.csv")
REP_OUT = os.path.join(OUT_TABLES, "r5_celltype_replication.csv")
os.makedirs(CACHE, exist_ok=True)

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


def stream_cell(ocell, genes, rsids):
    """Return {(gene, rsid): dict(rho, p, a1, a2)} using pure-python streaming."""
    cache = os.path.join(CACHE, f"onek1k_full_{ocell}.csv")
    out = {}
    if os.path.exists(cache):
        with open(cache, newline="") as f:
            for r in csv.DictReader(f):
                out[(r["gene"], r["snp"])] = r
        log(f"  {ocell}: {len(out)} gene rows from cache")
        return out
    path = os.path.join(ONEK, f"{ocell}_eqtl_table.tsv.gz")
    if not os.path.exists(path):
        log(f"  {ocell}: full table missing -> skip")
        return out
    n_lines = 0
    with gzip.open(path, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        gi, ri = idx["GENE"], idx["RSID"]
        ia1, ia2 = idx["A1"], idx["A2"]
        irho, ip = idx["SPEARMANS_RHO"], idx["P_VALUE"]
        for line in f:
            n_lines += 1
            if n_lines % 4_000_000 == 0:
                log(f"    {ocell}: {n_lines/1e6:.0f}M lines, {len(out)} kept")
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= irho:
                continue
            g = parts[gi]
            if g not in genes:
                continue
            s = parts[ri]
            if s not in rsids:
                continue  # keep only exact replication SNPs
            out[(g, s)] = {"gene": g, "snp": s, "a1": parts[ia1],
                           "a2": parts[ia2], "rho": parts[irho], "p": parts[ip]}
    # cache for reuse
    tmp = cache + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "snp", "a1", "a2", "rho", "p"])
        w.writeheader()
        w.writerows(out.values())
    os.replace(tmp, cache)
    log(f"  {ocell}: streamed {n_lines} lines, kept {len(out)} exact-SNP rows")
    return out


def main():
    with open(MATRIX, newline="") as f:
        rows = list(csv.DictReader(f))
    mr = [r for r in rows if r["outcome"].strip()]
    sig = [r for r in mr if r["significant"].strip() == "True"]
    # DICE exposure beta per (cell, gene, snp) from grid rows
    dice_beta = {}
    for r in rows:
        if r["instrument"].strip() == "pass" and r["beta"].strip():
            dice_beta[(r["cell"], r["gene"], r["snp_lead"])] = float(r["beta"])
    pairs = []
    seen = set()
    for r in sig:
        key = (r["cell"], r["gene"], r["snp_mr"], r["outcome"])
        if key in seen:
            continue
        seen.add(key)
        pairs.append({"cell": r["cell"], "gene": r["gene"], "snp": r["snp_mr"],
                      "ea": r["ea"], "oa": r["oa"],
                      "wald_beta": float(r["wald_beta"]), "outcome": r["outcome"]})
    log(f"significant pairs (pair x outcome rows): {len(pairs)}")
    genes = {p["gene"] for p in pairs}
    rsids = {p["snp"] for p in pairs}
    log(f"unique genes={sorted(genes)}  unique SNPs={len(rsids)}")

    tables = {}
    for ocell in sorted({DICE2ONEK[p["cell"]] for p in pairs if p["cell"] in DICE2ONEK}):
        tables[ocell] = stream_cell(ocell, genes, rsids)

    out_rows = []
    for p in pairs:
        ocell = DICE2ONEK.get(p["cell"], "")
        hit = tables.get(ocell, {}).get((p["gene"], p["snp"]))
        base = {**p, "onek1k_cell": ocell, "onek1k_snp": "", "onek1k_beta": "",
                "onek1k_p": "", "onek1k_beta_aligned": "", "allele_aligned": "",
                "same_direction": "", "exact_snp": False,
                "note": "exact SNP absent in OneK1K full table"}
        if hit is None:
            out_rows.append(base)
            continue
        rho = float(hit["rho"])
        ea = p["ea"].upper()
        a1, a2 = hit["a1"].upper(), hit["a2"].upper()
        if ea and ea == a2:
            aligned_beta, aligned = rho, True
        elif ea and ea == a1:
            aligned_beta, aligned = -rho, True
        else:
            aligned_beta, aligned = rho, False
        db = dice_beta.get((p["cell"], p["gene"], p["snp"]))
        same = (np.sign(aligned_beta) == np.sign(p["wald_beta"])) if aligned else None
        base.update({"onek1k_snp": p["snp"], "onek1k_beta": rho,
                     "onek1k_p": hit["p"], "onek1k_beta_aligned": aligned_beta,
                     "allele_aligned": aligned,
                     "same_direction": same if same is not None else "",
                     "exact_snp": True,
                     "dice_beta": db if db is not None else "",
                     "note": "" if aligned else "allele mismatch"})
        out_rows.append(base)

    with open(REP_OUT, "w", newline="") as f:
        cols = ["cell", "gene", "snp", "ea", "oa", "wald_beta", "outcome",
                "onek1k_cell", "onek1k_snp", "onek1k_beta", "onek1k_p",
                "onek1k_beta_aligned", "allele_aligned", "same_direction",
                "exact_snp", "dice_beta", "note"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # pair-level sign test (unique cell-gene-snp, direction vs DICE eQTL beta)
    pair_level = {}
    for r in out_rows:
        if r["exact_snp"] and r["allele_aligned"] is True and r["dice_beta"] != "":
            pair_level[(r["cell"], r["gene"], r["snp"])] = (
                np.sign(float(r["onek1k_beta_aligned"])) ==
                np.sign(float(r["dice_beta"])))
    n_test = len(pair_level)
    n_same = int(sum(pair_level.values()))
    p_pair = stats.binomtest(n_same, n_test, 0.5).pvalue if n_test else float("nan")

    # row-level (pair x outcome, direction vs wald_beta)
    valid = [r for r in out_rows if r["same_direction"] in (True, False)]
    nr_same = int(sum(1 for r in valid if r["same_direction"]))
    p_row = stats.binomtest(nr_same, len(valid), 0.5).pvalue if valid else float("nan")

    log(f"exact-SNP testable rows: {len(valid)}/{len(out_rows)}; "
        f"same-direction {nr_same}/{len(valid)} (binomial p={p_row:.4f})")
    log(f"unique testable pairs: {n_test}; same-direction {n_same}/{n_test} "
        f"(binomial p={p_pair:.4f})")
    for (c, g, s), v in sorted(pair_level.items()):
        log(f"    {'SAME ' if v else 'OPPOSITE'}  {c:14s} {g:8s} {s}")

    with open(os.path.join(OUT_LOGS, "r5_12d_replication.json"), "w") as f:
        json.dump({"log": LOG, "n_rows": len(out_rows),
                   "n_exact_testable_rows": len(valid),
                   "row_level_same": nr_same, "row_level_binom_p": p_row,
                   "n_unique_testable_pairs": n_test,
                   "pair_level_same": n_same, "pair_level_binom_p": p_pair},
                  f, indent=2, default=str)
    log("DONE 12d")


if __name__ == "__main__":
    main()
