#!/usr/bin/env python3
"""Shared helpers for the R5 revision (R5-02/03/04).

Reuses src/r1/utils.py and src/r4/r4_utils.py; does not re-implement MR methods.
All heavy scans are streaming and can be cached as small tables under
results/r5/tables/_cache/ so interrupted runs resume cheaply.
"""

import os
import sys
import json
import warnings
from typing import Dict, List, Optional, Tuple

# Keep BLAS thread pools tiny: several R5 agents share this machine.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src", "r1"))
sys.path.insert(0, os.path.join(ROOT, "src", "r3"))
sys.path.insert(0, os.path.join(ROOT, "src", "r4"))

import utils as r1_utils          # noqa: E402
import r3_utils                   # noqa: E402
import r4_utils                   # noqa: E402

DATA_DIR = os.path.join(ROOT, "data", "real")
R5_OUT = os.path.join(ROOT, "results", "r5")
R5_TABLES = os.path.join(R5_OUT, "tables")
R5_FIGURES = os.path.join(R5_OUT, "figures")
R5_LOGS = os.path.join(R5_OUT, "logs")
CACHE_DIR = os.path.join(R5_TABLES, "_cache")
for _d in (R5_TABLES, R5_FIGURES, R5_LOGS, CACHE_DIR):
    os.makedirs(_d, exist_ok=True)

RNG_SEED = 42

# ---------------------------------------------------------------------------
# Study definitions
# ---------------------------------------------------------------------------
FINNGEN_OUTCOMES = {
    "ALLERG_ASTHMA": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_ALLERG_ASTHMA.gz"),
    "ALLERG_RHINITIS": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_ALLERG_RHINITIS.gz"),
    "L12_ATOPIC": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_L12_ATOPIC.gz"),
}
# FinnGen R12 manifest sample sizes (cases + controls)
FINNGEN_N = {
    "ALLERG_ASTHMA": 13450 + 270290,
    "ALLERG_RHINITIS": 15569 + 474650,
    "L12_ATOPIC": 31245 + 432874,
}

# Pan-UKBB blood components (labels verified against phenotype_manifest.tsv.bgz).
# NOTE: the R5 task brief listed 30120/30130/30180/30190 for the four fractions;
# the manifest shows 30120=lymphocyte COUNT, 30130=monocyte COUNT, 30180=lymphocyte
# pct, 30190=monocyte pct, 30200=neutrophil pct, 30210=eosinophil pct.  We follow
# the manifest (identical to src/r3/r3_utils.PANUKBB_COMPONENTS used to build the
# R4 index) so R5 stays consistent with the R4 index construction.
PANUKBB_COMPONENTS = {
    "crp":            ("biomarkers-30710-both_sexes-irnt.tsv.bgz", 400094),
    "wbc":            ("continuous-30000-both_sexes-irnt.tsv.bgz", 407990),
    "neutrophil_pct": ("continuous-30200-both_sexes-irnt.tsv.bgz", 407282),
    "lymphocyte_pct": ("continuous-30180-both_sexes-irnt.tsv.bgz", 407282),
    "monocyte_pct":   ("continuous-30190-both_sexes-irnt.tsv.bgz", 407270),
    "eosinophil_pct": ("continuous-30210-both_sexes-irnt.tsv.bgz", 407270),
}

ASTLE_TRAITS = {
    "astle_eosinophil": os.path.join(DATA_DIR, "gwas_catalog", "eo_N172275_ukbb_ukbil_meta.tsv.gz"),
    "astle_wbc": os.path.join(DATA_DIR, "gwas_catalog", "wbc_N172435_ukbb_ukbil_meta.tsv.gz"),
}
# Astle component each trait replicates
ASTLE_COMPONENT_MAP = {"astle_eosinophil": "eosinophil_pct", "astle_wbc": "wbc"}


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Small-table caching
# ---------------------------------------------------------------------------
def cache_get(name: str) -> Optional[pd.DataFrame]:
    path = os.path.join(CACHE_DIR, name)
    if os.path.exists(path):
        return pd.read_csv(path)
    return None


def cache_put(name: str, df: pd.DataFrame) -> None:
    path = os.path.join(CACHE_DIR, name)
    df.to_csv(path, index=False)


# ---------------------------------------------------------------------------
# FinnGen instrument selection (allergy as exposure)
# ---------------------------------------------------------------------------
def select_finngen_instruments(outcome: str, p_thresh: float = 5e-8) -> pd.DataFrame:
    """Genome-wide significant SNPs for a FinnGen allergy outcome (as exposure)."""
    cache_name = f"r5_cache_finngen_gwas_{outcome}_{p_thresh:.0e}.csv"
    cached = cache_get(cache_name)
    if cached is not None:
        return cached
    path = FINNGEN_OUTCOMES[outcome]
    usecols = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
               "pval", "beta", "sebeta", "af_alt"]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000, low_memory=False,
                             usecols=lambda c: c in usecols):
        chunk = chunk.dropna(subset=["pval", "beta", "sebeta", "rsids"])
        chunk = chunk[chunk["pval"] < p_thresh]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    df = df.rename(columns={"#chrom": "chrom", "rsids": "snp", "pval": "p",
                            "sebeta": "se"})
    # multi-allelic rsid lists -> take the first rsid alias
    df["snp"] = df["snp"].astype(str).str.split(",").str[0]
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df["ea"] = df["alt"]
    df["oa"] = df["ref"]
    df = df.drop_duplicates("snp", keep="first").reset_index(drop=True)
    df = df.sort_values("p").reset_index(drop=True)
    cache_put(cache_name, df)
    return df


def _proximity_thin(df: pd.DataFrame, kb: int) -> pd.DataFrame:
    """Greedy p-ordered thinning by genomic proximity (no LD matrix).

    Keeps a SNP only if no already-kept SNP on the same chromosome is within
    `kb` kilobases.  Conservative pre-filter that makes the exact LD clump
    memory-feasible for ultra-polygenic traits (tens of thousands of hits).
    """
    import bisect
    window = kb * 1000
    kept_pos: Dict[str, List[int]] = {}
    keep_idx = []
    for i, row in df.sort_values("p").iterrows():
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        positions = kept_pos.setdefault(chrom, [])
        j = bisect.bisect_left(positions, pos)
        clash = False
        if j < len(positions) and abs(positions[j] - pos) <= window:
            clash = True
        if j > 0 and abs(positions[j - 1] - pos) <= window:
            clash = True
        if not clash:
            bisect.insort(positions, pos)
            keep_idx.append(i)
    return df.loc[keep_idx].reset_index(drop=True)


def clump_instruments(df: pd.DataFrame, r2_thresh: float = 0.001,
                      kb: int = 10000, max_for_ld: int = 12000) -> pd.DataFrame:
    """LD clump via r1_utils.ld_clump; preserves original (GRCh38) coordinates.

    For very large candidate sets (>max_for_ld), a proximity-only greedy
    pre-thin (same window) is applied first so that the genotype-matrix LD
    clump stays memory-feasible.
    """
    if df.empty:
        return df
    coord = df[["snp", "chrom", "pos"]].copy()
    work = df
    if len(df) > max_for_ld:
        work = _proximity_thin(df, kb)
        log(f"  clump_instruments: {len(df)} candidates proximity-thinned to "
            f"{len(work)} (window={kb}kb) before exact LD clump")
    clumped = r1_utils.ld_clump(work, r2_thresh=r2_thresh, kb=kb)
    if clumped.empty:
        return clumped
    # ld_clump overwrites chrom/pos with 1KG (GRCh37) values; restore originals.
    clumped = clumped.drop(columns=["chrom", "pos"], errors="ignore")
    clumped = clumped.merge(coord, on="snp", how="left")
    return clumped.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Pan-UKBB component scan (parallel over files, streaming per file)
# ---------------------------------------------------------------------------
def scan_panukbb_components(keys: List[str], components: Optional[List[str]] = None,
                            workers: int = 2, cache_name: Optional[str] = None) -> pd.DataFrame:
    """Scan Pan-UKBB component files for chr:pos:ref:alt keys.

    Spawns stdlib-only child processes (_panukbb_scanner.py), at most `workers`
    concurrently, each writing JSONL to the cache dir.  Robust to memory
    pressure because children never import pandas/numpy.
    """
    if cache_name is not None:
        cached = cache_get(cache_name)
        if cached is not None:
            return cached
    import json as _json
    import subprocess as _sp
    import time as _time
    comps = components or list(PANUKBB_COMPONENTS.keys())
    # unique per-call tag so concurrent runs never clobber each other's
    # keys file or intermediate jsonl outputs.
    tag = f"{os.getpid()}_{int(_time.time())}"
    keys_file = os.path.join(CACHE_DIR, f"r5_cache_panukbb_query_keys_{tag}.json")
    with open(keys_file, "w") as fh:
        _json.dump(sorted(set(keys)), fh)
    py = sys.executable
    scanner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_panukbb_scanner.py")

    def run_one(comp: str) -> str:
        out_jsonl = os.path.join(CACHE_DIR, f"r5_cache_scan_{comp}_{tag}.jsonl")
        fname = PANUKBB_COMPONENTS[comp][0]
        proc = _sp.run([py, scanner, comp, fname, keys_file, out_jsonl],
                       capture_output=True, text=True)
        if proc.returncode != 0:
            log(f"    [scanner {comp}] FAILED rc={proc.returncode}: "
                f"{proc.stderr.strip()[-400:]}")
            return ""
        log(f"    [scanner {comp}] {proc.stdout.strip()}")
        return out_jsonl

    out_files: List[str] = []
    if workers and workers > 1 and len(comps) > 1:
        procs = {}
        for comp in comps:
            out_jsonl = os.path.join(CACHE_DIR, f"r5_cache_scan_{comp}_{tag}.jsonl")
            fname = PANUKBB_COMPONENTS[comp][0]
            p = _sp.Popen([py, scanner, comp, fname, keys_file, out_jsonl],
                          stdout=_sp.PIPE, stderr=_sp.PIPE, text=True)
            procs[comp] = (p, out_jsonl)
            while len(procs) >= workers:
                for c, (pp, of) in list(procs.items()):
                    if pp.poll() is not None:
                        so, se = pp.communicate()
                        if pp.returncode == 0:
                            log(f"    [scanner {c}] {so.strip()}")
                            out_files.append(of)
                        else:
                            log(f"    [scanner {c}] FAILED: {se.strip()[-400:]}")
                        del procs[c]
                if procs:
                    import time as _t
                    _t.sleep(5)
        for c, (pp, of) in procs.items():
            so, se = pp.communicate()
            if pp.returncode == 0:
                log(f"    [scanner {c}] {so.strip()}")
                out_files.append(of)
            else:
                log(f"    [scanner {c}] FAILED: {se.strip()[-400:]}")
    else:
        for comp in comps:
            of = run_one(comp)
            if of:
                out_files.append(of)

    records: List[dict] = []
    for of in out_files:
        if not os.path.exists(of):
            continue
        with open(of) as fh:
            for line in fh:
                records.append(_json.loads(line))
    df = pd.DataFrame(records)
    if df.empty:
        log("    WARNING: Pan-UKBB scan returned 0 rows")
        return df
    if cache_name is not None:
        cache_put(cache_name, df)
    return df


def instrument_keys(df: pd.DataFrame) -> Tuple[List[str], Dict[str, str]]:
    """Both-orientation coordinate keys for an instrument df + key->rsid map."""
    keys, k2r = [], {}
    for _, r in df.iterrows():
        kf = f"{r['chrom']}:{int(r['pos'])}:{r['ref']}:{r['alt']}"
        kr = f"{r['chrom']}:{int(r['pos'])}:{r['alt']}:{r['ref']}"
        keys.extend([kf, kr])
        k2r[kf] = r["snp"]
        k2r[kr] = r["snp"]
    return keys, k2r


def bim_rsid_map(rsids) -> Dict[str, Tuple[str, int, str, str]]:
    """rsid -> (chrom, pos GRCh37, a1, a2) from the 1000G EUR bim."""
    want = set(map(str, rsids))
    out = {}
    with open(r1_utils.BIM_PATH) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            if p[1] in want:
                out[p[1]] = (p[0].replace("chr", ""), int(p[3]), p[4].upper(), p[5].upper())
                if len(out) == len(want):
                    break
    return out


def instrument_keys_grch37(df: pd.DataFrame, rsid_col: str = "snp") -> Tuple[List[str], Dict[str, str], int]:
    """Both-orientation GRCh37 keys via the 1KG bim (Pan-UKBB is GRCh37).

    Pan-UKBB summary statistics are in GRCh37 coordinates (verified: rs429358
    at 19:45411941), whereas FinnGen R12 is GRCh38.  The 1KG bim is GRCh37, so
    rsid -> bim coordinate -> Pan-UKBB key is the correct bridge.  Allele
    letters are build-independent, so harmonise_pair still aligns effects by
    allele.  Returns (keys, key->rsid, n_snps_with_bim_coord).
    """
    rsids = df[rsid_col].astype(str).tolist()
    bmap = bim_rsid_map(rsids)
    keys, k2r = [], {}
    for rsid, (chrom, pos, a1, a2) in bmap.items():
        kf = f"{chrom}:{pos}:{a1}:{a2}"
        kr = f"{chrom}:{pos}:{a2}:{a1}"
        keys.extend([kf, kr])
        k2r[kf] = rsid
        k2r[kr] = rsid
    return keys, k2r, len(bmap)


def map_panukbb_to_rsid(panu_df: pd.DataFrame, k2r: Dict[str, str]) -> pd.DataFrame:
    out = panu_df.copy()
    out["rsid"] = out["snp"].map(k2r)
    out = out.dropna(subset=["rsid"])
    out = out.drop_duplicates(["component", "rsid"], keep="first")
    return out


# ---------------------------------------------------------------------------
# Astle scans
# ---------------------------------------------------------------------------
ASTLE_COLS = ["SNP", "CHR", "POS", "A1", "A2", "EAF", "Beta", "se", "P", "N"]


def scan_astle_rsids(trait: str, rsids: set, cache_name: Optional[str] = None) -> pd.DataFrame:
    """Extract Astle rows for a set of rsids (as OUTCOME: alt=A1 effect allele)."""
    if cache_name is not None:
        cached = cache_get(cache_name)
        if cached is not None:
            return cached
    path = ASTLE_TRAITS[trait]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000, low_memory=False,
                             usecols=lambda c: c in ASTLE_COLS):
        chunk = chunk[chunk["SNP"].isin(rsids)]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    df = standardise_astle(df)
    if cache_name is not None:
        cache_put(cache_name, df)
    return df


def standardise_astle(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(columns={"SNP": "snp", "CHR": "chrom", "POS": "pos",
                             "Beta": "beta", "P": "p", "EAF": "af_alt", "N": "n"})
    out["chrom"] = out["chrom"].astype(str)
    out["A1"] = out["A1"].str.upper()
    out["A2"] = out["A2"].str.upper()
    # effect allele = A1
    out["alt"] = out["A1"]
    out["ref"] = out["A2"]
    out["ea"] = out["A1"]
    out["oa"] = out["A2"]
    out = out.drop_duplicates("snp", keep="first").reset_index(drop=True)
    return out


def select_astle_instruments(trait: str, p_thresh: float = 5e-8,
                             cache_name: Optional[str] = None) -> pd.DataFrame:
    """Genome-wide significant Astle SNPs (as EXPOSURE)."""
    if cache_name is not None:
        cached = cache_get(cache_name)
        if cached is not None:
            return cached
    path = ASTLE_TRAITS[trait]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000, low_memory=False,
                             usecols=lambda c: c in ASTLE_COLS):
        chunk = chunk.dropna(subset=["P", "Beta", "se"])
        chunk = chunk[chunk["P"] < p_thresh]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    df = standardise_astle(df)
    df = df.sort_values("p").reset_index(drop=True)
    if cache_name is not None:
        cache_put(cache_name, df)
    return df


# ---------------------------------------------------------------------------
# FinnGen outcome extraction by rsid
# ---------------------------------------------------------------------------
def finngen_outcome_by_rsids(outcome: str, rsids: set,
                             cache_name: Optional[str] = None) -> pd.DataFrame:
    if cache_name is not None:
        cached = cache_get(cache_name)
        if cached is not None:
            return cached
    df = r1_utils.load_finngen(FINNGEN_OUTCOMES[outcome], rsids=set(rsids))
    if df.empty:
        return df
    df["snp"] = df["snp"].astype(str).str.split(",").str[0]
    df = df.drop_duplicates("snp", keep="first").reset_index(drop=True)
    if cache_name is not None:
        cache_put(cache_name, df)
    return df


def finngen_coords_for_rsids(rsids: set, outcome: str = "ALLERG_ASTHMA",
                             cache_name: Optional[str] = None) -> pd.DataFrame:
    """Use FinnGen (GRCh38) as an rsid->coordinate dictionary (no betas used)."""
    if cache_name is not None:
        cached = cache_get(cache_name)
        if cached is not None:
            return cached
    path = FINNGEN_OUTCOMES[outcome]
    usecols = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes"]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000, low_memory=False,
                             usecols=lambda c: c in usecols):
        chunk = chunk.dropna(subset=["rsids"])
        chunk["rsid1"] = chunk["rsids"].astype(str).str.split(",").str[0]
        chunk = chunk[chunk["rsid1"].isin(rsids)]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    df = df.rename(columns={"#chrom": "chrom", "rsid1": "snp"})
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df = df.drop_duplicates("snp", keep="first").reset_index(drop=True)
    if cache_name is not None:
        cache_put(cache_name, df)
    return df


def nearest_genes_for_rsids(rsids: set, outcome: str = "ALLERG_ASTHMA") -> Dict[str, str]:
    df = finngen_coords_for_rsids(rsids, outcome=outcome,
                                  cache_name=None)
    if df.empty or "nearest_genes" not in df.columns:
        return {}
    return dict(zip(df["snp"], df["nearest_genes"].fillna("").astype(str)))


# ---------------------------------------------------------------------------
# MR battery (reuses r1 utils implementations)
# ---------------------------------------------------------------------------
def mr_battery(har: pd.DataFrame, do_presso: bool = True) -> List[dict]:
    """IVW fixed/random, MR-Egger, weighted median, PRESSO-corrected IVW."""
    bx = har["beta"].values.astype(float)
    by = har["beta_outcome"].values.astype(float)
    sy = har["se_outcome"].values.astype(float)
    k = len(har)
    rows = []

    def add(method, b, se, p, **extra):
        row = {"method": method, "n_snps": k, "beta": b, "se": se, "p": p,
               "or": np.exp(b), "or_lower": np.exp(b - 1.96 * se),
               "or_upper": np.exp(b + 1.96 * se)}
        row.update(extra)
        rows.append(row)

    if k == 1:
        b, se, p = r1_utils.mr_wald_ratio(bx[0], by[0], sy[0])
        add("Wald_ratio", b, se, p)
        return rows
    b, se, p = r1_utils.mr_ivw(bx, by, sy, random=False)
    add("IVW_fixed", b, se, p)
    b, se, p = r1_utils.mr_ivw(bx, by, sy, random=True)
    add("IVW_random", b, se, p)
    Q, Qp = r1_utils.cochran_q(bx, by, sy, b)
    slope, se_slope, p_slope, intercept, p_int = r1_utils.mr_egger(bx, by, sy)
    add("MR_Egger", slope, se_slope, p_slope,
        egger_intercept=intercept, egger_intercept_p=p_int, cochran_Q=Q, cochran_Q_p=Qp)
    b, se, p = r1_utils.weighted_median(bx, by, sy)
    add("weighted_median", b, se, p)
    if do_presso:
        keep, n_iter = r1_utils.mr_presso_outliers(bx, by, sy, alpha=0.05)
        n_removed = int(k - keep.sum())
        if 2 <= keep.sum() < k:
            b2, se2, p2 = r1_utils.mr_ivw(bx[keep], by[keep], sy[keep], random=False)
            add("IVW_PRESSO_corrected", b2, se2, p2,
                presso_outliers_removed=n_removed)
        elif n_removed == 0:
            add("IVW_PRESSO_corrected", b, se, p, presso_outliers_removed=0)
    return rows


def harmonise_for_mr(exposure_df: pd.DataFrame, outcome_df: pd.DataFrame,
                     n_exp: int, n_out: int, apply_steiger: bool = True) -> Tuple[pd.DataFrame, dict]:
    """harmonise_pair + optional Steiger filter. Returns (har, diagnostics)."""
    exp = exposure_df.copy()
    out = outcome_df.copy().drop_duplicates("snp", keep="first")
    har = r1_utils.harmonise_pair(exp, out)
    diag = {"n_harmonised": 0 if har is None else len(har),
            "n_after_steiger": 0, "steiger_applied": False}
    if har is None or har.empty:
        return pd.DataFrame(), diag
    if apply_steiger:
        har_s = r1_utils.steiger_filter(har, n_exp=n_exp, n_out=n_out)
        har_pass = har_s[har_s["steiger_pass"]]
        diag["n_after_steiger"] = len(har_pass)
        if len(har_pass) >= 3:
            diag["steiger_applied"] = True
            return har_pass.reset_index(drop=True), diag
        diag["steiger_fallback_unfiltered"] = True
    diag["n_after_steiger"] = len(har)
    return har.reset_index(drop=True), diag
