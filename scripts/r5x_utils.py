#!/usr/bin/env python3
"""Shared helpers for the R5-02/03/04 MR tasks (self-contained copy).

NOTE: this module is intentionally named uniquely to avoid clashing with any
other r5_common.py in the same directory.
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for sub in ("r1", "r3", "r4"):
    p = os.path.join(ROOT, "src", sub)
    if p not in sys.path:
        sys.path.insert(0, p)

import utils as r1_utils          # noqa: E402
import r3_utils                   # noqa: E402
import r4_utils                   # noqa: E402

DATA_DIR = os.path.join(ROOT, "data", "real")
OUT_DIR = os.path.join(ROOT, "results", "r5")
TABLES_DIR = os.path.join(OUT_DIR, "tables")
FIGURES_DIR = os.path.join(OUT_DIR, "figures")
LOGS_DIR = os.path.join(OUT_DIR, "logs")

for d in (TABLES_DIR, FIGURES_DIR, LOGS_DIR):
    os.makedirs(d, exist_ok=True)

RNG = np.random.default_rng(42)

FINNGEN_OUTCOMES = {
    "ALLERG_ASTHMA": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_ALLERG_ASTHMA.gz"),
    "ALLERG_RHINITIS": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_ALLERG_RHINITIS.gz"),
    "L12_ATOPIC": os.path.join(DATA_DIR, "finngen_full", "finngen_R12_L12_ATOPIC.gz"),
}
# FinnGen R12 sample sizes (cases + controls from manifest)
FINNGEN_N = {"ALLERG_ASTHMA": 283740, "ALLERG_RHINITIS": 490219, "L12_ATOPIC": 464119}

# Pan-UKBB components.  Mapping verified against phenotype_manifest.tsv.bgz:
#   30000 WBC count; 30180 lymphocyte %; 30190 monocyte %;
#   30200 neutrophil %; 30210 eosinophil %; 30710 CRP.
# (NB: 30120/30130 are lymphocyte/monocyte *counts*, not percentages.)
PANUKBB_COMPONENTS = {
    "crp": ("biomarkers-30710-both_sexes-irnt.tsv.bgz", 400094),
    "wbc": ("continuous-30000-both_sexes-irnt.tsv.bgz", 407990),
    "neutrophil_pct": ("continuous-30200-both_sexes-irnt.tsv.bgz", 407282),
    "lymphocyte_pct": ("continuous-30180-both_sexes-irnt.tsv.bgz", 407282),
    "monocyte_pct": ("continuous-30190-both_sexes-irnt.tsv.bgz", 407270),
    "eosinophil_pct": ("continuous-30210-both_sexes-irnt.tsv.bgz", 407270),
}

ASTLE_FILES = {
    "eosinophil_astle": os.path.join(DATA_DIR, "gwas_catalog", "eo_N172275_ukbb_ukbil_meta.tsv.gz"),
    "wbc_astle": os.path.join(DATA_DIR, "gwas_catalog", "wbc_N172435_ukbb_ukbil_meta.tsv.gz"),
}


# ---------------------------------------------------------------------------
# Retry helpers for transient OOM under parallel load
# ---------------------------------------------------------------------------
import time as _time


def retry_call(fn, tries: int = 4, base_sleep: float = 15.0, label: str = ""):
    """Retry fn() on MemoryError/ParserError (transient under parallel load)."""
    last = None
    for t in range(tries):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001
            last = exc
            msg = str(exc)
            if ("MemoryError" not in msg) and ("memory" not in msg.lower()) and t == tries - 1:
                raise
            _time.sleep(base_sleep * (t + 1))
    raise last


def harmonise_index_with_outcome_robust(inst, out_path, sample_sizes, n_exp: int = 400000,
                                        chunksize: int = 100000):
    """Replica of r4_utils.harmonise_index_with_outcome with retry + smaller chunks."""
    rsids = set(inst["rsid"])

    def _load():
        return r1_utils.load_finngen(out_path, rsids=rsids, chunksize=chunksize)

    out_df = retry_call(_load, label="load_finngen")
    if out_df is None or out_df.empty:
        return pd.DataFrame()
    exp_df = inst[["rsid", "ref", "alt", "beta_I", "se_I"]].copy().rename(
        columns={"rsid": "snp", "ref": "oa", "alt": "ea", "beta_I": "beta", "se_I": "se"})
    har = r1_utils.harmonise_pair(exp_df, out_df)
    if har.empty:
        return pd.DataFrame()
    out_name = os.path.basename(out_path).replace("finngen_R12_", "").replace(".gz", "")
    n_out = sample_sizes.get(out_name, 300000)
    har = r1_utils.steiger_filter(har, n_exp=n_exp, n_out=n_out)
    return har[har["steiger_pass"]]


# ---------------------------------------------------------------------------
# FinnGen instrument selection
# ---------------------------------------------------------------------------
def select_finngen_instruments(outcome: str, p_thresh: float = 5e-8) -> pd.DataFrame:
    """Genome-wide significant SNPs for a FinnGen outcome (exposure frame)."""
    return retry_call(lambda: _select_finngen_instruments_impl(outcome, p_thresh),
                      label=f"select_{outcome}")


def _select_finngen_instruments_impl(outcome: str, p_thresh: float) -> pd.DataFrame:
    path = FINNGEN_OUTCOMES[outcome]
    usecols = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
               "pval", "beta", "sebeta", "af_alt"]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=500000,
                             usecols=lambda c: c in usecols):
        chunk = chunk.rename(columns={"#chrom": "chrom", "rsids": "snp",
                                      "pval": "p", "sebeta": "se"})
        chunk = chunk.dropna(subset=["snp", "beta", "se", "p"])
        chunk = chunk[chunk["p"] < p_thresh]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    # Some FinnGen rows carry comma-separated rsid aliases; keep the first alias.
    df["snp"] = df["snp"].astype(str).str.split(",").str[0]
    df["chrom"] = df["chrom"].astype(str)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    df["ea"] = df["alt"]
    df["oa"] = df["ref"]
    df["af_alt"] = pd.to_numeric(df["af_alt"], errors="coerce")
    df = df.drop_duplicates("snp", keep="first")
    return df.sort_values("p").reset_index(drop=True)


def ld_clump_windowed(df: pd.DataFrame, r2_thresh: float = 0.001,
                      kb: int = 10000) -> pd.DataFrame:
    """Memory-safe greedy LD clumping for very large candidate sets.

    Same semantics as r1_utils.ld_clump (sort by p; keep a SNP only if no already
    kept SNP within the kb window has r2 > thresh on the 1KG EUR panel) but LD is
    computed on demand per (candidate, kept) pair inside the window instead of
    materialising the full N x N r2 matrix (which is infeasible when a single
    chromosome holds >20k candidates, e.g. the MHC region for blood traits).
    """
    if df.empty or len(df) == 1:
        return df.copy()
    df = df.copy().sort_values("p").reset_index(drop=True)
    df["snp"] = df["snp"].astype(str)
    lookup = r1_utils.build_snp_lookup(df["snp"].tolist())
    rows = []
    for i, row in df.iterrows():
        info = lookup.get(row["snp"])
        if info is not None:
            rows.append((i, info[0], info[1], info[2]))
    if not rows:
        return pd.DataFrame()
    df2 = df.loc[[r[0] for r in rows]].copy().reset_index(drop=True)
    df2["bed_idx"] = [r[1] for r in rows]
    df2["chrom_bim"] = [r[2] for r in rows]
    df2["pos_bim"] = [r[3] for r in rows]
    kept_snps = []
    win = kb * 1000
    for chrom, sub in df2.groupby("chrom_bim"):
        sub = sub.reset_index(drop=True)
        G = r1_utils.read_bed_indices(sub["bed_idx"].tolist())
        G = r1_utils.impute_mean(G).astype(np.float64)
        mean = G.mean(axis=0)
        std = G.std(axis=0)
        std[std == 0] = 1.0
        Z = (G - mean) / std
        n_ind = Z.shape[0]
        kept = []
        kept_pos = []
        for i in range(len(sub)):
            pos_i = sub.loc[i, "pos_bim"]
            ok = True
            for kk, k in enumerate(kept):
                if abs(pos_i - kept_pos[kk]) <= win:
                    r = float(np.dot(Z[:, i], Z[:, k]) / n_ind)
                    if r * r > r2_thresh:
                        ok = False
                        break
            if ok:
                kept.append(i)
                kept_pos.append(pos_i)
        kept_snps.extend(sub.loc[kept, "snp"].tolist())
    return df[df["snp"].isin(set(kept_snps))].reset_index(drop=True)


def clump_instruments(df: pd.DataFrame, r2: float = 0.001, kb: int = 10000) -> pd.DataFrame:
    """LD-clump (greedy by p, r2 < thresh within kb window, 1KG EUR panel).

    Uses ld_clump_windowed, a memory-safe implementation of the r1_utils.ld_clump
    algorithm that (a) only compares SNPs on the same chromosome (r1_utils.ld_clump
    lacks a chromosome check and spuriously drops SNPs on different chromosomes
    with coincidentally close 1KG positions, because r2 sampling noise at n=503
    exceeds r2=0.001) and (b) computes LD on demand within the window instead of
    building the full N x N matrix (OOM for >20k candidates on one chromosome,
    e.g. MHC for blood traits).  Verified to reproduce r1_utils.ld_clump exactly
    once the same-chromosome check is added.  Original (GRCh38) coordinates are
    preserved in the output.
    """
    if df.empty:
        return df
    return ld_clump_windowed(df, r2_thresh=r2, kb=kb)


# ---------------------------------------------------------------------------
# Pan-UKBB outcome scanning (parallel over component files)
# ---------------------------------------------------------------------------
def _scan_one_component(args):
    comp, fname, keys = args
    df = r3_utils.load_panukbb(fname, pop="EUR", p_thresh=None, rsids=set(keys))
    if df is None or df.empty:
        return comp, pd.DataFrame()
    df["component"] = comp
    return comp, df


def scan_panukbb_components(keys, components=None, n_jobs: int = 6):
    """Scan component Pan-UKBB files for chr:pos:ref:alt keys. Returns dict comp->df."""
    import multiprocessing as mp
    components = components or list(PANUKBB_COMPONENTS)
    args = [(c, PANUKBB_COMPONENTS[c][0], list(keys)) for c in components]
    if n_jobs > 1 and len(args) > 1:
        with mp.Pool(min(n_jobs, len(args))) as pool:
            results = pool.map(_scan_one_component, args)
    else:
        results = [_scan_one_component(a) for a in args]
    return {comp: df for comp, df in results}


def coord_key(chrom, pos, a, b) -> str:
    return f"{chrom}:{pos}:{a}:{b}"


def panukbb_keys_via_bim(inst_df: pd.DataFrame):
    """Map instrument rsids to 1KG (GRCh37) coordinates and build Pan-UKBB query keys.

    Pan-UKBB summary stats are GRCh37 (verified: rs429358 @ 19:45411941), FinnGen
    R12 is GRCh38, so the two cannot be matched by raw coordinates.  The bridge is
    rsid -> 1KG-bim GRCh37 (chrom, pos) -> Pan-UKBB chr:pos:ref:alt key, using the
    FinnGen allele letters (build-independent).  Returns (keys set, key -> rsid).
    """
    lookup = r1_utils.build_snp_lookup(inst_df["snp"].astype(str).tolist())
    keys = set()
    key2rsid = {}
    for _, r in inst_df.iterrows():
        info = lookup.get(str(r["snp"]))
        if info is None:
            continue
        _, chrom, pos = info
        ref, alt = str(r["ref"]).upper(), str(r["alt"]).upper()
        for k in (coord_key(chrom, pos, ref, alt), coord_key(chrom, pos, alt, ref)):
            keys.add(k)
            key2rsid[k] = r["snp"]
    return keys, key2rsid


def attach_rsid_via_keys(pan_df: pd.DataFrame, key2rsid: dict) -> pd.DataFrame:
    """Annotate scanned Pan-UKBB rows with the instrument rsid via the key map."""
    pan_df = pan_df.copy()
    pan_df["rsid"] = pan_df["snp"].map(key2rsid)
    return pan_df.dropna(subset=["rsid"]).drop_duplicates("rsid", keep="first")


def attach_rsid_by_coord(pan_df: pd.DataFrame, inst_df: pd.DataFrame) -> pd.DataFrame:
    """Map Pan-UKBB key rows to FinnGen rsids by coordinate + allele pair (either order)."""
    m = {}
    for _, r in inst_df.iterrows():
        k1 = coord_key(r["chrom"], r["pos"], r["ref"].upper(), r["alt"].upper())
        k2 = coord_key(r["chrom"], r["pos"], r["alt"].upper(), r["ref"].upper())
        m[k1] = r["snp"]
        m[k2] = r["snp"]
    pan_df = pan_df.copy()
    pan_df["rsid"] = pan_df["snp"].map(m)
    pan_df = pan_df.dropna(subset=["rsid"]).drop_duplicates("rsid", keep="first")
    return pan_df


def panukbb_outcome_frame(pan_df: pd.DataFrame) -> pd.DataFrame:
    """Build a harmonise_pair outcome frame from scanned Pan-UKBB rows.

    The scanned frame carries both 'snp' (the chr:pos:ref:alt query key) and
    'rsid' (instrument rsid from the GRCh37 bridge).  Build the frame
    explicitly from 'rsid' -- a naive rename(rsid->snp) would duplicate the
    'snp' column and silently index on the key string instead of the rsid,
    yielding zero harmonised SNPs.
    """
    out = pd.DataFrame({
        "snp": pan_df["rsid"].astype(str),
        "chrom": pan_df["chrom"].astype(str),
        "pos": pd.to_numeric(pan_df["pos"], errors="coerce"),
        "ref": pan_df["ref"].astype(str).str.upper(),
        "alt": pan_df["alt"].astype(str).str.upper(),
        "beta": pd.to_numeric(pan_df["beta"], errors="coerce"),
        "se": pd.to_numeric(pan_df["se"], errors="coerce"),
        "p": pd.to_numeric(pan_df["p"], errors="coerce"),
        "af_alt": pd.to_numeric(pan_df["af_alt"], errors="coerce"),
    })
    out = out.dropna(subset=["snp", "beta", "se", "af_alt"])
    return out.drop_duplicates("snp", keep="first").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Astle (GRCh37, rsid-based) loaders
# ---------------------------------------------------------------------------
def _astle_standardise(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={"SNP": "snp", "CHR": "chrom", "POS": "pos",
                            "A1": "alt", "A2": "ref", "EAF": "af_alt",
                            "Beta": "beta", "P": "p", "N": "n"})
    df["snp"] = df["snp"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()
    for c in ("beta", "se", "p", "af_alt"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["snp", "beta", "se"])
    return df.drop_duplicates("snp", keep="first").reset_index(drop=True)


def scan_astle_rsids(rsids) -> dict:
    """Subset both Astle files to a set of rsids (effect allele = A1 = alt).

    Uses a streaming zcat|awk pre-filter (near-zero memory, robust to the OOM
    killer under parallel load) with a pandas fallback.
    """
    rsids = set(map(str, rsids))
    import hashlib
    tag = hashlib.md5("\n".join(sorted(rsids)).encode()).hexdigest()[:10]
    out = {}
    os.makedirs(os.path.join(LOGS_DIR, "cache"), exist_ok=True)
    for name, path in ASTLE_FILES.items():
        cache = os.path.join(LOGS_DIR, "cache", f"astle_subset_{name}_{tag}.tsv")
        df = None
        if not os.path.exists(cache):
            try:
                import subprocess
                import tempfile
                with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
                    tf.write("\n".join(sorted(rsids)))
                    idfile = tf.name
                cmd = (f"zcat '{path}' 2>/dev/null | awk -F'\\t' "
                       f"'NR==FNR{{ids[$1];next}} FNR==1{{print;next}} ($1 in ids)' "
                       f"'{idfile}' - > '{cache}'")
                subprocess.run(["bash", "-c", cmd], check=True, timeout=600)
                os.unlink(idfile)
            except Exception:
                cache = None
        if cache is not None and os.path.exists(cache):
            try:
                df = pd.read_csv(cache, sep="\t")
            except Exception:
                df = None
        if df is None:
            rows = []
            for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=200000):
                chunk = chunk[chunk["SNP"].astype(str).isin(rsids)]
                if not chunk.empty:
                    rows.append(chunk)
            df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        out[name] = _astle_standardise(df) if not df.empty else pd.DataFrame()
    return out


def astle_significant(name: str, p_thresh: float = 5e-8) -> pd.DataFrame:
    """All p < threshold SNPs from one Astle file (as an exposure frame).

    Reads a pre-built streaming-filter cache (results/r5/logs/cache) when present
    to avoid the memory-heavy full-file pandas scan under parallel load.
    """
    cache = os.path.join(LOGS_DIR, "cache", f"astle_sig_{name}.tsv")
    if os.path.exists(cache):
        df = pd.read_csv(cache, sep="\t")
        df = df[pd.to_numeric(df["P"], errors="coerce") < p_thresh]
        if df.empty:
            return pd.DataFrame()
        df = _astle_standardise(df)
        df["ea"] = df["alt"]
        df["oa"] = df["ref"]
        return df.sort_values("p").reset_index(drop=True)
    return retry_call(lambda: _astle_significant_impl(name, p_thresh),
                      label=f"astle_sig_{name}")


def _astle_significant_impl(name: str, p_thresh: float) -> pd.DataFrame:
    path = ASTLE_FILES[name]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=500000):
        chunk = chunk[pd.to_numeric(chunk["P"], errors="coerce") < p_thresh]
        if not chunk.empty:
            rows.append(chunk)
    if not rows:
        return pd.DataFrame()
    df = _astle_standardise(pd.concat(rows, ignore_index=True))
    df["ea"] = df["alt"]
    df["oa"] = df["ref"]
    return df.sort_values("p").reset_index(drop=True)


# ---------------------------------------------------------------------------
# MR batteries
# ---------------------------------------------------------------------------
def mr_battery(har: pd.DataFrame) -> list:
    """IVW fixed/random + Egger + weighted median on a harmonise_pair output."""
    bx = har["beta"].values.astype(float)
    by = har["beta_outcome"].values.astype(float)
    sy = har["se_outcome"].values.astype(float)
    n = len(har)
    rows = []

    def add(method, b, se, p, extra=None):
        row = {"method": method, "n": n, "beta": b, "se": se, "p": p,
               "or": np.exp(b), "or_lower": np.exp(b - 1.96 * se),
               "or_upper": np.exp(b + 1.96 * se)}
        if extra:
            row.update(extra)
        rows.append(row)

    b, se, p = r1_utils.mr_ivw(bx, by, sy, random=False)
    add("IVW_fixed", b, se, p)
    b, se, p = r1_utils.mr_ivw(bx, by, sy, random=True)
    add("IVW_random", b, se, p)
    try:
        slope, se_s, p_s, icpt, p_i = r1_utils.mr_egger(bx, by, sy)
        add("MR_Egger", slope, se_s, p_s,
            {"egger_intercept": icpt, "egger_intercept_p": p_i})
    except Exception as exc:
        rows.append({"method": "MR_Egger", "n": n, "beta": np.nan, "se": np.nan,
                     "p": np.nan, "or": np.nan, "or_lower": np.nan,
                     "or_upper": np.nan, "error": str(exc)})
    try:
        b, se, p = r1_utils.weighted_median(bx, by, sy)
        add("weighted_median", b, se, p)
    except Exception as exc:
        rows.append({"method": "weighted_median", "n": n, "beta": np.nan,
                     "se": np.nan, "p": np.nan, "or": np.nan, "or_lower": np.nan,
                     "or_upper": np.nan, "error": str(exc)})
    return rows


def mr_raps_huber(bx, by, sy, k: float = 1.345, max_iter: int = 200, tol: float = 1e-10):
    """Robust IVW (MR-RAPS style) with Huber loss via IRLS.

    Minimises sum_i rho_k((by_i - b*bx_i)/sy_i) with weights 1/sy_i^2.
    Overdispersion phi = max(1, robust residual chi2 / df); SE inflated by sqrt(phi).
    """
    bx = np.asarray(bx, float)
    by = np.asarray(by, float)
    sy = np.asarray(sy, float)
    w = 1.0 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)  # IVW start
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        r = (by - b * bx) / sy
        u = np.where(np.abs(r) <= k, 1.0, k / np.abs(r))
        wu = w * u
        denom = np.sum(wu * bx ** 2)
        if denom <= 0:
            break
        b_new = np.sum(wu * bx * by) / denom
        if abs(b_new - b) < tol:
            b = b_new
            break
        b = b_new
    r = (by - b * bx) / sy
    u = np.where(np.abs(r) <= k, 1.0, k / np.abs(r))
    wu = w * u
    df = max(len(bx) - 1, 1)
    # robust residual chi-square on the standardised scale (r already includes 1/sy)
    chi2 = float(np.sum(u * r ** 2))
    phi = max(1.0, chi2 / df)
    se = float(np.sqrt(phi / np.sum(wu * bx ** 2)))
    p = 2 * (1 - stats.norm.cdf(abs(b / se))) if se > 0 else np.nan
    return float(b), se, float(p), phi, n_iter
