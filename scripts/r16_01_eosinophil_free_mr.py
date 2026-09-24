#!/usr/bin/env python3
"""R16-01 (editorial M2): eosinophil-free genetic instrument MR.

Rebuilds the composite haematological type-2 index instrument without the
eosinophil-percentage component SNPs and re-runs two-sample MR against the
three FinnGen R12 allergic outcomes, alongside the full instrument and an
eosinophil-only instrument under identical harmonisation (same-pipeline
comparison).

Exposure : results/r3/tables/r3_index_instruments.csv (per-SNP index effect
           beta_I / se_I; component column identifies eosinophil_pct SNPs)
Outcomes : data/real/finngen_full/finngen_R12_{ALLERG_ASTHMA,
           ALLERG_RHINITIS, L12_ATOPIC}.gz (stream-filtered by rsid,
           alt = effect allele; palindromic SNPs with outcome af_alt in
           0.42-0.58 dropped, as in src/r7/20_gwas_phewas_mr.py)
Estimators: IVW (fixed + random), weighted median, MR-Egger (+ intercept)

Output: results/r16/tables/ST133_eosinophil_free_instrument_mr.csv
        (cross-checked against archived r3_index_no_eosinophil_mr.csv)
"""
import gzip
import os
import re
from io import StringIO

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\衰老研究\v3_pipeline"
FIN = os.path.join(BASE, "data", "real", "finngen_full")
IDX = os.path.join(BASE, "results", "r3", "tables", "r3_index_instruments.csv")
ARCH = os.path.join(BASE, "results", "r3", "tables",
                    "r3_index_no_eosinophil_mr.csv")
OUTT = os.path.join(BASE, "results", "r16", "tables")
os.makedirs(OUTT, exist_ok=True)

OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
FG_COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval",
           "mlogp", "beta", "se", "af_alt", "af_alt_cases",
           "af_alt_controls"]


def load_instruments():
    d = pd.read_csv(IDX).drop_duplicates("snp")
    d["ea"] = d["ea"].str.upper()
    d["oa"] = d["oa"].str.upper()
    arms = {
        "full_index": d,
        "eosinophil_free_index": d[d["component"] != "eosinophil_pct"],
        "eosinophil_only": d[d["component"] == "eosinophil_pct"],
    }
    out = []
    for name, sub in arms.items():
        sub = sub.copy()
        sub["instrument"] = name
        out.append(sub)
    inst = pd.concat(out, ignore_index=True)
    inst = inst.rename(columns={"beta_I": "bx", "se_I": "se_x"})
    for name, sub in inst.groupby("instrument"):
        mean_f = float((sub["bx"] / sub["se_x"]).pow(2).mean())
        print(f"[instr] {name}: n={len(sub)}  mean F={mean_f:.1f}")
    return inst


def stream_endpoint(path, rsids):
    rs = set(rsids)
    pat = re.compile("|".join(re.escape(x) for x in sorted(rs)))
    hits = []
    with gzip.open(path, "rt", errors="replace") as fh:
        fh.readline()
        buf = ""
        while True:
            block = fh.read(64 * 1024 * 1024)
            if not block:
                break
            buf += block
            cut = buf.rfind("\n")
            if cut == -1:
                continue
            seg, buf = buf[:cut], buf[cut + 1:]
            if pat.search(seg):
                hits.extend(ln for ln in seg.split("\n") if pat.search(ln))
        if buf and pat.search(buf):
            hits.append(buf)
    if not hits:
        return pd.DataFrame()
    raw = pd.read_csv(StringIO("\n".join(hits)), sep="\t", header=None,
                      names=FG_COLS, low_memory=False)
    tok = raw["rsids"].astype(str).str.split(",")
    exact = tok.apply(lambda ts: any(t in rs for t in ts))
    raw = raw[exact].copy()
    raw["snp"] = tok[exact].apply(lambda ts: next(t for t in ts if t in rs))
    raw = raw.rename(columns={"beta": "beta_out", "se": "se_out"})
    return raw.drop_duplicates("snp")


def harmonise(sub, out):
    # instrument ref/alt columns collide with outcome ref/alt on merge;
    # only ea/oa are needed from the instrument side
    sub = sub.drop(columns=["ref", "alt", "af_alt"], errors="ignore")
    m = sub.merge(out, on="snp", how="inner")
    if m.empty:
        return m
    for c in ("ea", "oa", "alt", "ref"):
        m[c] = m[c].str.upper()
    pal = (m["ea"].isin(["A", "T"]) & m["oa"].isin(["A", "T"])) | \
          (m["ea"].isin(["C", "G"]) & m["oa"].isin(["C", "G"]))
    same = m["ea"] == m["alt"]
    flip = m["ea"] == m["ref"]
    amb = pal & (m["af_alt"].between(0.42, 0.58))
    m = m[(same | flip) & ~amb].copy()
    m["by"] = np.where(m["ea"] == m["alt"], m["beta_out"], -m["beta_out"])
    return m


def estimators(m):
    x, y = m["bx"].values, m["by"].values
    sy = m["se_out"].values
    w = 1.0 / sy ** 2
    denom = np.sum(w * x ** 2)
    if denom <= 0 or len(m) < 3:
        return None
    b = np.sum(w * x * y) / denom
    se = np.sqrt(1.0 / denom)
    Q = float(np.sum(w * (y - b * x) ** 2))
    df = len(m) - 1
    c = Q / df if Q > df else 1.0
    se_re = se * np.sqrt(c)
    X = np.column_stack([np.ones(len(x)), x])
    W = np.diag(w)
    try:
        cov = np.linalg.inv(X.T @ W @ X)
        coef = cov @ (X.T @ W @ y)
        b_eg, i_eg = coef[1], coef[0]
        se_eg, se_i = np.sqrt(cov[1, 1]), np.sqrt(cov[0, 0])
    except np.linalg.LinAlgError:
        b_eg = se_eg = i_eg = se_i = np.nan
    ratio = y / x
    wr = (x ** 2) / (sy ** 2)
    o = np.argsort(ratio)
    cw = np.cumsum(wr[o]) / np.sum(wr)
    b_wm = float(ratio[o][min(int(np.searchsorted(cw, 0.5)), len(x) - 1)])
    return dict(n_snp=len(m),
                mean_F=float((m["bx"] / m["se_x"]).pow(2).mean()),
                ivw_beta=float(b), ivw_se=float(se),
                ivw_or=float(np.exp(b)),
                or_lo=float(np.exp(b - 1.96 * se_re)),
                or_hi=float(np.exp(b + 1.96 * se_re)),
                ivw_p=float(2 * stats.norm.sf(abs(b / se))),
                ivw_re_p=float(2 * stats.norm.sf(abs(b / se_re))),
                wm_or=float(np.exp(b_wm)),
                egger_or=float(np.exp(b_eg)) if np.isfinite(b_eg) else np.nan,
                egger_int=float(i_eg) if np.isfinite(i_eg) else np.nan,
                egger_int_p=float(2 * stats.norm.sf(abs(i_eg / se_i)))
                if np.isfinite(se_i) and se_i > 0 else np.nan,
                Q=Q, Q_p=float(stats.chi2.sf(Q, df)),
                I2_pct=float(max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0))


def main():
    inst = load_instruments()
    rsids = inst["snp"].unique().tolist()
    rows = []
    for oc in OUTCOMES:
        path = os.path.join(FIN, f"finngen_R12_{oc}.gz")
        out = stream_endpoint(path, rsids)
        print(f"[{oc}] matched {len(out)} SNPs")
        for name, sub in inst.groupby("instrument"):
            h = harmonise(sub, out)
            e = estimators(h)
            if e is None:
                print(f"  {name}: insufficient SNPs")
                continue
            rows.append(dict(instrument=name, outcome=oc, **e))
            print(f"  {name}: k={e['n_snp']} OR={e['ivw_or']:.3f} "
                  f"P={e['ivw_p']:.2e} (RE P={e['ivw_re_p']:.2e}) "
                  f"WM OR={e['wm_or']:.3f}")
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTT, "ST133_eosinophil_free_instrument_mr.csv"),
               index=False)

    # cross-check against archived r3 sensitivity (IVW_fixed rows)
    if os.path.exists(ARCH):
        a = pd.read_csv(ARCH)
        a = a[a["method"] == "IVW_fixed"]
        print("\n[cross-check vs r3 archived IVW_fixed]")
        for _, r in a.iterrows():
            cur = res[(res.instrument.str.replace("eosinophil_free_index",
                                                  "no_eosinophil_index")
                       == r["index"]) & (res.outcome == r["outcome"])]
            if len(cur):
                c = cur.iloc[0]
                print(f"  {r['index']:22s} {r['outcome']:16s} "
                      f"r3: n={r['n_snps']} OR={r['or']:.3f} | "
                      f"r16: n={c['n_snp']} OR={c['ivw_or']:.3f}")
    print("\nDONE ->", OUTT)


if __name__ == "__main__":
    main()
