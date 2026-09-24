#!/usr/bin/env python3
"""R16-04 (editorial M3): independent-cohort replication of the genetic-level
eosinophil decomposition (ST133) in Pan-UK Biobank.

ST133 showed, in FinnGen R12, that the composite type-2 index's asthma and
rhinitis associations are carried by its eosinophil-percentage component while
the atopic-dermatitis association survives eosinophil removal.  That analysis
used a single outcome cohort (FinnGen).  Here the identical three instruments
(full / eosinophil-free / eosinophil-only) are tested against independent
Pan-UK Biobank European endpoints:
    asthma                 icd10-J45
    allergic rhinitis      phecode-476
    urticaria (control)    icd10-L50

Exposure : results/r3/tables/r3_index_instruments.csv (per-SNP index effect
           beta_I / se_I; component identifies eosinophil_pct SNPs).
Outcomes : data/real/panukbb/*.tsv.bgz, keyed by GRCh37 chr:pos, effect
           allele = alt.
Harmonisation: allele-aware (ea == alt keeps sign; ea == ref flips); palindromic
           SNPs with control ALT frequency in 0.42-0.58 dropped (mirrors ST133).
Estimators: IVW (fixed + random), weighted median, MR-Egger (+ intercept).

Output: results/r16/tables/ST136_eosinophil_free_panukbb_mr.csv
"""
import gzip
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\衰老研究\v3_pipeline"
PAN = os.path.join(BASE, "data", "real", "panukbb")
IDX = os.path.join(BASE, "results", "r3", "tables", "r3_index_instruments.csv")
OUTT = os.path.join(BASE, "results", "r16", "tables")
os.makedirs(OUTT, exist_ok=True)

OUTCOMES = {
    "UKB_asthma": "icd10-J45-both_sexes.tsv.bgz",
    "UKB_allergic_rhinitis": "phecode-476-both_sexes.tsv.bgz",
    "UKB_urticaria_control": "icd10-L50-both_sexes.tsv.bgz",
}


def load_instruments():
    d = pd.read_csv(IDX).drop_duplicates("snp")
    d["key"] = d.chrom.astype(str) + ":" + d.pos.astype(str)
    d["ea"] = d.ea.str.upper()
    d["oa"] = d.oa.str.upper()
    arms = {
        "full_index": d,
        "eosinophil_free_index": d[d["component"] != "eosinophil_pct"],
        "eosinophil_only": d[d["component"] == "eosinophil_pct"],
    }
    out = []
    for name, sub in arms.items():
        s = sub.copy()
        s["instrument"] = name
        out.append(s)
    print("[instr] full={} eos_free={} eos_only={}".format(
        len(arms["full_index"]), len(arms["eosinophil_free_index"]),
        len(arms["eosinophil_only"])))
    return pd.concat(out, ignore_index=True)


def read_outcome(path, keymap):
    rows = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        def pick(cands):
            for c in cands:
                if c in ci:
                    return c
            raise KeyError(list(ci))
        bcol = pick(["beta_meta_hq", "beta_meta", "beta_EUR"])
        scol = pick(["se_meta_hq", "se_meta", "se_EUR"])
        acol = pick(["af_controls_meta_hq", "af_controls_meta", "af_controls_EUR"])
        for line in f:
            p = line.rstrip("\n").split("\t")
            key = p[0] + ":" + p[1]
            if key in keymap and p[ci[bcol]] not in ("NA", "") and p[ci[scol]] not in ("NA", ""):
                try:
                    af = float(p[ci[acol]])
                except (ValueError, IndexError):
                    af = np.nan
                rows.append({"key": key, "ref_o": p[ci["ref"]].upper(),
                             "alt_o": p[ci["alt"]].upper(),
                             "beta_o": float(p[ci[bcol]]),
                             "se_o": float(p[ci[scol]]), "af_alt_o": af})
    return pd.DataFrame(rows)


def estimators(m):
    x, y, sy = m["bx"].values, m["by"].values, m["se_out"].values
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
                egger_int_p=float(2 * stats.norm.sf(abs(i_eg / se_i)))
                if np.isfinite(se_i) and se_i > 0 else np.nan,
                Q=Q, Q_p=float(stats.chi2.sf(Q, df)),
                I2_pct=float(max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0))


def main():
    inst = load_instruments()
    inst["key"] = inst.chrom.astype(str) + ":" + inst.pos.astype(str)
    keymap = dict(zip(inst.key, inst.snp))
    rows = []
    for oc, fn in OUTCOMES.items():
        path = os.path.join(PAN, fn)
        out = read_outcome(path, keymap)
        print(f"[{oc}] matched {len(out)} instrument SNPs")
        for name, sub in inst.groupby("instrument"):
            m = sub.merge(out, on="key", how="inner")
            for c in ("ea", "oa", "alt_o", "ref_o"):
                m[c] = m[c].str.upper()
            pal = (m.ea.isin(["A", "T"]) & m.oa.isin(["A", "T"])) | \
                  (m.ea.isin(["C", "G"]) & m.oa.isin(["C", "G"]))
            same = m.ea == m.alt_o
            flip = m.ea == m.ref_o
            amb = pal & (m.af_alt_o.between(0.42, 0.58))
            m = m[(same | flip) & ~amb].copy()
            if m.empty:
                print(f"  {name}: no SNPs")
                continue
            m["bx"] = m.beta_I
            m["se_x"] = m.se_I
            m["by"] = np.where(m.ea == m.alt_o, m.beta_o, -m.beta_o)
            m["se_out"] = m.se_o
            e = estimators(m)
            if e is None:
                continue
            rows.append(dict(instrument=name, outcome=oc, **e))
            print(f"  {name:22s} k={e['n_snp']:3d} OR={e['ivw_or']:.3f} "
                  f"P={e['ivw_p']:.2e} (RE P={e['ivw_re_p']:.2e}) "
                  f"WM={e['wm_or']:.3f}")
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTT, "ST136_eosinophil_free_panukbb_mr.csv"),
               index=False)
    print("\nDONE ->", OUTT)


if __name__ == "__main__":
    main()
