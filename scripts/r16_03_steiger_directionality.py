#!/usr/bin/env python3
"""R16-03 (editorial M3): formal SNP-level Steiger directionality tests for
the forward MR (ageing instruments and composite index -> allergic outcomes).

For each instrument SNP, variance explained is computed on both sides as
R2 = 2 * EAF * (1 - EAF) * beta^2 (orientation-free, so no allele
harmonisation is needed), and the directionality of the SNP is tested by a
Fisher-z contrast of r_x = sqrt(R2_x) versus r_y = sqrt(R2_y):

    z = (atanh(r_x) - atanh(r_y)) / sqrt(1/(N_x-3) + 1/(N_y-3))

(as in TwoSampleMR's directionality_test; treats the two correlations as
independent, which is conservative here because both GWAS are independent
studies). Per-SNP N: exposure N from the instrument table (ST02 n_snp) or,
for the index, effective N derived from se and EAF; outcome effective N
derived per SNP from FinnGen sebeta and af_alt (N = 1/(2*se^2*p*(1-p))).

Exposures: ST02 five ageing instruments (GrimAge, Hannum, IEAA, PhenoAge,
Frailty) + composite type-2 index (r3_index_instruments.csv, beta_I/se_I).
Outcomes : FinnGen R12 ALLERG_ASTHMA / ALLERG_RHINITIS / L12_ATOPIC
           (data/real/finngen_full, streamed by rsid).

Aggregate per pair: n SNPs, proportion with R2_x > R2_y (binomial P vs 0.5),
number Steiger-significant in the correct direction, and the
sum(R2_x)/sum(R2_y) ratio.

Output: results/r16/tables/ST135_steiger_directionality.csv (pair summary)
        results/r16/tables/ST135b_steiger_directionality_persnp.csv
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
ST02 = os.path.join(BASE, "Genome Medicine", "04_补充数据_supplementary_data",
                    "ST02_r5_aging_instruments_clumped.csv")
IDX = os.path.join(BASE, "results", "r3", "tables", "r3_index_instruments.csv")
OUTT = os.path.join(BASE, "results", "r16", "tables")
os.makedirs(OUTT, exist_ok=True)

OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
FG_COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval",
           "mlogp", "beta", "se", "af_alt", "af_alt_cases",
           "af_alt_controls"]


def load_exposures():
    st02 = pd.read_csv(ST02)
    st02 = st02.rename(columns={"beta": "bx", "se": "se_x", "eaf": "eaf_x",
                                "n_snp": "n_x"})
    st02["exposure"] = st02["trait"]
    idx = pd.read_csv(IDX).drop_duplicates("snp")
    idx = idx.rename(columns={"beta_I": "bx", "se_I": "se_x",
                              "af_alt": "eaf_x"})
    # effective exposure N for index SNPs from se and EAF
    idx["n_x"] = 1.0 / (2 * idx["se_x"] ** 2 * idx["eaf_x"] * (1 - idx["eaf_x"]))
    idx["exposure"] = "composite_index"
    cols = ["exposure", "snp", "bx", "se_x", "eaf_x", "n_x"]
    return pd.concat([st02[cols], idx[cols]], ignore_index=True)


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
    raw = pd.read_csv(StringIO("\n".join(hits)), sep="\t", header=None,
                      names=FG_COLS, low_memory=False)
    tok = raw["rsids"].astype(str).str.split(",")
    exact = tok.apply(lambda ts: any(t in rs for t in ts))
    raw = raw[exact].copy()
    raw["snp"] = tok[exact].apply(lambda ts: next(t for t in ts if t in rs))
    return raw.drop_duplicates("snp")


def steiger(bx, eaf_x, n_x, by, se_y, eaf_y):
    r2x = 2 * eaf_x * (1 - eaf_x) * bx ** 2
    r2y = 2 * eaf_y * (1 - eaf_y) * by ** 2
    n_y = 1.0 / (2 * se_y ** 2 * eaf_y * (1 - eaf_y))
    rx = np.sqrt(np.clip(r2x, 0, 0.999))
    ry = np.sqrt(np.clip(r2y, 0, 0.999))
    se_z = np.sqrt(1.0 / (n_x - 3) + 1.0 / (n_y - 3))
    z = (np.arctanh(rx) - np.arctanh(ry)) / se_z
    p = 2 * stats.norm.sf(np.abs(z))
    return r2x, r2y, z, p


def main():
    exp = load_exposures()
    rsids = exp["snp"].unique().tolist()
    # canonical allergy-locus windows documented in the manuscript (GRCh37)
    HLA = (6, 25_000_000, 34_000_000)
    Q31 = (5, 131_660_000, 133_680_000)
    # frailty immune-locus membership established by the R6 pipeline (ST63)
    st63_path = os.path.join(BASE, "Genome Medicine",
                             "04_补充数据_supplementary_data",
                             "ST63_r6_frailty_immune_excluded_snps.csv")
    fr_imm = set(pd.read_csv(st63_path)["snp"]) if os.path.exists(st63_path) else set()
    idx_pos = pd.read_csv(IDX).drop_duplicates("snp")[["snp", "chrom", "pos"]]
    rows, snp_rows = [], []
    for oc in OUTCOMES:
        out = stream_endpoint(os.path.join(FIN, f"finngen_R12_{oc}.gz"), rsids)
        print(f"[{oc}] matched {len(out)} instrument SNPs")
        for name, sub in exp.groupby("exposure"):
            m = sub.merge(out[["snp", "chrom", "pos", "beta", "se",
                                "af_alt"]], on="snp", how="inner").dropna()
            if len(m) < 3:
                continue
            r2x, r2y, z, p = steiger(m["bx"].values, m["eaf_x"].values,
                                     m["n_x"].values, m["beta"].values,
                                     m["se"].values, m["af_alt"].values)
            correct = r2x > r2y
            n = len(m)
            n_cor = int(correct.sum())
            binom_p = float(stats.binomtest(n_cor, n, 0.5).pvalue)
            n_sig_cor = int((correct & (p < 0.05)).sum())
            n_sig_wrong = int((~correct & (p < 0.05)).sum())
            rows.append(dict(
                exposure=name, outcome=oc, n_snp=n,
                median_r2_exposure=float(np.median(r2x)),
                median_r2_outcome=float(np.median(r2y)),
                sumR2_ratio=float(r2x.sum() / r2y.sum()),
                n_correct_direction=n_cor, prop_correct=n_cor / n,
                binom_p=binom_p,
                n_steiger_sig_correct=n_sig_cor,
                n_steiger_sig_wrong=n_sig_wrong))
            ps = pd.DataFrame({
                "exposure": name, "outcome": oc, "snp": m["snp"].values,
                "chrom": m["chrom"].values, "pos": m["pos"].values,
                "r2_exposure": r2x, "r2_outcome": r2y,
                "steiger_z": z, "steiger_p": p,
                "direction": np.where(correct, "correct", "wrong")})
            ps["steiger_sig"] = ps["steiger_p"] < 0.05
            ps["in_HLA"] = ((ps.chrom == HLA[0]) &
                            ps.pos.between(HLA[1], HLA[2]))
            ps["in_5q31"] = ((ps.chrom == Q31[0]) &
                             ps.pos.between(Q31[1], Q31[2]))
            ps["in_immune_locus_r6"] = ps["snp"].isin(fr_imm)
            snp_rows.append(ps)
            print(f"  {name:18s} n={n} correct={n_cor}/{n} "
                  f"sig-correct={n_sig_cor} sig-wrong={n_sig_wrong} "
                  f"sumR2 ratio={r2x.sum()/r2y.sum():.1f}")
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTT, "ST135_steiger_directionality.csv"),
               index=False)
    pd.concat(snp_rows, ignore_index=True).to_csv(
        os.path.join(OUTT, "ST135b_steiger_directionality_persnp.csv"),
        index=False)

    # overlap diagnostic: are Steiger-failing SNPs concentrated in canonical
    # allergy loci (HLA / 5q31) or the R6 immune-locus set?
    allps = pd.concat(snp_rows, ignore_index=True)
    print("\n[overlap] sig-wrong SNPs in canonical loci:")
    for (ex, oc), g in allps.groupby(["exposure", "outcome"]):
        sw = g[(g.direction == "wrong") & g.steiger_sig]
        if len(sw) == 0:
            continue
        flag = sw.in_HLA | sw.in_5q31 | sw.in_immune_locus_r6
        print(f"  {ex:18s} {oc:16s} sig-wrong={len(sw)}, "
              f"in canonical/immune loci: {int(flag.sum())} "
              f"({100*flag.mean():.0f}%)")
    print("\nDONE ->", OUTT)


if __name__ == "__main__":
    main()
