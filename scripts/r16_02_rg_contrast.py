#!/usr/bin/env python3
"""R16-02 (editorial M1): formal contrast of frailty vs epigenetic-clock
genetic correlations with allergic outcomes.

For each allergic outcome, tests rg(Frailty, outcome) - rg(clock, outcome)
against zero: z = diff / sqrt(se_f^2 + se_c^2 - 2*rho*se_f*se_c).

The two rg estimates share the same outcome GWAS, so their sampling errors
are positively correlated; the primary analysis assumes rho = 0, which
overestimates the SE of the difference and is therefore conservative.
A rho = 0.5 sensitivity (stronger contrasts) is also reported. 12 contrasts
(4 clocks x 3 outcomes), BH-FDR across the family.

Primary input : results/r6/tables/r6_ldsc_rg.csv (in-house LDSC, ST60)
Sensitivity   : Genome Medicine package ST126 (independent Python-3 LDSC
                reimplementation) where the same pairs exist.

Output: results/r16/tables/ST134_frailty_vs_clock_rg_contrast.csv
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\衰老研究\v3_pipeline"
RG = os.path.join(BASE, "results", "r6", "tables", "r6_ldsc_rg.csv")
ST126 = os.path.join(BASE, "Genome Medicine", "04_补充数据_supplementary_data",
                     "ST126_independent_ldsc_rg.csv")
OUTT = os.path.join(BASE, "results", "r16", "tables")
os.makedirs(OUTT, exist_ok=True)

OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
CLOCKS = ["GrimAge", "Hannum", "IEAA", "PhenoAge"]


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o]
    m = len(r)
    adj = np.minimum.accumulate((r * m / (np.arange(m) + 1))[::-1])[::-1]
    tmp = np.empty(m)
    tmp[o] = np.clip(adj, 0, 1)
    q[ok] = tmp
    return q


def contrast(rg_f, se_f, rg_c, se_c, rho=0.0):
    diff = rg_f - rg_c
    se_d = np.sqrt(se_f ** 2 + se_c ** 2 - 2 * rho * se_f * se_c)
    z = diff / se_d
    p = 2 * stats.norm.sf(abs(z))
    return diff, se_d, z, p


def main():
    rg = pd.read_csv(RG)
    rows = []
    for oc in OUTCOMES:
        f = rg[(rg.trait1 == "Frailty") & (rg.trait2 == oc)]
        if f.empty:
            f = rg[(rg.trait2 == "Frailty") & (rg.trait1 == oc)]
        f = f.iloc[0]
        for ck in CLOCKS:
            c = rg[(rg.trait1 == ck) & (rg.trait2 == oc)]
            if c.empty:
                c = rg[(rg.trait2 == ck) & (rg.trait1 == oc)]
            c = c.iloc[0]
            diff, se_d, z, p = contrast(f["rg"], f["rg_se"], c["rg"],
                                        c["rg_se"], rho=0.0)
            _, se_d5, _, p5 = contrast(f["rg"], f["rg_se"], c["rg"],
                                       c["rg_se"], rho=0.5)
            rows.append(dict(
                outcome=oc, clock=ck, source="in_house_ST60",
                rg_frailty=f["rg"], se_frailty=f["rg_se"],
                rg_clock=c["rg"], se_clock=c["rg_se"],
                clock_rg_hi95=c["rg"] + 1.96 * c["rg_se"],
                frailty_rg_lo95=f["rg"] - 1.96 * f["rg_se"],
                ci_nonoverlap=bool((f["rg"] - 1.96 * f["rg_se"])
                                   > (c["rg"] + 1.96 * c["rg_se"])),
                diff=diff, se_diff_rho0=se_d, z_rho0=z, p_rho0=p,
                se_diff_rho0p5=se_d5, p_rho0p5=p5))

    # sensitivity on the independent LDSC implementation (ST126)
    if os.path.exists(ST126):
        s = pd.read_csv(ST126)
        for oc in OUTCOMES:
            f = s[(s.exposure == "Frailty") & (s.outcome == oc)]
            if f.empty:
                continue
            f = f.iloc[0]
            for ck in CLOCKS:
                c = s[(s.exposure == ck) & (s.outcome == oc)]
                if c.empty:
                    continue
                c = c.iloc[0]
                diff, se_d, z, p = contrast(
                    f["rg_independent"], f["se_independent"],
                    c["rg_independent"], c["se_independent"], rho=0.0)
                rows.append(dict(
                    outcome=oc, clock=ck, source="independent_ST126",
                    rg_frailty=f["rg_independent"],
                    se_frailty=f["se_independent"],
                    rg_clock=c["rg_independent"],
                    se_clock=c["se_independent"],
                    clock_rg_hi95=c["rg_independent"]
                    + 1.96 * c["se_independent"],
                    frailty_rg_lo95=f["rg_independent"]
                    - 1.96 * f["se_independent"],
                    ci_nonoverlap=bool(
                        (f["rg_independent"] - 1.96 * f["se_independent"])
                        > (c["rg_independent"] + 1.96 * c["se_independent"])),
                    diff=diff, se_diff_rho0=se_d, z_rho0=z, p_rho0=p,
                    se_diff_rho0p5=np.nan, p_rho0p5=np.nan))

    res = pd.DataFrame(rows)
    for src, sub in res.groupby("source"):
        res.loc[sub.index, "q_bh_rho0"] = bh(sub["p_rho0"].values)
    res.to_csv(os.path.join(
        OUTT, "ST134_frailty_vs_clock_rg_contrast.csv"), index=False)
    pd.set_option("display.width", 200)
    print(res[["source", "outcome", "clock", "rg_frailty", "rg_clock",
               "diff", "se_diff_rho0", "z_rho0", "p_rho0", "q_bh_rho0",
               "ci_nonoverlap"]].to_string(index=False))
    print("\nDONE ->", OUTT)


if __name__ == "__main__":
    main()
