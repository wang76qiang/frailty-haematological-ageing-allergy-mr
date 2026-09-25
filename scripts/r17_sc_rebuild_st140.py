# -*- coding: utf-8 -*-
"""Rebuild ST140 as the primary Treg age-analysis result table with explicit tiers:
primary   = equal-donor-weighted WLS on donor log expression rates (12 donors)
secondary = count-weighted Poisson GLM (pseudo-replication-sensitive; sensitivity only)
plus Spearman (equal donor weight) and BH-FDR within the pre-declared 5-gene Treg family.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

DATA = r"D:\衰老研究\v3_pipeline\Nature子刊\data"
OUTT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"

obs = pd.read_pickle(DATA + r"\immune_obs.pkl")
counts = pd.read_pickle(DATA + r"\immune_target_counts.pkl")
tot = np.load(DATA + r"\immune_total_counts.npy")
obs["total_counts"] = tot
obs["age"] = obs["development_stage"].str.extract(r"(\d+)-year-old").astype(float)
dage = obs.groupby("donor_id")["age"].first()
long = counts.merge(obs[["donor_id", "cell_type", "age", "total_counts"]],
                    left_on="cell", right_index=True, how="left")

TREG_GENES = ["IL4", "IL13", "GATA3", "FOXP3", "IL5"]
sub = long[long["cell_type"] == "regulatory T cell"]
tot_d = obs[obs["cell_type"] == "regulatory T cell"].groupby("donor_id")["total_counts"].sum()

def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q

rows = []
for g in TREG_GENES:
    gc = sub[sub["gene"] == g].groupby("donor_id")["count"].sum()
    d = pd.DataFrame({"counts": gc}).reindex(tot_d.index, fill_value=0.0)
    d["total"] = tot_d
    d["age"] = dage[d.index]
    d["rate"] = d["counts"] / d["total"]
    d = d[d["total"] > 0]
    # primary: equal-weight WLS on log rate
    lr = np.log(d["rate"] + 0.5 / d["total"])
    X = sm.add_constant(d["age"] - d["age"].mean())
    w = sm.WLS(lr, X, weights=np.ones(len(d))).fit()
    beta, se = w.params.iloc[1], w.bse.iloc[1]
    # secondary: pooled Poisson
    Xp = sm.add_constant(d["age"] - d["age"].mean())
    po = sm.GLM(d["counts"], Xp, family=sm.families.Poisson(), offset=np.log(d["total"])).fit()
    rho, sp = stats.spearmanr(d["age"], d["rate"])
    rows.append(dict(gene=g, n_donors=len(d), total_gene_counts=float(d["counts"].sum()),
                     primary_rr_per_decade=np.exp(beta * 10),
                     primary_ci_lo=np.exp((beta - 1.96 * se) * 10),
                     primary_ci_hi=np.exp((beta + 1.96 * se) * 10),
                     primary_p=w.pvalues.iloc[1],
                     spearman_rho=rho, spearman_p=sp,
                     pooled_poisson_rr_per_decade=np.exp(po.params.iloc[1] * 10),
                     pooled_poisson_p=po.pvalues.iloc[1]))
res = pd.DataFrame(rows)
res["q_primary_family"] = bh(res["primary_p"].values)
res = res[["gene", "n_donors", "total_gene_counts", "primary_rr_per_decade", "primary_ci_lo", "primary_ci_hi",
           "primary_p", "q_primary_family", "spearman_rho", "spearman_p",
           "pooled_poisson_rr_per_decade", "pooled_poisson_p"]]
res.to_csv(OUTT + r"\ST140_treg_age_analysis.csv", index=False)
pd.set_option("display.width", 250)
print(res.to_string())
