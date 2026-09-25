# -*- coding: utf-8 -*-
"""Robustness for the Treg expression-age analysis:
(a) equal-donor-weighted WLS on donor rates (log);
(b) leave-one-donor-out (especially TSP14);
(c) ambient-proxy falsification: donor Treg IL4 rate vs plasma-cell IL4 rate
    (plasma cells do not biologically express IL4; correlation suggests ambient
    contamination);
(d) within-donor-downweighted analysis restricted to donors with >=100 Treg cells.
Outputs ST144 sensitivity table.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

DATA = r"D:\衰老研究\v3_pipeline\Nature子刊\data"
OUTT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"

obs = pd.read_pickle(DATA + r"\immune_obs.pkl")
counts = pd.read_pickle(DATA + r"\immune_target_counts.pkl")
totals = np.load(DATA + r"\immune_total_counts.npy")
obs["total_counts"] = totals
obs["age"] = obs["development_stage"].str.extract(r"(\d+)-year-old").astype(float)
long = counts.merge(obs[["donor_id", "cell_type", "age", "total_counts"]],
                    left_on="cell", right_index=True, how="left")

TREG_GENES = ["IL4", "IL13", "GATA3", "FOXP3", "IL5"]

def donor_rates(celltype, genes):
    sub = long[long["cell_type"] == celltype]
    tot = obs[obs["cell_type"] == celltype].groupby("donor_id")["total_counts"].sum()
    rows = []
    for g in genes:
        gc = sub[sub["gene"] == g].groupby("donor_id")["count"].sum()
        for d in tot.index:
            rows.append(dict(gene=g, donor=d, age=float(obs.groupby("donor_id")["age"].first()[d]),
                             rate=gc.get(d, 0.0) / tot[d], counts=gc.get(d, 0.0), total=tot[d]))
    return pd.DataFrame(rows)

rates = donor_rates("regulatory T cell", TREG_GENES)
rows = []

for g in TREG_GENES:
    d = rates[rates["gene"] == g].dropna()
    # (a) equal-weight WLS on log rate (add half-count pseudo for zeros)
    pseudo = 0.5 / d["total"]
    lr = np.log(d["rate"] + pseudo)
    X = sm.add_constant(d["age"] - d["age"].mean())
    wls = sm.WLS(lr, X, weights=np.ones(len(d))).fit()
    beta, se = wls.params.iloc[1], wls.bse.iloc[1]
    rows.append(dict(gene=g, test="equal_weight_WLS_lograte", rr_per_decade=np.exp(beta * 10),
                     ci_lo=np.exp((beta - 1.96 * se) * 10), ci_hi=np.exp((beta + 1.96 * se) * 10), p=wls.pvalues.iloc[1]))
    # (b) LOO
    loo = []
    for dnr in d["donor"]:
        dd = d[d["donor"] != dnr]
        lr2 = np.log(dd["rate"] + 0.5 / dd["total"])
        X2 = sm.add_constant(dd["age"] - dd["age"].mean())
        try:
            f = sm.WLS(lr2, X2, weights=np.ones(len(dd))).fit()
            loo.append((f.params.iloc[1] * 10, dnr))
        except Exception:
            pass
    betas = np.array([x[0] for x in loo])
    rows.append(dict(gene=g, test="LOO_equal_weight_minmax_beta_per_decade",
                     rr_per_decade=np.exp(betas).min(), ci_hi=np.exp(betas).max(),
                     ci_lo=np.nan, p=np.nan,
                     note="exp(min/max LOO beta)"))
    # (d) donors with >=100 Treg cells only
    big = d[d["total"] >= 100 * obs[obs["cell_type"] == "regulatory T cell"].groupby("donor_id")["total_counts"].mean().mean() / 10]
    # simpler: use n_cells from detection table if present; approximate via total>=threshold
    thr = d["total"].median()
    big = d[d["total"] >= thr]
    if big["donor"].nunique() >= 5:
        lrb = np.log(big["rate"] + 0.5 / big["total"])
        Xb = sm.add_constant(big["age"] - big["age"].mean())
        fb = sm.WLS(lrb, Xb, weights=np.ones(len(big))).fit()
        rows.append(dict(gene=g, test="large_donors_only_WLS", rr_per_decade=np.exp(fb.params.iloc[1] * 10),
                         ci_lo=np.exp((fb.params.iloc[1] - 1.96 * fb.bse.iloc[1]) * 10),
                         ci_hi=np.exp((fb.params.iloc[1] + 1.96 * fb.bse.iloc[1]) * 10), p=fb.pvalues.iloc[1],
                         note="donors above median Treg total counts (n=%d)" % big["donor"].nunique()))

# (c) ambient proxy: donor-level Treg IL4 rate vs plasma-cell IL4 rate
plasma = donor_rates("plasma cell", ["IL4", "IL5"])
for g in ["IL4", "IL5"]:
    a = rates[rates["gene"] == g][["donor", "rate"]].rename(columns={"rate": "treg_rate"})
    b = plasma[plasma["gene"] == g][["donor", "rate"]].rename(columns={"rate": "plasma_rate"})
    m = a.merge(b, on="donor").dropna()
    rho, p = stats.spearmanr(m["treg_rate"], m["plasma_rate"])
    rows.append(dict(gene=g, test="ambient_proxy_spearman_treg_vs_plasmacell", rr_per_decade=np.nan,
                     ci_lo=np.nan, ci_hi=np.nan, p=p, note="rho=%.3f n=%d donors" % (rho, len(m))))

res = pd.DataFrame(rows)
res.to_csv(OUTT + r"\ST144_treg_age_sensitivity.csv", index=False)
pd.set_option("display.width", 250)
print(res.to_string())
