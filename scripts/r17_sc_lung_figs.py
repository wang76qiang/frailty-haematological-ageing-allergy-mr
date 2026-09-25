# -*- coding: utf-8 -*-
"""Lung cell-source table + all figures for the Nature-tier strengthening package."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = r"D:\衰老研究\v3_pipeline\Nature子刊\data"
OUTT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"
FIG = r"D:\衰老研究\v3_pipeline\Nature子刊\figures"

obs = pd.read_pickle(DATA + r"\lung_obs.pkl")
counts = pd.read_pickle(DATA + r"\lung_target_counts.pkl")
long = counts.merge(obs[["donor_id", "cell_type", "development_stage"]],
                    left_on="cell", right_index=True, how="left")

rows = []
for g in ["IL18", "TSLP", "IL33", "TNFSF14", "IL4"]:
    s = long[long["gene"] == g]
    npos = s[s["count"] > 0].groupby("cell_type")["cell"].nunique()
    ntot = obs.groupby("cell_type").size()
    gc = s.groupby("cell_type")["count"].sum()
    for ct in ntot.index:
        rows.append(dict(gene=g, cell_type=ct, n_cells=int(ntot[ct]),
                         n_expressing=int(npos.get(ct, 0)),
                         pct_expressing=100 * npos.get(ct, 0) / ntot[ct],
                         total_counts=float(gc.get(ct, 0.0))))
src = pd.DataFrame(rows).sort_values(["gene", "pct_expressing"], ascending=[True, False])
src.to_csv(OUTT + r"\ST145_lung_cell_source_expression.csv", index=False)
pd.set_option("display.width", 250)
for g in ["IL18", "TSLP", "IL33"]:
    print("--- Lung", g)
    print(src[(src["gene"] == g) & (src["n_cells"] > 200)].head(6).to_string())

# ---------------- Figures ----------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})

# Fig A: donor-level Treg rates vs age (IL5 primary, IL4 secondary)
obs_i = pd.read_pickle(DATA + r"\immune_obs.pkl")
counts_i = pd.read_pickle(DATA + r"\immune_target_counts.pkl")
tot_i = np.load(DATA + r"\immune_total_counts.npy")
obs_i["total_counts"] = tot_i
obs_i["age"] = obs_i["development_stage"].str.extract(r"(\d+)-year-old").astype(float)
dage = obs_i.groupby("donor_id")["age"].first()
long_i = counts_i.merge(obs_i[["donor_id", "cell_type", "age", "total_counts"]],
                        left_on="cell", right_index=True, how="left")

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), sharex=True)
for ax, g in zip(axes, ["IL5", "IL4"]):
    sub = long_i[(long_i["cell_type"] == "regulatory T cell") & (long_i["gene"] == g)]
    tot = obs_i[obs_i["cell_type"] == "regulatory T cell"].groupby("donor_id")["total_counts"].sum()
    gc = sub.groupby("donor_id")["count"].sum()
    rate = (gc / tot).dropna()
    age = dage[rate.index]
    ax.scatter(age, 1e4 * rate, s=28, color="#b2182b", edgecolor="k", linewidth=0.4, zorder=3)
    if g == "IL5":
        X = sm_add = np.vstack([np.ones(len(age)), age - age.mean()]).T
        beta = np.linalg.lstsq(X, np.log(rate + 0.5 / tot[rate.index]), rcond=None)[0]
        xs = np.linspace(age.min(), age.max(), 50)
        ax.plot(xs, 1e4 * np.exp(beta[0] + beta[1] * (xs - age.mean())), color="grey", lw=1, ls="--")
    ax.set_xlabel("Donor age (years)")
    ax.set_ylabel("%s expression rate in Treg\n(counts per 10$^4$)" % g)
    n_cells_treg = obs_i[obs_i["cell_type"] == "regulatory T cell"].groupby("donor_id").size()
    note = "n=12 donors, %d cells" % n_cells_treg.sum()
    ax.set_title("%s (%s)" % (g, note), fontsize=9)
fig.tight_layout()
fig.savefig(FIG + r"\FigA_treg_il5_il4_vs_age.png", dpi=300)
fig.savefig(FIG + r"\FigA_treg_il5_il4_vs_age.pdf")
plt.close(fig)

# Fig B: immune cell source (IL18, TNFSF14, IL4, IL5)
imm = pd.read_csv(OUTT + r"\ST141_cell_source_expression.csv")
fig, axes = plt.subplots(1, 4, figsize=(11, 3.2))
for ax, g in zip(axes, ["IL4", "IL5", "IL18", "TNFSF14"]):
    d = imm[(imm["gene"] == g) & (imm["n_cells"] > 500)].head(8).iloc[::-1]
    ax.barh(d["cell_type"].str.wrap(28), d["pct_expressing"], color="#2166ac", edgecolor="k", linewidth=0.3)
    ax.set_title("%s" % g, fontsize=9)
    ax.set_xlabel("% cells expressing")
    ax.tick_params(labelsize=6)
fig.tight_layout()
fig.savefig(FIG + r"\FigB_immune_cell_source.png", dpi=300)
fig.savefig(FIG + r"\FigB_immune_cell_source.pdf")
plt.close(fig)

# Fig C: lung localization (IL18, TSLP, IL33)
fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
for ax, g in zip(axes, ["IL18", "TSLP", "IL33"]):
    d = src[(src["gene"] == g) & (src["n_cells"] > 200)].head(8).iloc[::-1]
    ax.barh(d["cell_type"].str.wrap(28), d["pct_expressing"], color="#b35806", edgecolor="k", linewidth=0.3)
    ax.set_title("Lung: %s" % g, fontsize=9)
    ax.set_xlabel("% cells expressing")
    ax.tick_params(labelsize=6)
fig.tight_layout()
fig.savefig(FIG + r"\FigC_lung_localization.png", dpi=300)
fig.savefig(FIG + r"\FigC_lung_localization.pdf")
plt.close(fig)
print("figures written")
