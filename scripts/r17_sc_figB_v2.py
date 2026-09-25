# -*- coding: utf-8 -*-
"""FigB v2: immune cell-source using expression rate (counts per 1e4) instead of
fraction expressing (ambient-robust for lowly-expressed genes)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = r"D:\衰老研究\v3_pipeline\Nature子刊\data"
OUTT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"
FIG = r"D:\衰老研究\v3_pipeline\Nature子刊\figures"

obs = pd.read_pickle(DATA + r"\immune_obs.pkl")
tot = np.load(DATA + r"\immune_total_counts.npy")
obs["total_counts"] = tot
counts = pd.read_pickle(DATA + r"\immune_target_counts.pkl")
long = counts.merge(obs[["donor_id", "cell_type", "total_counts"]], left_on="cell", right_index=True, how="left")

rows = []
for g in ["IL4", "IL5", "IL18", "TNFSF14"]:
    s = long[long["gene"] == g]
    gc = s.groupby("cell_type")["count"].sum()
    tc = obs.groupby("cell_type")["total_counts"].sum()
    n = obs.groupby("cell_type").size()
    for ct in tc.index:
        if n[ct] >= 500:
            rows.append(dict(gene=g, cell_type=ct, n_cells=int(n[ct]),
                             rate_per_1e4=1e4 * gc.get(ct, 0.0) / tc[ct]))
d = pd.DataFrame(rows)
d.to_csv(OUTT + r"\ST141_cell_source_expression_rate.csv", index=False)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig, axes = plt.subplots(1, 4, figsize=(11, 3.2))
for ax, g in zip(axes, ["IL4", "IL5", "IL18", "TNFSF14"]):
    dd = d[d["gene"] == g].sort_values("rate_per_1e4").tail(8)
    ax.barh(dd["cell_type"].str.wrap(28), dd["rate_per_1e4"], color="#2166ac", edgecolor="k", linewidth=0.3)
    ax.set_title("%s (rate per 10$^4$)" % g, fontsize=9)
    ax.set_xlabel("counts per 10$^4$")
    ax.tick_params(labelsize=6)
fig.tight_layout()
fig.savefig(FIG + r"\FigB_immune_cell_source.png", dpi=300)
fig.savefig(FIG + r"\FigB_immune_cell_source.pdf")
plt.close(fig)
print(d[d["gene"].isin(["IL4", "IL5"])].sort_values(["gene", "rate_per_1e4"], ascending=[True, False]).groupby("gene").head(5).to_string())
