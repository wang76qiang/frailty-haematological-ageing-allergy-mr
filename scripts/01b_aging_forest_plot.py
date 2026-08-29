"""R5-01 森林图：5 衰老工具 × 3 过敏结局的 IVW OR(95% CI)，按路径 C 着色。"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
df = pd.read_csv(os.path.join(BASE, "results", "r5", "tables", "r5_orthogonal_aging_mr.csv"))
os.makedirs(os.path.join(BASE, "results", "r5", "figures"), exist_ok=True)

traits = ["Frailty", "PhenoAge", "IEAA", "Hannum", "GrimAge"]
outcomes = [("ALLERG_ASTHMA", "Asthma"), ("L12_ATOPIC", "Atopic dermatitis"),
            ("ALLERG_RHINITIS", "Allergic rhinitis")]

fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), sharey=True)
for ax, (oc, label) in zip(axes, outcomes):
    for i, tr in enumerate(traits):
        row = df[(df["trait"] == tr) & (df["outcome"] == oc)]
        if row.empty or pd.isna(row["ivw_beta"].iloc[0]):
            continue
        b, se = row["ivw_beta"].iloc[0], row["ivw_se"].iloc[0]
        q = row["ivw_q"].iloc[0] if "ivw_q" in row else np.nan
        lo, hi = np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
        or_ = np.exp(b)
        sig = (not pd.isna(q)) and q < 0.05
        color = ("#c0392b" if b > 0 else "#2471a3") if sig else "#9aa0a6"
        ax.errorbar(or_, i, xerr=[[or_ - lo], [hi - or_]], fmt="o", color=color,
                    ecolor=color, capsize=3, markersize=7 if sig else 5,
                    zorder=3)
    ax.axvline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlim(0.7, 6)
    ax.set_xticks([0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
    ax.set_xticklabels(["0.8", "1.0", "1.5", "2.0", "3.0", "4.0", "6.0"], fontsize=8)
    ax.set_title(label, fontsize=12)
    ax.set_xlabel("IVW OR per 1-SD genetic liability (95% CI)", fontsize=9)
    ax.grid(axis="x", alpha=0.25, zorder=0)

axes[0].set_yticks(range(len(traits)))
axes[0].set_yticklabels(["Frailty index", "PhenoAge accel", "IEAA", "Hannum accel", "GrimAge accel"],
                        fontsize=10)
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=8, label="Risk-increasing (FDR<0.05)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#2471a3", markersize=8, label="Protective (FDR<0.05)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#9aa0a6", markersize=6, label="Not significant")]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.suptitle("R5-01  Orthogonal aging instruments vs allergic disease (FinnGen R12) — decision path C",
             fontsize=12, y=0.99)
fig.tight_layout(rect=[0, 0.07, 1, 0.96])
out = os.path.join(BASE, "results", "r5", "figures", "r5_01_orthogonal_aging_forest.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print("saved", out)
