# -*- coding: utf-8 -*-
"""Forest plot for the NHANES target-trial emulation (ST146)."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"
FIG = r"D:\衰老研究\v3_pipeline\Nature子刊\figures"

d = pd.read_csv(OUTT + r"\ST146_nhanes_target_trial.csv")
labels = ["High vs low index (IPTW)", "Per-SD index (IPTW, sensitivity)"]

fig, ax = plt.subplots(figsize=(6.2, 2.6))
y = [0, 1]
for yi, (_, r) in zip(y, d.iterrows()):
    ax.plot([r["ci_lo"], r["ci_hi"]], [yi, yi], color="#2166ac", lw=2)
    ax.plot(r["OR"], yi, "s", color="#b2182b", ms=8)
    ax.text(3.4, yi, "OR %.2f (%.2f–%.2f), P = %.1e" % (r["OR"], r["ci_lo"], r["ci_hi"], r["p"]),
            va="center", fontsize=9)
ax.axvline(1, color="grey", ls="--", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xscale("log")
ax.set_xlim(0.8, 3.2)
ax.set_xticks([0.8, 1.0, 1.5, 2.0, 3.0])
ax.set_xticklabels(["0.8", "1.0", "1.5", "2.0", "3.0"])
ax.set_xlabel("Odds ratio for the allergic composite (log scale)")
ax.set_title("NHANES 2017-18 target-trial emulation (n = 4,403 adults; 689 events;\nstabilized IPTW, income-adjusted propensity, survey-weighted, PSU-clustered)", fontsize=9)
fig.tight_layout()
fig.savefig(FIG + r"\FigD_target_trial_forest.png", dpi=300)
fig.savefig(FIG + r"\FigD_target_trial_forest.pdf")
plt.close(fig)
print("FigD written")
