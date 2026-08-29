"""Standalone rebuild of r5_bidirectional_matrix.csv + heatmap from the GOLDEN
(high-n, run3) r5_reverse_mr_components.csv.  Does NOT re-run any MR; pure
table/figure regeneration.  Mirrors build_matrix/make_figure in
02_bidirectional_mr.py.
"""
import os
import numpy as np
import pandas as pd

R5 = r"D:\衰老研究\v3_pipeline\results\r5"
TABLES = os.path.join(R5, "tables")
FIGS = os.path.join(R5, "figures")
OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]

comp = pd.read_csv(os.path.join(TABLES, "r5_reverse_mr_components.csv"))
print("components rows:", len(comp), "| asthma wbc n_snps:",
      comp.query("method=='IVW_fixed' and outcome=='ALLERG_ASTHMA' and trait=='wbc'")["n_snps"].tolist())

rows = []
r4_path = os.path.abspath(os.path.join(R5, "..", "r4", "tables", "r4_mr_diagnostics.csv"))
if os.path.exists(r4_path):
    r4 = pd.read_csv(r4_path).drop_duplicates("outcome").set_index("outcome")
    for outcome in OUTCOMES:
        if outcome not in r4.index:
            continue
        r = r4.loc[outcome]
        orv = r.get("ivw_fixed_or")
        rows.append({"direction": "forward_index_to_allergy", "outcome": outcome,
                     "trait": "immune_index", "method": "IVW_fixed",
                     "n_snps": r.get("n_snps_harmonised"),
                     "beta": np.log(orv) if pd.notna(orv) and orv > 0 else np.nan,
                     "se": np.nan, "p": r.get("ivw_fixed_p"),
                     "or": orv, "or_lower": np.nan, "or_upper": np.nan})
ivw = comp[comp["method"] == "IVW_fixed"]
for _, r in ivw.iterrows():
    rows.append({"direction": "reverse_allergy_to_component",
                 "outcome": r["outcome"], "trait": r["trait"],
                 "method": "IVW_fixed", "n_snps": r["n_snps"],
                 "beta": r["beta"], "se": r["se"], "p": r["p"],
                 "or": r["or"], "or_lower": r["or_lower"], "or_upper": r["or_upper"]})
matrix = pd.DataFrame(rows)
p2 = os.path.join(TABLES, "r5_bidirectional_matrix.csv")
matrix.to_csv(p2, index=False)
print(f"Saved matrix: {p2} ({len(matrix)} rows)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
rev = matrix[matrix["direction"] == "reverse_allergy_to_component"].copy()
piv = rev.pivot_table(index="trait", columns="outcome", values="beta")
sig = rev.pivot_table(index="trait", columns="outcome", values="p")
annot = piv.copy().astype(str)
for t in piv.index:
    for o in piv.columns:
        b = piv.loc[t, o]
        p = sig.loc[t, o]
        annot.loc[t, o] = "" if pd.isna(b) else f"{b:+.3f}{'*' if p < 0.05 else ''}"
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(piv, annot=annot, fmt="", cmap="RdBu_r", center=0, ax=ax,
            cbar_kws={"label": "reverse-MR beta (allergy -> component)"})
ax.set_title("R5-02 Reverse MR: genetically predicted allergy -> blood components\n"
             "(* p<0.05, IVW fixed; Astle = independent replication)")
fig.tight_layout()
out = os.path.join(FIGS, "r5_bidirectional_diagram.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved figure:", out)
