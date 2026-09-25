# -*- coding: utf-8 -*-
"""Extract obs metadata (small) + var gene index from Immune.h5ad; save locally."""
import h5py
import numpy as np
import pandas as pd
import io

p = r"G:\衰老\Immune.h5ad"
out = r"D:\衰老研究\v3_pipeline\Nature子刊\data"

f = h5py.File(p, "r")
obs = f["obs"]

def read_cat(grp):
    cats = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["categories"][:]]
    codes = grp["codes"][:]
    return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)

cols = {}
for k in ["donor_id", "cell_type", "development_stage", "broad_cell_class",
          "donor_tissue", "anatomical_position", "sex", "ethnicity_original",
          "n_genes_by_counts", "assay", "method"]:
    if k in obs:
        if isinstance(obs[k], h5py.Group):
            cols[k] = read_cat(obs[k])
        else:
            cols[k] = obs[k][:]

df = pd.DataFrame(cols, index=[x.decode() if isinstance(x, bytes) else str(x) for x in obs["_index"][:]])
df.to_pickle(out + r"\immune_obs.pkl")
print("obs saved:", df.shape)

var_index = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["_index"][:]]
pd.Series(var_index).to_csv(out + r"\immune_var_index.csv", index=False)
print("var genes:", len(var_index))

# also raw var (if different)
try:
    ridx = [x.decode() if isinstance(x, bytes) else str(x) for x in f["raw"]["var"]["_index"][:]]
    pd.Series(ridx).to_csv(out + r"\immune_raw_var_index.csv", index=False)
    print("raw var genes:", len(ridx))
except Exception as e:
    print("raw var err:", e)

for c in ["cell_type", "development_stage", "donor_id", "donor_tissue", "broad_cell_class"]:
    if c in df.columns:
        print("--- %s (%d):" % (c, df[c].nunique()))
        print(df[c].value_counts().head(50))
f.close()
