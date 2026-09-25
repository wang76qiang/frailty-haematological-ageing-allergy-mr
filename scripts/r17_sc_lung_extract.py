# -*- coding: utf-8 -*-
"""Lung.h5ad: same extraction for IL18 / TSLP / IL33 / TNFSF14 / IL4 (cell-source
localization in the barrier tissue; donors 33-69y)."""
import h5py
import numpy as np
import pandas as pd
import time

H5 = r"G:\衰老\Lung.h5ad"
OUT = r"D:\衰老研究\v3_pipeline\Nature子刊\data"

TARGETS = {
    "ENSG00000113520": "IL4",
    "ENSG00000150782": "IL18",
    "ENSG00000125735": "TNFSF14",
    "ENSG00000137033": "IL33",
    "ENSG00000162591": "TSLP",
    "ENSG00000169194": "IL13",
    "ENSG00000113525": "IL5",
}

f = h5py.File(H5, "r")
n_cells = f["obs"]["_index"].shape[0]
print("cells:", n_cells)

# obs
def read_cat(grp):
    cats = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["categories"][:]]
    codes = grp["codes"][:]
    return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)

cols = {}
for k in ["donor_id", "cell_type", "development_stage", "tissue", "broad_cell_class"]:
    if k in f["obs"]:
        cols[k] = read_cat(f["obs"][k]) if isinstance(f["obs"][k], h5py.Group) else f["obs"][k][:]
obs = pd.DataFrame(cols, index=[x.decode() if isinstance(x, bytes) else str(x) for x in f["obs"]["_index"][:]])
obs.to_pickle(OUT + r"\lung_obs.pkl")
print("obs cols:", list(obs.columns), "| cell types:", obs["cell_type"].nunique() if "cell_type" in obs else "n/a")
print(obs["development_stage"].value_counts() if "development_stage" in obs else obs["donor_id"].value_counts().head(10))

var_index = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["_index"][:]]
col_of = {var_index.index(e): n for e, n in TARGETS.items() if e in var_index}
print("mapped:", sorted(col_of.values()))

# X location
if isinstance(f["X"], h5py.Group):
    X = f["X"]
else:
    X = f  # dense
data, indices, indptr = X["data"], X["indices"], X["indptr"]
ip = indptr[:]
targets_arr = np.array(sorted(col_of.keys()))
tname = {c: col_of[c] for c in targets_arr}
BLOCK = 20000
records = []
totals = np.zeros(n_cells, dtype=np.float64)
t0 = time.time()
for start in range(0, n_cells, BLOCK):
    end = min(start + BLOCK, n_cells)
    s, e = int(ip[start]), int(ip[end])
    d = data[s:e]; ix = indices[s:e]
    cnt = np.diff(ip[start:end + 1])
    totals[start:end] = cnt
    mask = np.isin(ix, targets_arr)
    if mask.any():
        rows = np.repeat(np.arange(end - start), cnt)
        sub_r, sub_c, sub_v = rows[mask], ix[mask], d[mask]
        for c in np.unique(sub_c):
            m2 = sub_c == c
            records.append((sub_r[m2] + start, np.full(m2.sum(), tname[int(c)]), sub_v[m2]))
f.close()

cell_idx = np.concatenate([r[0] for r in records])
df = pd.DataFrame({"cell": obs.index.values[cell_idx],
                   "gene": np.concatenate([r[1] for r in records]),
                   "count": np.concatenate([r[2] for r in records]).astype(np.float64)})
df.to_pickle(OUT + r"\lung_target_counts.pkl")
np.save(OUT + r"\lung_total_counts.npy", totals)
print("lung extraction done, %.0fs; cells with target:", time.time() - t0, df["cell"].nunique())
