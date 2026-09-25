# -*- coding: utf-8 -*-
"""Stream-extract per-cell counts for target genes from Immune.h5ad
layers/decontXcounts (CSR, int32) + per-cell total counts. Row-block wise."""
import h5py
import numpy as np
import pandas as pd
import time

H5 = r"G:\衰老\Immune.h5ad"
OBS = r"D:\衰老研究\v3_pipeline\Nature子刊\data\immune_obs.pkl"
OUT = r"D:\衰老研究\v3_pipeline\Nature子刊\data"

TARGETS = {
    "ENSG00000113520": "IL4",
    "ENSG00000169194": "IL13",
    "ENSG00000107485": "GATA3",
    "ENSG00000049768": "FOXP3",
    "ENSG00000113525": "IL5",
    "ENSG00000150782": "IL18",
    "ENSG00000125735": "TNFSF14",
    "ENSG00000137033": "IL33",
    "ENSG00000162591": "TSLP",
    "ENSG00000145625": "IL6",
    "ENSG00000232810": "TNF",
    "ENSG00000107477": "RORC",
    "ENSG00000129559": "ACTB",
    "ENSG00000075624": "GAPDH",
}

var_index = pd.read_csv(OUT + r"\immune_var_index.csv", header=None)[0].tolist()
col_of = {}
for ens, name in TARGETS.items():
    if ens in var_index:
        col_of[var_index.index(ens)] = name
    else:
        print("MISSING in var:", ens, name)
print("mapped genes:", sorted(col_of.values()))

obs = pd.read_pickle(OBS)
assert len(obs) == 592317

t0 = time.time()
f = h5py.File(H5, "r")
X = f["layers"]["decontXcounts"]
data, indices, indptr = X["data"], X["indices"], X["indptr"]
n_cells = indptr.shape[0] - 1
ip = indptr[:]  # 592318 int64

targets_arr = np.array(sorted(col_of.keys()))
tname = {c: col_of[c] for c in targets_arr}

BLOCK = 20000
records = []  # (cell_idx, gene_name, count)
totals = np.zeros(n_cells, dtype=np.float64)

for start in range(0, n_cells, BLOCK):
    end = min(start + BLOCK, n_cells)
    s, e = int(ip[start]), int(ip[end])
    d = data[s:e]
    ix = indices[s:e]
    # per-cell totals
    cnt = np.diff(ip[start:end + 1])
    totals[start:end] = cnt
    # target hits
    mask = np.isin(ix, targets_arr)
    if mask.any():
        rows = np.repeat(np.arange(end - start), cnt)
        sub_r = rows[mask]
        sub_c = ix[mask]
        sub_v = d[mask]
        for c in np.unique(sub_c):
            m2 = sub_c == c
            records.append((sub_r[m2] + start, np.full(m2.sum(), tname[int(c)]), sub_v[m2]))
    if (start // BLOCK) % 15 == 0:
        print("rows %d/%d  elapsed %.0fs" % (end, n_cells, time.time() - t0))

f.close()

cell_idx = np.concatenate([r[0] for r in records])
gene = np.concatenate([r[1] for r in records])
val = np.concatenate([r[2] for r in records]).astype(np.float64)

df = pd.DataFrame({"cell": obs.index.values[cell_idx], "gene": gene, "count": val})
df.to_pickle(OUT + r"\immune_target_counts.pkl")
np.save(OUT + r"\immune_total_counts.npy", totals)
print("cells with any target:", df['cell'].nunique(), "| total nnz:", len(df))
print("elapsed %.0fs" % (time.time() - t0))
