# -*- coding: utf-8 -*-
"""LINCS Level 5 alarmin-signature connectivity (P0-4 re-computation, v2).

GSE70138 layout: matrix[signature, gene]; META/ROW = Entrez gene ids,
META/COL = sig ids. Alarmin signature: IL33, TSLP, IL25, IL4, IL5, IL13.
For each drug: signature-level mean z across the 6 genes, aggregated across
its signatures (mean, SE, z), BH-FDR within the 11-drug family.
"""
import h5py
import numpy as np
import pandas as pd
from scipy import stats

GCTX = r"E:\lincs\lincs_level5.gctx"
SIGINFO = r"D:\衰老研究\v3_pipeline\data\real\lincs\GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"
OUT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"

ALARMIN_ENTREZ = {"3565": "IL4", "3567": "IL5", "3596": "IL13",
                  "4057": "IL33", "85480": "TSLP", "64806": "IL25"}
DRUGS = ["quercetin", "tofacitinib", "navitoclax", "ruxolitinib", "azithromycin",
         "dasatinib", "azathioprine", "everolimus", "sirolimus", "temsirolimus",
         "baricitinib"]

f = h5py.File(GCTX, "r")
row_ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["0"]["META"]["ROW"]["id"][:]]
col_ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["0"]["META"]["COL"]["id"][:]]
gene_rows = sorted(row_ids.index(e) for e in ALARMIN_ENTREZ if e in row_ids)
present = [ALARMIN_ENTREZ[row_ids[i]] for i in gene_rows]
print("alarmin genes present:", present)
mat = f["0"]["DATA"]["0"]["matrix"]
print("matrix:", mat.shape)

# read the 6 gene columns for all signatures: matrix[:, gene_rows]
sub = mat[:, gene_rows]
print("sub shape:", sub.shape)
f.close()

sig = pd.read_csv(SIGINFO, sep="\t")
col_index = {c: i for i, c in enumerate(col_ids)}
sig["col"] = sig["sig_id"].map(col_index)
sig = sig[sig["pert_iname"].isin(DRUGS)].dropna(subset=["col"])
sig["sig_mean_z"] = np.array(sub)[sig["col"].astype(int).values].mean(axis=1)

agg = sig.groupby("pert_iname")["sig_mean_z"].agg(n_sig="count", mean_z="mean", std="std").reset_index()
agg["se_z"] = agg["std"] / np.sqrt(agg["n_sig"])
agg["z_score"] = agg["mean_z"] / agg["se_z"]
agg["p"] = 2 * stats.norm.sf(np.abs(agg["z_score"]))
agg["connectivity_score"] = -agg["mean_z"]
agg["signature"] = "alarmin_type2_genes"
agg = agg[["pert_iname", "signature", "n_sig", "mean_z", "se_z", "z_score", "p", "connectivity_score"]]

p = agg["p"].values
order = np.argsort(p)
q = np.empty(len(p))
prev = 1.0
for rank in range(len(p), 0, -1):
    i = order[rank - 1]
    prev = min(prev, p[i] * len(p) / rank)
    q[i] = prev
agg["q_alarmin_family"] = q
agg.to_csv(OUT + r"\ST151_lincs_alarmin_signature_connectivity.csv", index=False)
pd.set_option("display.width", 220)
print(agg.sort_values("p").to_string())
