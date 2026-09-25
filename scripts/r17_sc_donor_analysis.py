# -*- coding: utf-8 -*-
"""P0-1 donor-level expression-vs-age analysis, Immune atlas (Tabula Sapiens 2.0).

Primary family (Treg, donors with >=30 Treg cells):
  IL4 (primary), IL13, GATA3, FOXP3, IL5 (context)
  Model: Poisson GLM gene_counts ~ donor_age, offset log(total_counts), per donor.
  Rate ratio reported per decade.
Secondary family (translational):
  IL18 in dominant myeloid compartment; TNFSF14 in dominant lymphoid compartment.
Descriptive: TSLP/IL33 detection; cell-source table.
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

print("donor -> age:")
dage = obs.groupby("donor_id")["age"].first().astype(int)
print(dage.sort_values().to_dict())

# join counts onto obs
long = counts.merge(obs[["donor_id", "cell_type", "donor_tissue", "age", "total_counts"]],
                    left_on="cell", right_index=True, how="left")

def donor_pseudobulk(celltype, genes):
    sub = long[long["cell_type"] == celltype]
    tot = obs[obs["cell_type"] == celltype].groupby("donor_id")["total_counts"].sum()
    rows = []
    for g in genes:
        gc = sub[sub["gene"] == g].groupby("donor_id")["count"].sum()
        for d in tot.index:
            rows.append(dict(cell_type=celltype, gene=g, donor=d,
                             gene_counts=gc.get(d, 0.0), total_counts=tot[d],
                             age=float(dage[d]),
                             n_cells=int((obs["cell_type"] == celltype).groupby(obs["donor_id"]).size().get(d, 0))))
    return pd.DataFrame(rows)

def poisson_age(df):
    """Poisson GLM counts ~ age with offset log(total). Returns RR per decade (95% CI), P."""
    d = df[(df["total_counts"] > 0)]
    if d["gene_counts"].sum() == 0 or d["donor"].nunique() < 6:
        return dict(rr_per_decade=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan, n_donors=d["donor"].nunique(),
                    total_gene_counts=d["gene_counts"].sum())
    X = sm.add_constant(d["age"] - d["age"].mean())
    m = sm.GLM(d["gene_counts"], X, family=sm.families.Poisson(),
               offset=np.log(d["total_counts"])).fit()
    beta = m.params.iloc[1]
    se = m.bse.iloc[1]
    rr = np.exp(beta * 10)
    return dict(rr_per_decade=rr, ci_lo=np.exp((beta - 1.96 * se) * 10), ci_hi=np.exp((beta + 1.96 * se) * 10),
                p=m.pvalues.iloc[1], n_donors=int(d["donor"].nunique()), total_gene_counts=float(d["gene_counts"].sum()))

# ---------- Treg primary family ----------
treg = obs[obs["cell_type"] == "regulatory T cell"]
print("\nTreg cells:", len(treg), "| donors with Treg:", treg["donor_id"].nunique())
print(treg.groupby("donor_id").size().sort_values(ascending=False))

TREG_GENES = ["IL4", "IL13", "GATA3", "FOXP3", "IL5"]
pb = donor_pseudobulk("regulatory T cell", TREG_GENES)
pb.to_csv(OUTT + r"\ST139_treg_donor_pseudobulk.csv", index=False)

res = []
for g in TREG_GENES:
    r = poisson_age(pb[pb["gene"] == g])
    r["gene"] = g
    res.append(r)
resdf = pd.DataFrame(res)

# descriptive: detection fraction by donor (cells with count>0 among Treg)
det = []
sub = long[long["cell_type"] == "regulatory T cell"]
ncell_treg = treg.groupby("donor_id").size()
for g in TREG_GENES:
    pos = sub[(sub["gene"] == g) & (sub["count"] > 0)].groupby("cell").size()
    for d in ncell_treg.index:
        det.append(dict(gene=g, donor=d, age=float(dage[d]), n_cells=int(ncell_treg[d]),
                        n_positive=int(pos.reindex(treg[treg["donor_id"] == d].index, fill_value=0).sum())))
detdf = pd.DataFrame(det)
detdf["detection"] = detdf["n_positive"] / detdf["n_cells"]
detdf.to_csv(OUTT + r"\ST139b_treg_detection_by_donor.csv", index=False)

# Spearman on donor rates
sp = []
for g in TREG_GENES:
    d = pb[pb["gene"] == g]
    rate = d["gene_counts"] / d["total_counts"]
    rho, p = stats.spearmanr(d["age"], rate)
    sp.append(dict(gene=g, spearman_rho=rho, spearman_p=p))
resdf = resdf.merge(pd.DataFrame(sp), on="gene")

# BH-FDR within Treg family
def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q
resdf["q_treg_family"] = bh(resdf["p"].values)
resdf.to_csv(OUTT + r"\ST140_treg_age_glm_results.csv", index=False)
print("\nTreg family results:")
print(resdf.to_string())

# ---------- IL18 / TNFSF14 cell source ----------
src = []
for g in ["IL18", "TNFSF14", "IL4", "TSLP", "IL33"]:
    s = long[long["gene"] == g]
    npos = s[s["count"] > 0].groupby("cell_type")["cell"].nunique()
    ntot = obs.groupby("cell_type").size()
    gc = s.groupby("cell_type")["count"].sum()
    for ct in ntot.index:
        src.append(dict(gene=g, cell_type=ct, n_cells=int(ntot[ct]),
                        n_expressing=int(npos.get(ct, 0)),
                        pct_expressing=100 * npos.get(ct, 0) / ntot[ct],
                        total_counts=float(gc.get(ct, 0.0))))
srcdf = pd.DataFrame(src).sort_values(["gene", "pct_expressing"], ascending=[True, False])
srcdf.to_csv(OUTT + r"\ST141_cell_source_expression.csv", index=False)
print("\nTop expressing cell types:")
print(srcdf[srcdf["n_cells"] > 500].groupby("gene").head(4).to_string())

# ---------- IL18 / TNFSF14 age slopes in dominant compartments ----------
il18_ct = srcdf[(srcdf["gene"] == "IL18") & (srcdf["n_cells"] > 500)].head(2)["cell_type"].tolist()
lt_ct = srcdf[(srcdf["gene"] == "TNFSF14") & (srcdf["n_cells"] > 500)].head(2)["cell_type"].tolist()
res2 = []
for ct in il18_ct + lt_ct:
    g = "IL18" if ct in il18_ct else "TNFSF14"
    pb2 = donor_pseudobulk(ct, [g])
    pb2.to_csv(OUTT + (r"\ST142_%s_%s_pseudobulk.csv" % (g, ct.replace(" ", "_").replace(",", "")[:30])), index=False)
    r = poisson_age(pb2)
    r.update(gene=g, cell_type=ct)
    res2.append(r)
res2df = pd.DataFrame(res2)
res2df["q_translational_family"] = bh(res2df["p"].values)
res2df.to_csv(OUTT + r"\ST143_translational_age_glm_results.csv", index=False)
print("\nTranslational family:")
print(res2df.to_string())

print("\nDONE")
