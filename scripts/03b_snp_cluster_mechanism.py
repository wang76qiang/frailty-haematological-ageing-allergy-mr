#!/usr/bin/env python3
"""
R5-03b SNP-cluster mechanism analysis.

Cluster the index-instrument SNPs (the 93 that harmonise with FinnGen) by their
standardised effect profile across the 6 Pan-UKBB index components
(GaussianMixture, BIC-selected 1-4 clusters, random_state=42).  For each
cluster: per-cluster index->allergy IVW for the 3 FinnGen outcomes, dominant
component effect directions, and high-frequency FinnGen nearest_genes.

Outputs:
  results/r5/tables/r5_snp_clusters.csv
  results/r5/tables/r5_cluster_annotation_enrichment.csv
  results/r5/figures/r5_cluster_mechanism.png
"""

import os
import sys
import traceback
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r5x_utils as C

OUTCOMES = list(C.FINNGEN_OUTCOMES)
COMPONENTS = ["crp", "wbc", "neutrophil_pct", "lymphocyte_pct", "monocyte_pct", "eosinophil_pct"]
N_EXP = 400000

log_lines = []


def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)


def load_harmonised_map(inst):
    """Harmonise index instruments with each FinnGen outcome."""
    sample_sizes = C.r1_utils.load_finngen_sample_sizes()
    out = {}
    for name in OUTCOMES:
        har = C.harmonise_index_with_outcome_robust(
            inst, C.FINNGEN_OUTCOMES[name], sample_sizes=sample_sizes, n_exp=N_EXP)
        if har is not None and not har.empty:
            out[name] = har.reset_index(drop=True)
            log(f"  {name}: harmonised n={len(har)}")
    return out


def load_nearest_genes(rsids):
    """rsid -> nearest_genes from a FinnGen full summary file (variant-level annotation)."""
    return C.retry_call(lambda: _load_nearest_genes_impl(rsids), label="nearest_genes")


def _load_nearest_genes_impl(rsids):
    rsids = set(map(str, rsids))
    path = C.FINNGEN_OUTCOMES["ALLERG_ASTHMA"]
    rows = []
    for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=200000,
                             usecols=lambda c: c in ["rsids", "nearest_genes"]):
        chunk["first"] = chunk["rsids"].astype(str).str.split(",").str[0]
        chunk = chunk[chunk["first"].isin(rsids)]
        if not chunk.empty:
            rows.append(chunk[["first", "nearest_genes"]])
    if not rows:
        return {}
    df = pd.concat(rows, ignore_index=True).drop_duplicates("first", keep="first")
    return dict(zip(df["first"], df["nearest_genes"].fillna("")))


def main():
    log("=" * 72)
    log("R5-03b SNP clustering by component effect profile")
    log("=" * 72)

    inst = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables",
                                    "r4_index_instrument_diagnostics.csv"))
    inst = inst[["rsid", "ref", "alt", "beta_I", "se_I"]].copy()
    log(f"Index instruments: {len(inst)} SNPs")

    log("Harmonising index instruments with FinnGen outcomes ...")
    har_map = load_harmonised_map(inst)
    if not har_map:
        log("No harmonised data; aborting.")
        return

    # SNP set for clustering = SNPs harmonised with FinnGen (the "93")
    harm_rsids = set(har_map[OUTCOMES[0]]["snp"])
    log(f"Clustering SNP set (harmonised with FinnGen): {len(harm_rsids)}")

    wide = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables",
                                    "r4_component_mvmr_betas_wide.csv"))
    wide = wide[wide["rsid"].isin(harm_rsids)].drop_duplicates("rsid").reset_index(drop=True)
    beta_cols = [f"{c}_beta" for c in COMPONENTS]
    beta_cols = [c for c in beta_cols if c in wide.columns]
    wide = wide.dropna(subset=beta_cols).reset_index(drop=True)
    log(f"SNPs with complete 6-component betas: {len(wide)}")

    # ---- clustering ----
    X = wide[beta_cols].values.astype(float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    best, best_bic = None, np.inf
    bics = {}
    for k in range(1, 5):
        try:
            gmm = GaussianMixture(n_components=k, covariance_type="full",
                                  random_state=42, n_init=5, max_iter=500)
            gmm.fit(Xz)
            bic = gmm.bic(Xz)
            bics[k] = bic
            log(f"  GMM k={k}: BIC={bic:.1f}")
            if bic < best_bic:
                best_bic, best = bic, gmm
        except Exception:
            log(f"  GMM k={k} failed:\n{traceback.format_exc()}")
    if best is None:
        log("GMM failed entirely; aborting.")
        return
    labels = best.predict(Xz)
    n_clusters = best.n_components
    log(f"Selected k={n_clusters} (BIC={best_bic:.1f}); cluster sizes: "
        f"{dict(Counter(labels))}")

    wide["cluster"] = labels
    snp_out = wide[["rsid", "cluster"] + beta_cols].copy()
    path = os.path.join(C.TABLES_DIR, "r5_snp_clusters.csv")
    snp_out.to_csv(path, index=False)
    log(f"Saved {path} ({len(snp_out)} rows)")

    # ---- gene annotation ----
    gene_map = load_nearest_genes(wide["rsid"].tolist())

    # ---- per-cluster characterisation ----
    comp_short = {f"{c}_beta": c for c in COMPONENTS if f"{c}_beta" in beta_cols}
    enrich_rows = []
    for cl in sorted(set(labels)):
        sub = wide[wide["cluster"] == cl]
        row = {"cluster": int(cl), "n_snps": len(sub)}
        # dominant component directions (mean z per component)
        zsub = Xz[labels == cl]
        means = zsub.mean(axis=0)
        ranked = sorted(zip(beta_cols, means), key=lambda t: -abs(t[1]))
        for bc, m in zip(beta_cols, means):
            row[f"mean_z_{comp_short[bc]}"] = float(m)
        row["dominant_components"] = "; ".join(
            f"{comp_short[bc]}({'+' if m >= 0 else '-'})" for bc, m in ranked[:3])
        # high-frequency nearest genes
        genes = []
        for rs in sub["rsid"]:
            g = gene_map.get(rs, "")
            genes.extend([x for x in str(g).split(",") if x and x != "nan"])
        top = Counter(genes).most_common(8)
        row["top_genes"] = "; ".join(f"{g}({n})" for g, n in top)
        # per-cluster IVW per outcome
        snps = set(sub["rsid"])
        for out in OUTCOMES:
            har = har_map.get(out)
            if har is None:
                continue
            h = har[har["snp"].isin(snps)]
            row[f"n_snps_{out}"] = len(h)
            if len(h) >= 3:
                b, se, p = C.r1_utils.mr_ivw(h["beta"].values, h["beta_outcome"].values,
                                             h["se_outcome"].values, random=False)
                row[f"ivw_beta_{out}"] = b
                row[f"ivw_se_{out}"] = se
                row[f"ivw_p_{out}"] = p
                row[f"ivw_or_{out}"] = np.exp(b)
            else:
                row[f"ivw_beta_{out}"] = np.nan
                row[f"ivw_p_{out}"] = np.nan
                row[f"ivw_or_{out}"] = np.nan
        enrich_rows.append(row)
        log(f"  cluster {cl}: n={len(sub)}, dominant={row['dominant_components']}, "
            f"top genes={row['top_genes'][:80]}")

    enrich = pd.DataFrame(enrich_rows)
    path = os.path.join(C.TABLES_DIR, "r5_cluster_annotation_enrichment.csv")
    enrich.to_csv(path, index=False)
    log(f"Saved {path}")

    # ---- figure: clustered heatmap ----
    try:
        order = wide.sort_values(["cluster", "rsid"]).index
        mat = pd.DataFrame(Xz, columns=[comp_short[c] for c in beta_cols], index=wide["rsid"])
        mat = mat.loc[wide.loc[order, "rsid"]]
        cl_vals = wide.loc[order, "cluster"].values
        palette = sns.color_palette("Set2", n_clusters)
        row_colors = pd.Series([palette[c] for c in cl_vals], index=mat.index, name="cluster")
        g = sns.clustermap(mat, row_cluster=False, col_cluster=True, cmap="RdBu_r",
                           center=0, row_colors=row_colors, figsize=(9, 12),
                           xticklabels=True, yticklabels=False,
                           cbar_kws={"label": "z-scored component beta"})
        g.fig.suptitle(f"Index-instrument SNP clusters (GMM BIC k={n_clusters})", y=1.01)
        handles = [plt.Rectangle((0, 0), 1, 1, color=palette[c]) for c in range(n_clusters)]
        g.ax_heatmap.legend(handles, [f"cluster {c} (n={int((labels == c).sum())})"
                                      for c in range(n_clusters)],
                            loc="upper left", bbox_to_anchor=(1.25, 1.0), title="cluster")
        out = os.path.join(C.FIGURES_DIR, "r5_cluster_mechanism.png")
        g.fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(g.fig)
        log(f"Saved {out}")
    except Exception:
        log(f"Figure FAILED\n{traceback.format_exc()}")

    with open(os.path.join(C.LOGS_DIR, "r5_03b_run.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    log("R5-03b complete.")


if __name__ == "__main__":
    main()
