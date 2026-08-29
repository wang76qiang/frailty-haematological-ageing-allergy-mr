#!/usr/bin/env python3
"""
R5-03 Heterogeneity governance: from defect to discovery.

(a) MR-PRESSO-style outlier governance on the 93-SNP index instruments
    (global Cochran-Q, per-SNP outlier test, distortion test) -> r5_presso_outliers.csv
(b) MR-RAPS (Huber robust IVW, overdispersion phi) vs IVW/Egger/weighted median
    -> r5_raps_vs_ivw.csv
(c) Mechanism: cluster the 93 SNPs by their 6-component effect vectors
    (Gaussian mixture, BIC-selected), annotate clusters with nearest genes,
    DICE cell-type eQTLs and Blueprint ChromHMM states, and run subgroup MR
    -> r5_heterogeneity_snp_clusters.csv, r5_cluster_annotation_enrichment.csv,
       r5_cluster_mechanism.png
"""

import os
import sys
import traceback
from collections import Counter

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import r5_common as C

log = C.log
OUTCOMES = C.FINNGEN_OUTCOMES


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_harmonised() -> dict:
    """Reproduce the R4 93-SNP harmonised index->allergy data (cached)."""
    inst = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables",
                                    "r4_index_instrument_diagnostics.csv"))
    sample_sizes = C.r1_utils.load_finngen_sample_sizes()
    out = {}
    for outcome, path in OUTCOMES.items():
        cache_name = f"r5_cache_harmonised_index_{outcome}.csv"
        har = C.cache_get(cache_name)
        if har is None:
            har = C.r4_utils.harmonise_index_with_outcome(
                inst[["rsid", "ref", "alt", "beta_I", "se_I"]],
                path, sample_sizes=sample_sizes, n_exp=400000)
            if har is not None and not har.empty:
                C.cache_put(cache_name, har)
        if har is None or har.empty:
            log(f"  WARNING: no harmonised data for {outcome}")
            continue
        out[outcome] = har
        log(f"  {outcome}: {len(har)} harmonised SNPs")
    return inst, out


# ---------------------------------------------------------------------------
# (a) PRESSO governance
# ---------------------------------------------------------------------------
def presso_governance(harmonised: dict) -> pd.DataFrame:
    rows = []
    for outcome, har in harmonised.items():
        bx = har["beta"].values.astype(float)
        by = har["beta_outcome"].values.astype(float)
        sy = har["se_outcome"].values.astype(float)
        k = len(har)
        b_full, se_full, _ = C.r1_utils.mr_ivw(bx, by, sy, random=False)
        Q, Qp = C.r1_utils.cochran_q(bx, by, sy, b_full)
        keep, n_iter = C.r1_utils.mr_presso_outliers(bx, by, sy, alpha=0.05)
        # per-SNP residual stats under the full IVW estimate
        z = (by - bx * b_full) / sy
        p_raw = 2 * (1 - stats.norm.cdf(np.abs(z)))
        p_bonf = np.minimum(p_raw * k, 1.0)
        # outlier-removed estimate + distortion test
        if 2 <= keep.sum() < k:
            b_c, se_c, p_c = C.r1_utils.mr_ivw(bx[keep], by[keep], sy[keep])
            d_z = (b_full - b_c) / np.sqrt(se_full ** 2 + se_c ** 2)
            distortion_p = 2 * (1 - stats.norm.cdf(abs(d_z)))
        else:
            b_c, distortion_p = b_full, np.nan
        log(f"  {outcome}: Q={Q:.1f} (p={Qp:.2g}), outliers={int(k - keep.sum())}/{k}, "
            f"IVW {b_full:.4f} -> PRESSO-corrected {b_c:.4f} "
            f"(distortion p={distortion_p:.3g})")
        for i in range(k):
            rows.append({
                "outcome": outcome, "snp": har["snp"].iloc[i],
                "beta_exp": bx[i], "se_exp": har["se"].iloc[i],
                "beta_out": by[i], "se_out": sy[i],
                "wald_ratio": by[i] / bx[i],
                "resid_z": z[i], "p_raw": p_raw[i], "p_bonferroni": p_bonf[i],
                "is_outlier": bool(not keep[i]),
                "cochran_Q": Q, "cochran_Q_p": Qp,
                "ivw_full_beta": b_full, "ivw_presso_corrected_beta": b_c,
                "distortion_p": distortion_p, "n_snps": k,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (b) MR-RAPS (Huber robust IVW)
# ---------------------------------------------------------------------------
def mr_raps(bx, by, sy, k_huber: float = 1.345, max_iter: int = 200, tol: float = 1e-10):
    """Robust IVW with Huber loss (IRLS) and overdispersion phi."""
    bx = np.asarray(bx, float); by = np.asarray(by, float); sy = np.asarray(sy, float)
    w0 = 1.0 / sy ** 2
    b = np.sum(w0 * bx * by) / np.sum(w0 * bx ** 2)
    for _ in range(max_iter):
        r = (by - b * bx) / sy
        u = np.where(np.abs(r) <= k_huber, 1.0, k_huber / np.abs(r))
        w = w0 * u
        b_new = np.sum(w * bx * by) / np.sum(w * bx ** 2)
        if abs(b_new - b) < tol:
            b = b_new
            break
        b = b_new
    r = (by - b * bx) / sy
    u = np.where(np.abs(r) <= k_huber, 1.0, k_huber / np.abs(r))
    w = w0 * u
    df = max(len(bx) - 1, 1)
    phi = max(1.0, float(np.sum(u * r ** 2)) / df)
    se = np.sqrt(phi / np.sum(w * bx ** 2))
    p = 2 * (1 - stats.norm.cdf(abs(b / se)))
    return b, se, p, phi


def raps_comparison(harmonised: dict) -> pd.DataFrame:
    rows = []
    for outcome, har in harmonised.items():
        bx = har["beta"].values.astype(float)
        by = har["beta_outcome"].values.astype(float)
        sy = har["se_outcome"].values.astype(float)
        k = len(har)
        for r in C.mr_battery(har):
            r.update({"outcome": outcome})
            rows.append(r)
        b, se, p, phi = mr_raps(bx, by, sy)
        rows.append({"outcome": outcome, "method": "MR_RAPS", "n_snps": k,
                     "beta": b, "se": se, "p": p,
                     "or": np.exp(b), "or_lower": np.exp(b - 1.96 * se),
                     "or_upper": np.exp(b + 1.96 * se), "phi": phi})
        log(f"  {outcome}: RAPS beta={b:.4f} (se={se:.4f}, p={p:.3g}, phi={phi:.3f})")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (c) SNP clustering by component-effect vectors
# ---------------------------------------------------------------------------
COMPONENTS = ["crp", "wbc", "neutrophil_pct", "lymphocyte_pct",
              "monocyte_pct", "eosinophil_pct"]


def load_component_matrix(rsids: set) -> pd.DataFrame:
    wide = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables",
                                    "r4_component_mvmr_betas_wide.csv"))
    wide = wide[wide["rsid"].isin(rsids)].drop_duplicates("rsid").reset_index(drop=True)
    return wide


def gmm_cluster(Z: np.ndarray, max_k: int = 4):
    from sklearn.mixture import GaussianMixture
    best = None
    for k in range(1, max_k + 1):
        gm = GaussianMixture(n_components=k, covariance_type="full",
                             random_state=C.RNG_SEED, n_init=5, max_iter=500)
        gm.fit(Z)
        bic = gm.bic(Z)
        log(f"    GMM k={k}: BIC={bic:.1f}")
        if best is None or bic < best[1]:
            best = (gm, bic, k)
    return best[0], best[2]


def cluster_stability(Z: np.ndarray, k: int, labels: np.ndarray, n_boot: int = 20) -> float:
    from sklearn.mixture import GaussianMixture
    if k < 2:
        return 1.0
    rng = np.random.default_rng(C.RNG_SEED)
    jaccs = []
    n = len(Z)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            gm = GaussianMixture(n_components=k, covariance_type="full",
                                 random_state=C.RNG_SEED, n_init=3)
            lab_b = gm.fit_predict(Z[idx])
        except Exception:
            continue
        # match bootstrap clusters to reference by majority overlap
        for c in range(k):
            ref_members = set(np.where(labels == c)[0])
            if not ref_members:
                continue
            best_j = 0.0
            for cb in range(k):
                boot_members = set(idx[np.where(lab_b == cb)[0]].tolist())
                if not boot_members:
                    continue
                inter = len(ref_members & boot_members)
                union = len(ref_members | boot_members)
                best_j = max(best_j, inter / union if union else 0.0)
            jaccs.append(best_j)
    return float(np.mean(jaccs)) if jaccs else 0.0


def run_clustering(wide: pd.DataFrame):
    beta_cols = [f"{c}_beta" for c in COMPONENTS]
    X = wide[beta_cols].values.astype(float)
    Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    log("  Gaussian mixture BIC selection:")
    gm, k = gmm_cluster(Z, max_k=4)
    labels = gm.predict(Z)
    stability = cluster_stability(Z, k, labels)
    log(f"  BIC-selected k={k}, bootstrap Jaccard stability={stability:.2f}")
    method = "gmm_bic"
    if k >= 2 and stability < 0.6:
        log("  Clusters unstable (Jaccard<0.6) -> fallback: eosinophil-effect dichotomy")
        labels = (wide["eosinophil_pct_beta"].values < 0).astype(int)
        k = 2
        method = "eosinophil_sign_fallback"
    out = wide.copy()
    out["cluster"] = labels
    for j, c in enumerate(COMPONENTS):
        out[f"{c}_z"] = Z[:, j]
    out["cluster_method"] = method
    return out, Z, k, method, stability


# ---------------------------------------------------------------------------
# Cluster annotation: nearest genes, DICE eQTL, Blueprint ChromHMM
# ---------------------------------------------------------------------------
DICE_DIR = os.path.join(C.DATA_DIR, "dice_eqtl")
DICE_P_THRESH = 1e-4  # DICE cis-eQTL suggestive threshold


def scan_dice(rsids: set) -> pd.DataFrame:
    """Line-scan the 15 DICE cell VCFs for the cluster SNPs (cached)."""
    cached = C.cache_get("r5_cache_dice_cluster_snps.csv")
    if cached is not None:
        return cached
    rows = []
    for fn in sorted(os.listdir(DICE_DIR)):
        if not fn.endswith(".vcf"):
            continue
        cell = fn.replace(".vcf", "")
        path = os.path.join(DICE_DIR, fn)
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8 or parts[2] not in rsids:
                    continue
                info = {}
                for kv in parts[7].split(";"):
                    if "=" in kv:
                        kk, vv = kv.split("=", 1)
                        info[kk] = vv
                try:
                    p = float(info.get("Pvalue", "nan"))
                except ValueError:
                    continue
                rows.append({"snp": parts[2], "cell": cell,
                             "gene": info.get("GeneSymbol", ""), "p": p})
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["snp", "cell", "gene", "p"])
    C.cache_put("r5_cache_dice_cluster_snps.csv", df)
    return df


def bim_grch37_coords(rsids: set) -> dict:
    """rsid -> (chrom, pos GRCh37) from the 1000G EUR bim."""
    lookup = C.r1_utils.build_snp_lookup(list(rsids))
    return {r: (v[1], v[2]) for r, v in lookup.items()}


def scan_chromhmm(rsids: set) -> pd.DataFrame:
    """Overlap SNPs (GRCh37 via bim) with Blueprint healthy-blood segmentations.

    States are reported as raw E1..E12 labels (Blueprint 12-state model legend
    is not shipped locally; interpretation left to the report). Cached.
    """
    cached = C.cache_get("r5_cache_chromhmm_cluster_snps.csv")
    if cached is not None:
        return cached
    import bisect
    coords = bim_grch37_coords(rsids)
    seg_dir = os.path.join(C.DATA_DIR, "blueprint", "SEGMENTATION",
                           "SEGMENTATION_healthy")
    rows = []
    files = [f for f in os.listdir(seg_dir) if f.endswith(".bed")]
    for fn in files:
        sample = fn.split("_")[0]
        intervals = {}
        with open(os.path.join(seg_dir, fn)) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4 or not p[0].startswith("chr"):
                    continue
                chrom = p[0].replace("chr", "")
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                intervals.setdefault(chrom, []).append((s, e, p[3]))
        # sort and build bisect index per chrom
        starts = {}
        for chrom, ivs in intervals.items():
            ivs.sort()
            starts[chrom] = [iv[0] for iv in ivs]
        for rsid, (chrom, pos) in coords.items():
            ivs = intervals.get(chrom)
            if not ivs:
                continue
            st = starts[chrom]
            i = bisect.bisect_right(st, pos) - 1
            if i >= 0 and ivs[i][0] <= pos < ivs[i][1]:
                rows.append({"snp": rsid, "sample": sample, "state": ivs[i][2]})
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["snp", "sample", "state"])
    C.cache_put("r5_cache_chromhmm_cluster_snps.csv", df)
    return df


def bh_fdr(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank, i in enumerate(order[::-1]):
        running = min(running, p[i] * n / (n - rank))
        q[i] = running
    return np.clip(q, 0, 1)


def annotate_clusters(clusters: pd.DataFrame, harmonised: dict) -> pd.DataFrame:
    rsids = set(clusters["rsid"])
    k = clusters["cluster"].nunique()
    rows = []

    # --- component effect profile per cluster ---
    for c in sorted(clusters["cluster"].unique()):
        sub = clusters[clusters["cluster"] == c]
        for comp in COMPONENTS:
            rows.append({"cluster": int(c), "category": "component_mean_z",
                         "item": comp, "value": sub[f"{comp}_z"].mean(),
                         "p": np.nan, "fdr": np.nan})

    # --- nearest genes (FinnGen annotation) ---
    genes = C.nearest_genes_for_rsids(rsids)
    C.cache_put("r5_cache_nearest_genes_cluster_snps.csv",
                pd.DataFrame([{"snp": s, "gene": g} for s, g in genes.items()]))
    for c in sorted(clusters["cluster"].unique()):
        sub = clusters[clusters["cluster"] == c]
        cnt = Counter()
        for s in sub["rsid"]:
            for g in str(genes.get(s, "")).split(","):
                g = g.strip()
                if g:
                    cnt[g] += 1
        for g, n in cnt.most_common(10):
            rows.append({"cluster": int(c), "category": "nearest_gene_top10",
                         "item": g, "value": n, "p": np.nan, "fdr": np.nan})

    # --- DICE eQTL enrichment (Fisher: cluster vs rest, per cell type) ---
    dice = scan_dice(rsids)
    tests = []
    if not dice.empty:
        sig = dice[dice["p"] < DICE_P_THRESH]
        eqtl_snps_by_cell = {cell: set(g["snp"]) for cell, g in sig.groupby("cell")}
        all_snps = rsids
        for c in sorted(clusters["cluster"].unique()):
            cl = set(clusters.loc[clusters["cluster"] == c, "rsid"])
            rest = all_snps - cl
            for cell, eset in eqtl_snps_by_cell.items():
                a = len(cl & eset)
                b = len(cl - eset)
                c2 = len(rest & eset)
                d = len(rest - eset)
                odds, p = stats.fisher_exact([[a, b], [c2, d]])
                tests.append({"cluster": int(c), "category": "dice_eqtl_enrichment",
                              "item": cell, "value": a,
                              "n_cluster_snps": len(cl), "p": p,
                              "odds_ratio": odds})
    # --- ChromHMM state enrichment (Fisher per state, cluster vs rest) ---
    chrom = scan_chromhmm(rsids)
    if not chrom.empty:
        state_by_snp = chrom.groupby("snp")["state"].agg(set).to_dict()
        states = sorted(chrom["state"].unique())
        all_snps = rsids
        for c in sorted(clusters["cluster"].unique()):
            cl = set(clusters.loc[clusters["cluster"] == c, "rsid"])
            rest = all_snps - cl
            for st in states:
                in_state = {s for s, ss in state_by_snp.items() if st in ss}
                a = len(cl & in_state)
                b = len(cl - in_state)
                c2 = len(rest & in_state)
                d = len(rest - in_state)
                odds, p = stats.fisher_exact([[a, b], [c2, d]])
                tests.append({"cluster": int(c), "category": "chromhmm_state_enrichment",
                              "item": st, "value": a,
                              "n_cluster_snps": len(cl), "p": p,
                              "odds_ratio": odds})
    if tests:
        tdf = pd.DataFrame(tests)
        tdf["fdr"] = bh_fdr(tdf["p"].values)
        for _, r in tdf.iterrows():
            rows.append({"cluster": r["cluster"], "category": r["category"],
                         "item": r["item"], "value": r["value"],
                         "p": r["p"], "fdr": r["fdr"]})

    # --- subgroup MR per cluster ---
    for c in sorted(clusters["cluster"].unique()):
        cl = set(clusters.loc[clusters["cluster"] == c, "rsid"])
        for outcome, har in harmonised.items():
            sub = har[har["snp"].isin(cl)]
            if len(sub) < 3:
                continue
            bx = sub["beta"].values.astype(float)
            by = sub["beta_outcome"].values.astype(float)
            sy = sub["se_outcome"].values.astype(float)
            b, se, p = C.r1_utils.mr_ivw(bx, by, sy, random=False)
            rows.append({"cluster": int(c), "category": f"subgroup_mr_{outcome}",
                         "item": "IVW_fixed", "value": b, "se": se,
                         "or": np.exp(b), "n_snps": len(sub), "p": p, "fdr": np.nan})
            log(f"  cluster {c} x {outcome}: n={len(sub)}, IVW beta={b:.4f} "
                f"(OR={np.exp(b):.3f}, p={p:.3g})")
    return pd.DataFrame(rows)


def make_cluster_figure(clusters: pd.DataFrame, Z: np.ndarray):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = clusters.sort_values(["cluster", "eosinophil_pct_z"]).index
    Zs = clusters.loc[order, [f"{c}_z" for c in COMPONENTS]].values
    lab = clusters.loc[order, "cluster"].values

    fig, axes = plt.subplots(1, 2, figsize=(13, 7),
                             gridspec_kw={"width_ratios": [2.4, 1]})
    ax = axes[0]
    vmax = np.nanmax(np.abs(Zs))
    im = ax.imshow(Zs, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto",
                   interpolation="nearest")
    ax.set_xticks(range(len(COMPONENTS)))
    ax.set_xticklabels([c.replace("_pct", "%") for c in COMPONENTS],
                       rotation=35, ha="right", fontsize=9)
    ax.set_ylabel(f"Index SNPs (n={len(Zs)}, sorted by cluster)", fontsize=9)
    # cluster boundaries
    bounds = np.where(np.diff(lab) != 0)[0]
    for b in bounds:
        ax.axhline(b + 0.5, color="black", lw=1.0)
    for c in np.unique(lab):
        idx = np.where(lab == c)[0]
        ax.text(-0.65, idx.mean(), f"C{c} (n={len(idx)})", fontsize=9,
                va="center", ha="right", fontweight="bold")
    ax.set_title("93-SNP index instruments: standardised 6-component effects", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.7, label="z-scored beta")

    ax = axes[1]
    prof = clusters.groupby("cluster")[[f"{c}_z" for c in COMPONENTS]].mean()
    x = np.arange(len(COMPONENTS))
    for c, row in prof.iterrows():
        ax.plot(x, row.values, marker="o", label=f"C{c}")
    ax.axhline(0, color="grey", lw=0.7, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_pct", "%") for c in COMPONENTS],
                       rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("cluster mean z", fontsize=9)
    ax.set_title("Cluster effect profiles", fontsize=10)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(C.R5_FIGURES, "r5_cluster_mechanism.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved figure: {out}")
    return out


def main():
    log("=" * 72)
    log("R5-03 heterogeneity governance (PRESSO / RAPS / SNP clusters)")
    log("=" * 72)

    log("\n[1/6] Loading harmonised 93-SNP index data (R4 reproduction)...")
    inst, harmonised = load_harmonised()

    log("\n[2/6] MR-PRESSO outlier governance...")
    try:
        presso = presso_governance(harmonised)
        p = os.path.join(C.R5_TABLES, "r5_presso_outliers.csv")
        presso.to_csv(p, index=False)
        log(f"Saved: {p} ({len(presso)} rows)")
    except Exception as exc:
        log(f"  PRESSO failed: {exc}")
        traceback.print_exc()

    log("\n[3/6] MR-RAPS vs IVW/Egger/weighted median...")
    try:
        raps = raps_comparison(harmonised)
        p = os.path.join(C.R5_TABLES, "r5_raps_vs_ivw.csv")
        raps.to_csv(p, index=False)
        log(f"Saved: {p} ({len(raps)} rows)")
    except Exception as exc:
        log(f"  RAPS failed: {exc}")
        traceback.print_exc()

    log("\n[4/6] Clustering 93 SNPs by 6-component effect vectors...")
    all_rsids = set()
    for har in harmonised.values():
        all_rsids.update(har["snp"].tolist())
    wide = load_component_matrix(all_rsids)
    log(f"  component matrix: {len(wide)} SNPs x {len(COMPONENTS)} components")
    clusters, Z, k, method, stability = run_clustering(wide)

    # link to eosinophil quartile assignment (the non-monotonic flip)
    q_path = os.path.join(C.ROOT, "results", "r4", "tables",
                          "r4_eosinophil_quartile_assignment.csv")
    if os.path.exists(q_path):
        q = pd.read_csv(q_path)[["rsid", "quartile"]]
        clusters = clusters.merge(q, on="rsid", how="left")

    p = os.path.join(C.R5_TABLES, "r5_heterogeneity_snp_clusters.csv")
    clusters.to_csv(p, index=False)
    log(f"Saved: {p} ({len(clusters)} rows, k={k}, method={method})")

    log("\n[5/6] Annotating clusters (genes / DICE eQTL / ChromHMM) + subgroup MR...")
    try:
        enrich = annotate_clusters(clusters, harmonised)
        for name in ("r5_cluster_annotation_enrichment.csv",
                     "r5_heterogeneity_cluster_enrichment.csv"):
            p = os.path.join(C.R5_TABLES, name)
            enrich.to_csv(p, index=False)
            log(f"Saved: {p} ({len(enrich)} rows)")
    except Exception as exc:
        log(f"  annotation failed: {exc}")
        traceback.print_exc()

    log("\n[6/6] Drawing cluster mechanism figure...")
    try:
        out1 = make_cluster_figure(clusters, Z)
        # collision-safe copy under a unique name
        import shutil
        out2 = os.path.join(C.R5_FIGURES, "r5_heterogeneity_cluster_mechanism.png")
        shutil.copyfile(out1, out2)
        log(f"Saved copy: {out2}")
    except Exception as exc:
        log(f"  figure failed: {exc}")
        traceback.print_exc()

    log("\nR5-03 complete.")


if __name__ == "__main__":
    main()
