#!/usr/bin/env python3
"""
R5-08: Age x index interaction in pooled NHANES D (2005-2006) + J (2017-2018) cycles.

(a) Rebuild index I in both cycles with the same fixed R3 component weights; components
    standardized against fixed J-cycle reference means/SDs (common absolute scale across
    cycles; J is the weight-derivation cycle). I is then z-scored on the pooled analytic
    sample. Weights: WTMEC4YR = WTMEC2YR / 2 (NHANES multi-cycle convention); cycle
    indicator added. (Original 1999-2018 all-cycle plan downgraded to P2: only D and J
    cycles are available locally.)
(b) Pooled weighted logistic: allergy ~ I * age(centered) + sex + bmi + cycle.
(c) Age-stratified ORs (<40, 40-60, >60 and per-10-year) + marginal effect curve of the
    I OR as a function of age.

Variance: statsmodels GLM(freq_weights).fit(cov_type='cluster', groups=strata x PSU) —
design-consistent approximation (samplics install failed; see R5-07 header / notes).
NHANES SDMVPSU is 1/2 within each SDMVSTRA, so cluster id = strata x PSU (30 clusters).
"""

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

SEED = 42

BASE_DIR = Path(__file__).resolve().parents[2]
NHANES_DIR = BASE_DIR / "data" / "real" / "nhanes"
WEIGHTS_JSON = BASE_DIR / "results" / "r3" / "tables" / "r3_index_weights.json"
NHANES_CSV = NHANES_DIR / "nhanes_merged.csv"

OUT_DIR = BASE_DIR / "results" / "r5"
TBL_DIR, FIG_DIR, LOG_DIR = OUT_DIR / "tables", OUT_DIR / "figures", OUT_DIR / "logs"
for d in (TBL_DIR, FIG_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "r5_08_age_interaction_pooled.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("r5_08")

COMPONENTS = ["crp", "wbc", "neutrophil_pct", "lymphocyte_pct", "monocyte_pct", "eosinophil_pct"]
CONSOLIDATED = []


def add_psu_cluster(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["psu_cluster"] = (df["SDMVSTRA"].astype(int) * 10 + df["SDMVPSU"].astype(int)).astype(int)
    return df


def load_fixed_index_weights() -> dict:
    if WEIGHTS_JSON.exists():
        with open(WEIGHTS_JSON) as f:
            w = json.load(f)["weights"]
        log.info(f"Loaded fixed index weights from {WEIGHTS_JSON}")
        return {c: w[c] for c in COMPONENTS}
    log.warning("r3_index_weights.json missing; recomputing from merged CSV")
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler
    df = pd.read_csv(NHANES_CSV)
    sub = df[COMPONENTS].dropna()
    X = StandardScaler().fit_transform(sub)
    lr = LinearRegression().fit(X, df.loc[sub.index, "I"].values)
    return dict(zip(COMPONENTS, lr.coef_.tolist()))


def j_cycle_component_stats() -> dict:
    df = pd.read_csv(NHANES_CSV)
    sub = df[COMPONENTS].dropna()
    return {c: {"mean": float(sub[c].mean()), "sd": float(sub[c].std(ddof=0))} for c in COMPONENTS}


def compute_fixed_index(df: pd.DataFrame, weights: dict, ref_stats: dict) -> pd.Series:
    out = pd.Series(np.nan, index=df.index)
    valid = df[COMPONENTS].dropna().index
    val = np.zeros(len(valid))
    for c in COMPONENTS:
        val += weights[c] * (df.loc[valid, c].values - ref_stats[c]["mean"]) / ref_stats[c]["sd"]
    out.loc[valid] = val
    return out


def fit_wglm(df, formula, family, weight_col, cluster_col="psu_cluster"):
    return smf.glm(formula=formula, data=df, family=family,
                   freq_weights=df[weight_col]).fit(
        cov_type="cluster", cov_kwds={"groups": df[cluster_col]})


def fit_uglm(df, formula, family):
    return smf.glm(formula=formula, data=df, family=family).fit()


def build_j_frame() -> pd.DataFrame:
    df = pd.read_csv(NHANES_CSV)
    demo = pd.read_sas(NHANES_DIR / "DEMO_J.xpt", format="xport")
    demo = demo[["SEQN", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]].copy()
    demo["SEQN"] = demo["SEQN"].astype(int)
    df["SEQN"] = df["SEQN"].astype(int)
    df = df.merge(demo, on="SEQN", how="left")
    df["sex_female"] = (df["sex"] == "F").astype(int) if df["sex"].dtype == object \
        else (df["sex"] == 2).astype(int)
    df["cycle"] = "J"
    return df[["SEQN", "age", "sex_female", "bmi", "I", "allergy_composite",
               "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "cycle"]]


def build_d_frame(weights: dict, ref_stats: dict) -> pd.DataFrame:
    cbc = pd.read_sas(NHANES_DIR / "CBC_D.xpt", format="xport")
    cbc = cbc[["SEQN", "LBXWBCSI", "LBXNEPCT", "LBXLYPCT", "LBXMOPCT", "LBXEOPCT"]].copy()
    cbc.columns = ["SEQN", "wbc", "neutrophil_pct", "lymphocyte_pct", "monocyte_pct", "eosinophil_pct"]
    crp = pd.read_sas(NHANES_DIR / "CRP_D.xpt", format="xport")[["SEQN", "LBXCRP"]].copy()
    crp.columns = ["SEQN", "crp"]
    demo = pd.read_sas(NHANES_DIR / "DEMO_D.xpt", format="xport")
    demo = demo[["SEQN", "RIAGENDR", "RIDAGEYR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]].copy()
    demo.columns = ["SEQN", "sex", "age", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]
    bmx = pd.read_sas(NHANES_DIR / "BMX_D.xpt", format="xport")[["SEQN", "BMXBMI"]].copy()
    bmx.columns = ["SEQN", "bmi"]
    mcq = pd.read_sas(NHANES_DIR / "MCQ_D.xpt", format="xport")[["SEQN", "MCQ010", "MCQ080"]].copy()
    mcq.columns = ["SEQN", "asthma_ever", "hay_fever_ever"]

    df = cbc.merge(crp, on="SEQN").merge(demo, on="SEQN").merge(bmx, on="SEQN").merge(mcq, on="SEQN", how="left")
    df["I"] = compute_fixed_index(df, weights, ref_stats)
    df["sex_female"] = (df["sex"] == 2).astype(int)

    def allergy(row):
        a, h = row["asthma_ever"], row["hay_fever_ever"]
        if a == 1 or h == 1:
            return 1.0
        if (a == 2 or pd.isna(a)) and (h == 2 or pd.isna(h)) and not (pd.isna(a) and pd.isna(h)):
            return 0.0
        return np.nan
    df["allergy_composite"] = df.apply(allergy, axis=1)
    df["cycle"] = "D"
    return df[["SEQN", "age", "sex_female", "bmi", "I", "allergy_composite",
               "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "cycle"]]


def extract_or(m, term, label, n, group=None) -> dict:
    try:
        b, se, p = float(m.params[term]), float(m.bse[term]), float(m.pvalues[term])
    except Exception:
        b, se, p = np.nan, np.nan, np.nan
    return {"model": label, "group": group, "predictor": term, "n": n,
            "beta": b, "se": se, "or": np.exp(b),
            "or_lower": np.exp(b - 1.96 * se), "or_upper": np.exp(b + 1.96 * se), "p": p}


def consolidate(analysis, group, predictor, model, r, weight, note=""):
    CONSOLIDATED.append({
        "analysis": analysis, "group": group, "predictor": predictor, "model": model,
        "beta": r["beta"], "se": r["se"], "p": r["p"], "n": r["n"],
        "n_events": r.get("n_events"),
        "effect": r["or"], "ci_lower": r["or_lower"], "ci_upper": r["or_upper"],
        "effect_scale": "or", "weight": weight, "note": note,
    })


def main():
    log.info("R5-08 age interaction + pooled D+J analysis starting (seed=%d)", SEED)
    weights = load_fixed_index_weights()
    ref_stats = j_cycle_component_stats()
    log.info("Standardization: components z-scored against fixed J-cycle reference means/SDs "
             "in BOTH cycles (common absolute scale; J is the weight-derivation cycle); "
             "I then z-scored on the pooled analytic sample for OR-per-SD reporting. "
             "All-cycle 1999-2018 pooling downgraded to P2 (only D+J available locally).")

    # (a) pooled frame
    try:
        jf = build_j_frame()
        dd = build_d_frame(weights, ref_stats)
        pooled = pd.concat([jf, dd], ignore_index=True)
        pooled = pooled.dropna(subset=["I", "allergy_composite", "age", "sex_female", "bmi",
                                       "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]).copy()
        pooled = pooled[pooled["WTMEC2YR"] > 0].copy()
        pooled["WTMEC4YR"] = pooled["WTMEC2YR"] / 2.0
        pooled["I_z"] = (pooled["I"] - pooled["I"].mean()) / pooled["I"].std(ddof=0)
        pooled["age_c"] = pooled["age"] - pooled["age"].mean()
        pooled["cycle_D"] = (pooled["cycle"] == "D").astype(int)
        pooled = add_psu_cluster(pooled)
        nJ, nD = int((pooled["cycle"] == "J").sum()), int((pooled["cycle"] == "D").sum())
        log.info(f"Pooled analytic sample: n={len(pooled)} (J={nJ}, D={nD}); "
                 f"allergy cases={int(pooled['allergy_composite'].sum())}; "
                 f"mean age={pooled['age'].mean():.1f}; clusters={pooled['psu_cluster'].nunique()}")
    except Exception as e:
        log.error(f"(a) pooled frame construction failed: {e}")
        return

    # (b) pooled weighted interaction model
    int_rows = []
    formula = "allergy_composite ~ I_z * age_c + sex_female + bmi + cycle_D"
    terms = ["I_z", "age_c", "I_z:age_c", "sex_female", "bmi", "cycle_D"]
    mu, mw = None, None
    try:
        mu = fit_uglm(pooled, formula, sm.families.Binomial())
        for t in terms:
            r = extract_or(mu, t, "unweighted", int(mu.nobs))
            int_rows.append(r)
            consolidate("age_interaction_pooled", "all", t, "unweighted", r, "none")
        log.info(f"  unweighted interaction beta={mu.params['I_z:age_c']:.4f}, "
                 f"p={mu.pvalues['I_z:age_c']:.3g}")
    except Exception as e:
        log.error(f"  unweighted interaction fit failed: {e}")
    try:
        mw = fit_wglm(pooled, formula, sm.families.Binomial(), "WTMEC4YR")
        for t in terms:
            r = extract_or(mw, t, "weighted", int(mw.nobs))
            int_rows.append(r)
            consolidate("age_interaction_pooled", "all", t, "weighted", r, "WTMEC4YR")
        log.info(f"  weighted interaction beta={mw.params['I_z:age_c']:.4f} "
                 f"(se {mw.bse['I_z:age_c']:.4f}), p={mw.pvalues['I_z:age_c']:.3g}")
    except Exception as e:
        log.error(f"  weighted interaction fit failed: {e}")
    int_df = pd.DataFrame(int_rows)
    int_df.to_csv(TBL_DIR / "r5_age_interaction.csv", index=False)
    log.info(f"  saved r5_age_interaction.csv ({len(int_df)} rows)")

    # (c) age-stratified ORs
    broad = [("<40", 0, 40), ("40-60", 40, 60), (">60", 60, 200)]
    decade = [(f"{a}-{a+10}", a, a + 10) for a in range(0, 70, 10)] + [("70+", 70, 200)]
    strata = [("broad", g, lo, hi) for g, lo, hi in broad] + \
             [("decade", g, lo, hi) for g, lo, hi in decade]
    strat_rows = []
    for scheme, glab, lo, hi in strata:
        sub = pooled[(pooled["age"] >= lo) & (pooled["age"] < hi)]
        if len(sub) < 100 or sub["allergy_composite"].nunique() < 2:
            log.info(f"  stratum {glab}: n={len(sub)} insufficient, skipped")
            continue
        f_stratum = "allergy_composite ~ I_z + age + sex_female + bmi + cycle_D"
        nev = int(sub["allergy_composite"].sum())
        try:
            msu = fit_uglm(sub, f_stratum, sm.families.Binomial())
            r = extract_or(msu, "I_z", "unweighted", int(msu.nobs), glab)
            r["scheme"] = scheme
            r["n_events"] = nev
            strat_rows.append(r)
            consolidate("age_stratified", glab, "I_z", "unweighted", r, "none",
                        note=f"scheme={scheme}")
        except Exception as e:
            log.error(f"  stratum {glab} unweighted failed: {e}")
        try:
            msw = fit_wglm(sub, f_stratum, sm.families.Binomial(), "WTMEC4YR")
            r = extract_or(msw, "I_z", "weighted", int(msw.nobs), glab)
            r["scheme"] = scheme
            r["n_events"] = nev
            strat_rows.append(r)
            consolidate("age_stratified", glab, "I_z", "weighted", r, "WTMEC4YR",
                        note=f"scheme={scheme}")
            log.info(f"  stratum {glab} (n={len(sub)}, cases={nev}): weighted OR={r['or']:.3f} "
                     f"({r['or_lower']:.3f}-{r['or_upper']:.3f}), p={r['p']:.3g}")
        except Exception as e:
            log.error(f"  stratum {glab} weighted failed: {e}")
    strat_df = pd.DataFrame(strat_rows)
    strat_df.to_csv(TBL_DIR / "r5_age_stratified_or.csv", index=False)
    log.info(f"  saved r5_age_stratified_or.csv ({len(strat_df)} rows)")

    # power note for the consolidated report
    try:
        fit_p = mw if mw is not None else mu
        se_int = float(fit_p.bse["I_z:age_c"])
        # detectable interaction per year of age at 80% power (approx, normal approx)
        det = 2.8 * se_int
        log.info(f"  power note: interaction SE={se_int:.4f}/yr -> ~80%-power detectable "
                 f"slope change |beta|>={det:.4f}/yr (OR ratio over 40 y span = "
                 f"{np.exp(det*40):.2f})")
    except Exception as e:
        log.error(f"  power note failed: {e}")

    # marginal effect curve
    try:
        fit_for_marginal = mw if mw is not None else mu
        mean_age = float(pooled["age"].mean())
        ages = np.arange(20, 81, 1)
        ac = ages - mean_age
        b_I = fit_for_marginal.params["I_z"]
        b_int = fit_for_marginal.params["I_z:age_c"]
        cov = fit_for_marginal.cov_params().loc[["I_z", "I_z:age_c"], ["I_z", "I_z:age_c"]].values
        lin = b_I + b_int * ac
        var = cov[0, 0] + (ac ** 2) * cov[1, 1] + 2 * ac * cov[0, 1]
        se = np.sqrt(np.maximum(var, 0))
        or_c, lo_c, hi_c = np.exp(lin), np.exp(lin - 1.96 * se), np.exp(lin + 1.96 * se)

        fig, ax = plt.subplots(figsize=(7.5, 5))
        ax.plot(ages, or_c, color="darkred", lw=2, label="Marginal OR of I (per 1 SD)")
        ax.fill_between(ages, lo_c, hi_c, color="darkred", alpha=0.15, label="95% CI")
        sw = strat_df[strat_df["model"] == "weighted"]
        plotted = set()
        for scheme, marker, color in [("decade", "o", "steelblue"), ("broad", "s", "darkgreen")]:
            for _, r in sw[sw["scheme"] == scheme].iterrows():
                g = r["group"]
                if g.endswith("+"):
                    mid = float(g[:-1]) + 5
                elif g.startswith("<"):
                    mid = (0 + float(g[1:])) / 2
                elif g.startswith(">"):
                    mid = float(g[1:]) + 5
                else:
                    a0, a1 = g.split("-")
                    mid = (float(a0) + float(a1)) / 2
                lab = f"Stratified OR ({scheme})" if scheme not in plotted else None
                plotted.add(scheme)
                ax.errorbar(mid, r["or"], yerr=[[r["or"] - r["or_lower"]], [r["or_upper"] - r["or"]]],
                            fmt=marker, color=color, ecolor="gray", capsize=3, markersize=6, label=lab)
        ax.axhline(1.0, color="grey", ls="--", lw=0.8)
        ax.set_xlabel("Age (years)")
        ax.set_ylabel("OR of allergy per 1 SD of I")
        ax.set_title("Age-dependent marginal effect of immunosenescence index\n"
                     "Pooled NHANES 2005-2006 + 2017-2018 (weighted)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "r5_age_interaction_marginal.png", dpi=300)
        plt.close(fig)
        log.info("  saved r5_age_interaction_marginal.png")
    except Exception as e:
        log.error(f"  marginal effect curve failed: {e}")

    # forest plot of stratified weighted ORs
    try:
        sw = strat_df[strat_df["model"] == "weighted"].copy()
        order = {g: i for i, (_, g, _, _) in enumerate(strata)}
        sw["ord"] = sw["group"].map(order)
        sw = sw.sort_values("ord")
        fig, ax = plt.subplots(figsize=(8, max(4.5, 0.4 * len(sw))))
        y_pos = np.arange(len(sw))
        colors = ["darkgreen" if s == "broad" else "steelblue" for s in sw["scheme"]]
        for i, (_, r) in enumerate(sw.iterrows()):
            ax.errorbar(r["or"], y_pos[i],
                        xerr=[[r["or"] - r["or_lower"]], [r["or_upper"] - r["or"]]],
                        fmt="o", color=colors[i], ecolor="gray", capsize=3)
        ax.axvline(1.0, color="red", ls="--", lw=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{r['group']}  (n={r['n']})" for _, r in sw.iterrows()], fontsize=8)
        ax.set_xlabel("Weighted OR of allergy per 1 SD of I (95% CI)")
        ax.set_title("Age-stratified association of I with allergy\nPooled NHANES D+J (weighted)")
        ax.set_xlim(left=0)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "r5_age_stratified_forest.png", dpi=300)
        plt.close(fig)
        log.info("  saved r5_age_stratified_forest.png")
    except Exception as e:
        log.error(f"  forest plot failed: {e}")

    cons = pd.DataFrame(CONSOLIDATED)
    cons.to_csv(TBL_DIR / "r5_age_interaction_results.csv", index=False)
    log.info(f"saved consolidated r5_age_interaction_results.csv ({len(cons)} rows)")
    log.info("R5-08 complete.")


if __name__ == "__main__":
    main()
