#!/usr/bin/env python3
"""
R5-07: NHANES complex-survey weighted analyses.

(a) J-cycle (2017-2018) weighted logistic main analysis:
        allergy_composite ~ I + age + sex + bmi   (WTMEC2YR; strata x PSU cluster SE)
(b) D-cycle (2005-2006) weighted IgE validation (index rebuilt with fixed R3 weights):
        log total IgE / high IgE(>100) / any specific-IgE+ / self-reported allergy among IgE+
(c) Weighted Cox model for all-cause mortality (NHANES 2017-2018 linked mortality 2019).
(d) Weighted threshold (hinge) analysis of I with stratified PSU cluster bootstrap CI.

Variance strategy
-----------------
samplics was attempted once (`pip install samplics`) but the install timed out, so all
weighted GLMs use the documented fallback:
    statsmodels GLM(..., freq_weights=w).fit(cov_type='cluster', cov_kwds={'groups': cl})
i.e. design-consistent weighted pseudo-likelihood point estimates with cluster-robust
sandwich SEs. IMPORTANT: NHANES SDMVPSU is coded 1/2 *within* each SDMVSTRA (verified:
15 strata x 2 PSUs), so the cluster id is the strata x PSU combination (30 clusters),
NOT SDMVPSU alone. A strata-fixed-effects sensitivity (C(SDMVSTRA) dummies + same
clustering) is reported for the main model. This is a design-consistent approximation
of Taylor linearization; the deviation from a full samplics/R-survey run is documented.
"""

import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
import statsmodels.api as sm
import statsmodels.formula.api as smf
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

SEED = 42
RNG = np.random.default_rng(SEED)

BASE_DIR = Path(__file__).resolve().parents[2]
NHANES_DIR = BASE_DIR / "data" / "real" / "nhanes"
MORT_FILE = BASE_DIR / "data" / "real" / "nhanes_mortality" / "NHANES_2017_2018_MORT_2019_PUBLIC.dat"
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
        logging.FileHandler(LOG_DIR / "r5_07_nhanes_survey_weighted.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("r5_07")

COMPONENTS = ["crp", "wbc", "neutrophil_pct", "lymphocyte_pct", "monocyte_pct", "eosinophil_pct"]

# IgE variable mapping (NHANES 2005-2006 AL_IGE_D codebook; verified via read_sas).
IGE_VARS = {
    "LBXIGE": "total_IgE",
    "LBXID1": "dust_mite_dp", "LBXID2": "dust_mite_df",
    "LBXIE1": "cat", "LBXIE5": "dog", "LBXII6": "cockroach",
    "LBXIM6": "alternaria", "LBXF13": "peanut", "LBXIF1": "egg",
    "LBXIF2": "milk", "LBXIW1": "ragweed", "LBXIG5": "rye_grass",
    "LBXIG2": "bermuda_grass", "LBXIT7": "oak", "LBXIT3": "birch",
    "LBXF24": "shrimp", "LBXIM3": "aspergillus", "LBXW11": "thistle",
    "LBXE72": "mouse", "LBXE74": "rat",
}

CONSOLIDATED = []  # rows for r5_nhanes_weighted_results.csv


def add_psu_cluster(df: pd.DataFrame) -> pd.DataFrame:
    """NHANES PSU (1/2) is nested within stratum -> combined cluster id."""
    df = df.copy()
    df["psu_cluster"] = (df["SDMVSTRA"].astype(int) * 10 + df["SDMVPSU"].astype(int)).astype(int)
    return df


def load_fixed_index_weights() -> dict:
    if WEIGHTS_JSON.exists():
        with open(WEIGHTS_JSON) as f:
            w = json.load(f)["weights"]
        log.info(f"Loaded fixed index weights from {WEIGHTS_JSON}")
        return {c: w[c] for c in COMPONENTS}
    log.warning("r3_index_weights.json missing; recomputing weights from merged CSV")
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


def extract_terms(m, terms, scale, model_label, n):
    rows = []
    for t in terms:
        try:
            b, se, p = float(m.params[t]), float(m.bse[t]), float(m.pvalues[t])
        except Exception:
            b, se, p = np.nan, np.nan, np.nan
        row = {"model": model_label, "predictor": t, "n": n, "beta": b, "se": se, "p": p}
        if scale == "or":
            row.update(or_=np.exp(b), or_lower=np.exp(b - 1.96 * se), or_upper=np.exp(b + 1.96 * se))
        rows.append(row)
    return rows


def consolidate(analysis, outcome, exposure, model, beta, se, p, n, effect_scale,
                n_events=None, weight="none", note=""):
    CONSOLIDATED.append({
        "analysis": analysis, "outcome": outcome, "exposure": exposure, "model": model,
        "beta": beta, "se": se, "p": p, "n": n, "n_events": n_events,
        "effect": np.exp(beta) if effect_scale in ("or", "hr") else beta,
        "ci_lower": np.exp(beta - 1.96 * se) if effect_scale in ("or", "hr") else beta - 1.96 * se,
        "ci_upper": np.exp(beta + 1.96 * se) if effect_scale in ("or", "hr") else beta + 1.96 * se,
        "effect_scale": effect_scale, "weight": weight, "note": note,
    })


def read_mortality_fixed_width(path: Path) -> pd.DataFrame:
    rows = []
    with open(path, "r") as f:
        for line in f:
            if len(line) < 48:
                continue
            rows.append({
                "SEQN": int(line[0:6].strip() or -1),
                "eligstat": int(line[14:15].strip() or -1),
                "mortstat": line[15:16].strip(),
                "permth_int": line[42:45].strip(),
            })
    df = pd.DataFrame(rows)
    for col in ["mortstat", "permth_int"]:
        df[col] = pd.to_numeric(df[col].replace("", np.nan).replace(".", np.nan), errors="coerce")
    return df


# --------------------------------------------------------------------------------------
# (a) J-cycle weighted main analysis
# --------------------------------------------------------------------------------------
def load_j_analytic() -> pd.DataFrame:
    df = pd.read_csv(NHANES_CSV)
    demo = pd.read_sas(NHANES_DIR / "DEMO_J.xpt", format="xport")
    ok = all(c in demo.columns for c in ["WTMEC2YR", "SDMVPSU", "SDMVSTRA"])
    log.info(f"DEMO_J columns verified: WTMEC2YR/SDMVPSU/SDMVSTRA present = {ok}")
    demo = demo[["SEQN", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]].copy()
    demo["SEQN"] = demo["SEQN"].astype(int)
    df["SEQN"] = df["SEQN"].astype(int)
    df = df.merge(demo, on="SEQN", how="left")
    ana = df.dropna(subset=["allergy_composite", "I", "age", "sex", "bmi",
                            "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]).copy()
    ana = ana[ana["WTMEC2YR"] > 0].copy()
    ana["I_z"] = (ana["I"] - ana["I"].mean()) / ana["I"].std(ddof=0)
    ana["sex_female"] = (ana["sex"] == "F").astype(int) if ana["sex"].dtype == object \
        else (ana["sex"] == 2).astype(int)
    ana = add_psu_cluster(ana)
    log.info(f"J-cycle analytic n={len(ana)}; I SD={ana['I'].std(ddof=0):.3f}; "
             f"clusters={ana['psu_cluster'].nunique()}, strata={ana['SDMVSTRA'].nunique()}")
    return ana


def task_a_main(ana: pd.DataFrame) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("(a) J-cycle weighted logistic main analysis")
    rows = []
    specs = [
        ("I_z", "I_z (per SD)"),
        ("I", "I (per unit, R4 scale)"),
    ]
    for xvar, xlabel in specs:
        formula = f"allergy_composite ~ {xvar} + age + sex_female + bmi"
        terms = [xvar, "age", "sex_female", "bmi"]
        for kind, fitter, wname in [
            ("unweighted", lambda f: fit_uglm(ana, f, sm.families.Binomial()), "none"),
            ("weighted", lambda f: fit_wglm(ana, f, sm.families.Binomial(), "WTMEC2YR"), "WTMEC2YR"),
        ]:
            try:
                m = fitter(formula)
                rows += extract_terms(m, terms, "or", f"{kind}|{xlabel}", int(m.nobs))
                b, se, p = m.params[xvar], m.bse[xvar], m.pvalues[xvar]
                log.info(f"  {kind:11s} {xlabel}: OR={np.exp(b):.3f} "
                         f"({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f}), p={p:.3g}")
                consolidate("main_logistic", "allergy_composite", xlabel, kind, b, se, p,
                            int(m.nobs), "or", n_events=int(ana["allergy_composite"].sum()),
                            weight=wname)
            except Exception as e:
                log.error(f"  {kind} {xlabel} fit failed: {e}")
    # strata fixed-effects sensitivity (SDMVSTRA dummies + strata x PSU cluster)
    try:
        ms = fit_wglm(ana, "allergy_composite ~ I_z + age + sex_female + bmi + C(SDMVSTRA)",
                      sm.families.Binomial(), "WTMEC2YR")
        b, se, p = ms.params["I_z"], ms.bse["I_z"], ms.pvalues["I_z"]
        rows += extract_terms(ms, ["I_z", "age", "sex_female", "bmi"], "or",
                              "weighted_strataFE|I_z (per SD)", int(ms.nobs))
        log.info(f"  weighted+strata-FE sensitivity: OR={np.exp(b):.3f} "
                 f"({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f}), p={p:.3g}")
        consolidate("main_logistic", "allergy_composite", "I_z (per SD)", "weighted_strataFE",
                    b, se, p, int(ms.nobs), "or",
                    n_events=int(ana["allergy_composite"].sum()), weight="WTMEC2YR",
                    note="SDMVSTRA fixed effects added; cluster=strata x PSU")
    except Exception as e:
        log.error(f"  strata-FE sensitivity failed: {e}")
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "r5_nhanes_weighted_vs_unweighted.csv", index=False)
    log.info(f"  saved r5_nhanes_weighted_vs_unweighted.csv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------------------
# (b) D-cycle weighted IgE validation
# --------------------------------------------------------------------------------------
def build_d_cycle(weights: dict, ref_stats: dict) -> pd.DataFrame:
    cbc = pd.read_sas(NHANES_DIR / "CBC_D.xpt", format="xport")
    cbc = cbc[["SEQN", "LBXWBCSI", "LBXNEPCT", "LBXLYPCT", "LBXMOPCT", "LBXEOPCT"]].copy()
    cbc.columns = ["SEQN", "wbc", "neutrophil_pct", "lymphocyte_pct", "monocyte_pct", "eosinophil_pct"]
    crp = pd.read_sas(NHANES_DIR / "CRP_D.xpt", format="xport")[["SEQN", "LBXCRP"]].copy()
    crp.columns = ["SEQN", "crp"]
    demo = pd.read_sas(NHANES_DIR / "DEMO_D.xpt", format="xport")
    ok = all(c in demo.columns for c in ["WTMEC2YR", "SDMVPSU", "SDMVSTRA"])
    log.info(f"DEMO_D columns verified: WTMEC2YR/SDMVPSU/SDMVSTRA present = {ok}")
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
    return df


def task_b_ige(weights: dict, ref_stats: dict) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("(b) D-cycle weighted IgE validation")
    rows = []
    try:
        dfd = build_d_cycle(weights, ref_stats)
        ige = pd.read_sas(NHANES_DIR / "AL_IGE_D.xpt", format="xport")
        has_wtsaf = "WTSAF2YR" in ige.columns
        log.info(f"AL_IGE_D verified: total IgE column = LBXIGE (not URXIGE); "
                 f"WTSAF2YR present = {has_wtsaf}")
        keep = ["SEQN"] + (["WTSAF2YR"] if has_wtsaf else [])
        lc_map = {}
        for col in ige.columns:
            if col.startswith("LBX") and col in IGE_VARS:
                lc = f"LBD{col[3:]}LC"
                if lc in ige.columns:
                    keep += [col, lc]
                    lc_map[col] = lc
        merged = dfd.merge(ige[keep].copy(), on="SEQN", how="inner")
        wcol = "WTSAF2YR" if has_wtsaf else "WTMEC2YR"
        if not has_wtsaf:
            log.warning("WTSAF2YR absent from local AL_IGE_D.xpt (no WT* columns at all); "
                        "using WTMEC2YR as documented approximation of the 1/2-sample IgE weight")
        merged = merged.rename(columns={k: v for k, v in IGE_VARS.items() if k in merged.columns})
        merged = add_psu_cluster(merged)
        log.info(f"  D-cycle merged IgE analytic frame: n={len(merged)}")

        merged["log_total_IgE"] = np.log1p(merged["total_IgE"])
        merged["high_total_IgE"] = np.where(merged["total_IgE"].notna(),
                                            (merged["total_IgE"] > 100).astype(float), np.nan)
        pos_cols = []
        for orig, label in IGE_VARS.items():
            if orig == "LBXIGE" or orig not in lc_map:
                continue
            if lc_map[orig] in merged.columns:
                c = f"pos_{label}"
                merged[c] = (merged[lc_map[orig]] < 0.5).astype(int)
                pos_cols.append(c)
        if pos_cols:
            merged["any_sensitization"] = (merged[pos_cols].sum(axis=1) > 0).astype(int)

        valid = merged.dropna(subset=["I", "age", "sex_female", "bmi", wcol, "psu_cluster"]).copy()
        valid = valid[valid[wcol] > 0].copy()
        valid["I_z"] = (valid["I"] - valid["I"].mean()) / valid["I"].std(ddof=0)
        log.info(f"  complete-case weighted frame: n={len(valid)}, weight={wcol}")

        specs = [
            ("log_total_IgE", sm.families.Gaussian(), "beta", "total_IgE_continuous"),
            ("high_total_IgE", sm.families.Binomial(), "or", "high_total_IgE_gt100"),
            ("any_sensitization", sm.families.Binomial(), "or", "any_specific_IgE_positive"),
        ]
        for y, fam, scale, label in specs:
            sub = valid.dropna(subset=[y])
            if sub[y].nunique() < 2:
                log.warning(f"  {label}: outcome degenerate, skipped")
                continue
            formula = f"{y} ~ I_z + age + sex_female + bmi"
            for kind, fitter, wname in [
                ("unweighted", lambda f, s: fit_uglm(s, f, fam), "none"),
                ("weighted", lambda f, s: fit_wglm(s, f, fam, wcol), wcol),
            ]:
                try:
                    m = fitter(formula, sub)
                    r = extract_terms(m, ["I_z"], scale, f"{label}|{kind}", int(m.nobs))
                    for rr in r:
                        rr["analysis"] = label
                        rr["weight"] = wname
                    rows += r
                    b, se, p = m.params["I_z"], m.bse["I_z"], m.pvalues["I_z"]
                    consolidate("ige_validation", label, "I_z (per SD)", kind, b, se, p,
                                int(m.nobs), scale, weight=wname)
                    if kind == "weighted":
                        eff = np.exp(b) if scale == "or" else b
                        log.info(f"  {label} weighted: effect={eff:.3f}, p={p:.3g}, n={int(m.nobs)}")
                except Exception as e:
                    log.error(f"  {label} {kind} failed: {e}")

        # self-reported allergy among IgE+
        if "any_sensitization" in valid.columns:
            sub = valid[valid["any_sensitization"] == 1].dropna(subset=["allergy_composite"])
            if sub["allergy_composite"].nunique() == 2 and len(sub) >= 30:
                formula = "allergy_composite ~ I_z + age + sex_female + bmi"
                for kind, fitter, wname in [
                    ("unweighted", lambda f, s: fit_uglm(s, f, sm.families.Binomial()), "none"),
                    ("weighted", lambda f, s: fit_wglm(s, f, sm.families.Binomial(), wcol), wcol),
                ]:
                    try:
                        m = fitter(formula, sub)
                        r = extract_terms(m, ["I_z"], "or",
                                          f"self_reported_allergy_in_IgE_pos|{kind}", int(m.nobs))
                        for rr in r:
                            rr["analysis"] = "self_reported_allergy_in_IgE_pos"
                            rr["weight"] = wname
                        rows += r
                        b, se, p = m.params["I_z"], m.bse["I_z"], m.pvalues["I_z"]
                        consolidate("ige_validation", "self_reported_allergy_in_IgE_pos",
                                    "I_z (per SD)", kind, b, se, p, int(m.nobs), "or",
                                    n_events=int(sub["allergy_composite"].sum()), weight=wname)
                        if kind == "weighted":
                            log.info(f"  allergy-in-IgE+ weighted: OR={np.exp(b):.3f}, p={p:.3g}, n={int(m.nobs)}")
                    except Exception as e:
                        log.error(f"  allergy-in-IgE+ {kind} failed: {e}")
            else:
                log.warning("  allergy-in-IgE+ skipped: insufficient variation/n")
    except Exception as e:
        log.error(f"(b) IgE validation failed: {e}")
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "r5_ige_validation_weighted.csv", index=False)
    log.info(f"  saved r5_ige_validation_weighted.csv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------------------
# (c) weighted Cox mortality
# --------------------------------------------------------------------------------------
def task_c_cox() -> pd.DataFrame:
    log.info("=" * 60)
    log.info("(c) Weighted Cox model for all-cause mortality")
    rows = []
    try:
        mort = read_mortality_fixed_width(MORT_FILE)
        nhanes = pd.read_csv(NHANES_CSV)
        nhanes["SEQN"] = nhanes["SEQN"].astype(int)
        demo = pd.read_sas(NHANES_DIR / "DEMO_J.xpt", format="xport")
        demo = demo[["SEQN", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]].copy()
        demo["SEQN"] = demo["SEQN"].astype(int)
        df = mort.merge(nhanes, on="SEQN").merge(demo, on="SEQN")

        df = df[df["eligstat"] == 1]
        df = df[df["age"] >= 18]
        df = df[df["permth_int"].notna() & (df["permth_int"] > 0)]
        df = df[df["crp"].notna() & (df["crp"] > 0)]
        df = df[df["mortstat"].isin([0, 1])]
        df["time_years"] = df["permth_int"] / 12.0
        df["event"] = df["mortstat"].astype(int)
        df["sex_female"] = (df["sex"] == "F").astype(int) if df["sex"].dtype == object \
            else (df["sex"] == 2).astype(int)

        cols = ["time_years", "event", "I", "age", "sex_female", "bmi",
                "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]
        sub = df[cols].dropna().copy()
        sub = sub[sub["WTMEC2YR"] > 0].copy()
        sub["I_z"] = (sub["I"] - sub["I"].mean()) / sub["I"].std(ddof=0)
        sub = add_psu_cluster(sub)
        covs = ["I_z", "age", "sex_female", "bmi"]
        log.info(f"  mortality analytic sample: n={len(sub)}, deaths={int(sub['event'].sum())}, "
                 f"mean follow-up={sub['time_years'].mean():.2f} y "
                 f"(R2 reference: n=4981, deaths=54 for logCRP+BMI model)")

        def cox_rows(use_weights: bool, label: str):
            out = []
            try:
                fit_cols = ["time_years", "event"] + covs + (["WTMEC2YR", "psu_cluster"] if use_weights else [])
                cph = CoxPHFitter(penalizer=0.0)
                kw = dict(duration_col="time_years", event_col="event")
                if use_weights:
                    kw.update(weights_col="WTMEC2YR", robust=True, cluster_col="psu_cluster")
                cph.fit(sub[fit_cols], **kw)
                summ = cph.summary
                for cov in covs:
                    b = float(summ.loc[cov, "coef"])
                    se = float(summ.loc[cov, "se(coef)"])
                    p = float(summ.loc[cov, "p"])
                    out.append({"model": label, "predictor": cov, "n": len(sub),
                                "n_events": int(sub["event"].sum()),
                                "beta": b, "se": se, "hr": float(np.exp(b)),
                                "hr_lower": float(np.exp(b - 1.96 * se)),
                                "hr_upper": float(np.exp(b + 1.96 * se)), "p": p})
                    consolidate("mortality_cox", "all_cause_mortality", cov, label, b, se, p,
                                len(sub), "hr", n_events=int(sub["event"].sum()),
                                weight="WTMEC2YR" if use_weights else "none")
            except Exception as e:
                log.error(f"  cox {label} failed: {e}")
            return out

        rows += cox_rows(False, "unweighted")
        rows += cox_rows(True, "weighted")
        for r in rows:
            if r["predictor"] == "I_z":
                log.info(f"  {r['model']}: HR(I_z)={r['hr']:.3f} ({r['hr_lower']:.3f}-{r['hr_upper']:.3f}), "
                         f"p={r['p']:.3g}")
    except Exception as e:
        log.error(f"(c) weighted Cox failed: {e}")
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "r5_mortality_cox_weighted.csv", index=False)
    log.info(f"  saved r5_mortality_cox_weighted.csv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------------------
# (d) weighted threshold (hinge) analysis with stratified PSU cluster bootstrap
# --------------------------------------------------------------------------------------
def grid_hinge_search(y, X_base, w, knots):
    best_k, best_ll = np.nan, -np.inf
    iz = X_base[:, 1]
    for k in knots:
        X = np.column_stack([X_base, np.maximum(iz - k, 0.0)])
        try:
            m = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=w).fit()
            if m.llf > best_ll:
                best_ll, best_k = m.llf, k
        except Exception:
            continue
    return best_k, best_ll


def task_d_threshold(ana: pd.DataFrame, n_boot: int = 200) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("(d) Weighted threshold (hinge) analysis with stratified PSU cluster bootstrap")
    knots = np.round(np.arange(-2.0, 2.0 + 1e-9, 0.05), 4)
    y = ana["allergy_composite"].values.astype(float)
    X_base = np.column_stack([
        np.ones(len(ana)), ana["I_z"].values, ana["age"].values,
        ana["sex_female"].values.astype(float), ana["bmi"].values,
    ])
    w = ana["WTMEC2YR"].values.astype(float)

    k_hat, ll_hat = grid_hinge_search(y, X_base, w, knots)
    log.info(f"  weighted best knot: {k_hat} (ll={ll_hat:.2f})")
    k_unw, ll_unw = grid_hinge_search(y, X_base, np.ones(len(ana)), knots)
    log.info(f"  unweighted best knot: {k_unw} (ll={ll_unw:.2f})")

    hinge_beta, hinge_p = np.nan, np.nan
    try:
        tmp = ana.copy()
        tmp["hinge"] = np.maximum(tmp["I_z"] - k_hat, 0.0)
        mf = fit_wglm(tmp, "allergy_composite ~ I_z + hinge + age + sex_female + bmi",
                      sm.families.Binomial(), "WTMEC2YR")
        hinge_beta, hinge_p = float(mf.params["hinge"]), float(mf.pvalues["hinge"])
        log.info(f"  hinge slope change at {k_hat}: beta={hinge_beta:.3f}, p={hinge_p:.3g}")
    except Exception as e:
        log.error(f"  final hinge model failed: {e}")

    # stratified cluster bootstrap: resample the 2 PSUs within each stratum (with replacement)
    strata = ana["SDMVSTRA"].values.astype(int)
    psu = ana["SDMVPSU"].values.astype(int)
    boot_knots, n_fail = [], 0
    idx_by_sp = {}
    for s in np.unique(strata):
        for p in (1, 2):
            idx_by_sp[(s, p)] = np.where((strata == s) & (psu == p))[0]
    for b in range(n_boot):
        parts = []
        for s in np.unique(strata):
            for _ in range(2):  # 2 PSU draws per stratum, with replacement
                p_draw = RNG.integers(1, 3)
                parts.append(idx_by_sp[(s, p_draw)])
        idx = np.concatenate(parts)
        kb, _ = grid_hinge_search(y[idx], X_base[idx], w[idx], knots)
        if np.isnan(kb):
            n_fail += 1
        else:
            boot_knots.append(kb)
        if (b + 1) % 50 == 0:
            log.info(f"  bootstrap {b + 1}/{n_boot} done (failures={n_fail})")
    boot_knots = np.array(boot_knots)
    ci_lo = float(np.percentile(boot_knots, 2.5)) if len(boot_knots) else np.nan
    ci_hi = float(np.percentile(boot_knots, 97.5)) if len(boot_knots) else np.nan
    log.info(f"  bootstrap knot 95% CI: [{ci_lo}, {ci_hi}] "
             f"(median={np.median(boot_knots):.2f}, n_ok={len(boot_knots)})")

    rows = [
        {"model": "weighted", "knot_hat": k_hat, "loglik": ll_hat,
         "ci_lower": ci_lo, "ci_upper": ci_hi,
         "boot_median": float(np.median(boot_knots)) if len(boot_knots) else np.nan,
         "boot_mean": float(np.mean(boot_knots)) if len(boot_knots) else np.nan,
         "n_boot_ok": len(boot_knots), "n_boot_fail": n_fail,
         "hinge_beta": hinge_beta, "hinge_p": hinge_p,
         "n": len(ana), "n_events": int(y.sum())},
        {"model": "unweighted", "knot_hat": k_unw, "loglik": ll_unw,
         "ci_lower": np.nan, "ci_upper": np.nan, "boot_median": np.nan,
         "boot_mean": np.nan, "n_boot_ok": 0, "n_boot_fail": 0,
         "hinge_beta": np.nan, "hinge_p": np.nan,
         "n": len(ana), "n_events": int(y.sum())},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "r5_threshold_weighted_bootstrap.csv", index=False)
    consolidate("threshold_hinge", "allergy_composite", "I_z hinge knot", "weighted",
                k_hat, (ci_hi - ci_lo) / (2 * 1.96) if not np.isnan(ci_lo) else np.nan,
                np.nan, len(ana), "beta", n_events=int(y.sum()), weight="WTMEC2YR",
                note=f"knot on I_z scale; bootstrap 95% CI [{ci_lo}, {ci_hi}]; "
                     f"stratified PSU resampling; hinge slope p={hinge_p:.3g}")
    log.info("  saved r5_threshold_weighted_bootstrap.csv")
    return out


def main():
    log.info("R5-07 NHANES survey-weighted analysis starting (seed=%d)", SEED)
    log.info("Method: statsmodels freq_weights + strata x PSU cluster-robust SE "
             "(fallback after samplics install timeout; design-consistent approximation)")
    weights = load_fixed_index_weights()
    ref_stats = j_cycle_component_stats()
    ana = load_j_analytic()

    task_a_main(ana)
    task_b_ige(weights, ref_stats)
    task_c_cox()
    task_d_threshold(ana, n_boot=200)

    cons = pd.DataFrame(CONSOLIDATED)
    cons.to_csv(TBL_DIR / "r5_nhanes_weighted_results.csv", index=False)
    log.info(f"saved consolidated r5_nhanes_weighted_results.csv ({len(cons)} rows)")
    log.info("R5-07 complete.")


if __name__ == "__main__":
    main()
