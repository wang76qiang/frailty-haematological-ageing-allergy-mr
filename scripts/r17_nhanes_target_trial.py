# -*- coding: utf-8 -*-
"""P0-2: NHANES (cycle J 2017-18) target-trial emulation for the blood-cell
type-2-like index and allergic outcomes.

Design (emulated trial): adults >=18 with complete index, outcome and
confounders; 'assignment' to high index (>= median). Propensity model:
age, sex, BMI, smoking status, race/ethnicity. Stabilized IPTW.
Analysis: survey-weighted logistic outcome ~ high index, cluster-robust SE
by SDMVPSU. E-value for the OR. Continuous per-SD association as sensitivity.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

NH = r"D:\衰老研究\v3_pipeline\data\real\nhanes"
OUT = r"D:\衰老研究\v3_pipeline\Nature子刊\tables"

merged = pd.read_csv(NH + r"\nhanes_merged.csv")

demo = pd.read_sas(NH + r"\DEMO_J.xpt", format="xport")[[
    "SEQN", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "RIDRETH1", "RIDAGEYR", "INDFMPIR"]]
smq = pd.read_sas(NH + r"\SMQ_J.xpt", format="xport")[["SEQN", "SMQ020", "SMQ040"]]

df = merged.merge(demo, on="SEQN", how="left").merge(smq, on="SEQN", how="left")
print("merged:", df.shape)

df = df[(df["age"] >= 18)].copy()
df["smoking"] = np.where(df["SMQ020"] == 1,
                         np.where(df["SMQ040"].isin([1, 2]), 1, 0),
                         0)
df.loc[df["SMQ020"].isna(), "smoking"] = np.nan
conf = ["age", "sex", "bmi", "smoking", "RIDRETH1", "INDFMPIR"]

med = df["I"].median()
df["high_index"] = (df["I"] >= med).astype(int)

results = []
# NB: in this extract allergy_composite is asthma-anchored and coincides with
# asthma_ever; hay_fever_year has too few observations (n=477) for a second arm.
for outcome in ["allergy_composite"]:
    d = df.dropna(subset=[outcome, "I", "WTMEC2YR"] + conf).copy()
    n = len(d)
    cases = int(d[outcome].sum())
    Xp = pd.get_dummies(d[conf].astype({"sex": str, "RIDRETH1": str}), drop_first=True).astype(float)
    Xp = sm.add_constant(Xp)
    ps_model = sm.Logit(d["high_index"], Xp).fit(disp=0)
    ps = ps_model.predict(Xp)
    pexp = d["high_index"].mean()
    d["sw"] = np.where(d["high_index"] == 1, pexp / ps, (1 - pexp) / (1 - ps))
    d["w"] = d["WTMEC2YR"] * d["sw"]
    Xo = sm.add_constant(d[["high_index"]].astype(float))
    m = sm.GLM(d[outcome], Xo, family=sm.families.Binomial(), freq_weights=d["w"])
    r = m.fit(cov_type="cluster", cov_kwds={"groups": d["SDMVPSU"]})
    ci = r.conf_int()
    orv = np.exp(r.params["high_index"])
    lo = np.exp(ci.loc["high_index", 0])
    hi = np.exp(ci.loc["high_index", 1])
    p = r.pvalues["high_index"]
    # E-values on the RR scale (VanderWeele & Ding): convert OR to RR at the
    # outcome prevalence because the outcome is not rare (15.6%)
    prev = d[outcome].mean()
    rr = orv / ((1 - prev) + prev * orv)
    rr_lo = lo / ((1 - prev) + prev * lo)
    ev = rr + np.sqrt(rr * (rr - 1)) if rr > 1 else np.nan
    ev_ci = rr_lo + np.sqrt(rr_lo * (rr_lo - 1)) if rr_lo > 1 else np.nan
    results.append(dict(outcome=outcome, n=n, cases=cases, estimand="high_index_vs_low_IPTW",
                        OR=orv, ci_lo=lo, ci_hi=hi, p=p, prevalence=prev, RR=rr,
                        E_value_point=ev, E_value_ci_limit=ev_ci))
    d["I_sd"] = (d["I"] - d["I"].mean()) / d["I"].std()
    Xc = sm.add_constant(d[["I_sd"]].astype(float))
    mc = sm.GLM(d[outcome], Xc, family=sm.families.Binomial(), freq_weights=d["w"])
    rc = mc.fit(cov_type="cluster", cov_kwds={"groups": d["SDMVPSU"]})
    cic = rc.conf_int()
    results.append(dict(outcome=outcome, n=n, cases=cases, estimand="per_SD_IPTW_sensitivity",
                        OR=np.exp(rc.params["I_sd"]), ci_lo=np.exp(cic.loc["I_sd", 0]),
                        ci_hi=np.exp(cic.loc["I_sd", 1]), p=rc.pvalues["I_sd"],
                        prevalence=np.nan, RR=np.nan,
                        E_value_point=np.nan, E_value_ci_limit=np.nan))

res = pd.DataFrame(results)
res.to_csv(OUT + r"\ST146_nhanes_target_trial.csv", index=False)
pd.set_option("display.width", 200)
print(res.to_string())
