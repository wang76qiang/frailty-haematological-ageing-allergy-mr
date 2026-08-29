#!/usr/bin/env python3
"""
R5-04 Winner's curse correction and effect-size credibility.

Three-sample design:
  * instrument SELECTION  : Astle 2016 blood GWAS (eo_N172275 / wbc_N172435),
                            p<5e-8 + LD clump (r2<0.001, 10 Mb);
  * exposure ESTIMATION   : Pan-UKBB component betas (independent of the
                            selection p-values -> no winner's-curse inflation
                            of the IVW denominator);
  * outcome               : FinnGen R12 allergy GWAS (independent cohort).

Both Astle and Pan-UKBB are GRCh37.  Pan-UKBB carries no rsids, so Astle rsids
are bridged to Pan-UKBB chr:pos:ref:alt keys through the 1000G EUR bim (also
GRCh37).  Allele letters are build-independent, so harmonisation is unaffected
by the build difference.

Compares: naive (Astle-selected + Astle-estimated) vs corrected (Astle-selected
+ Pan-UKBB-estimated) vs the R4 composite-index MR.  Methodological point
recorded in the log: with selection on the EXPOSURE and an independent outcome
GWAS, winner's curse inflates the exposure betas and therefore biases the IVW
ratio toward the null -- it cannot explain R4's large ORs; a corrected estimate
that remains strongly positive is the credibility check.
"""
import bisect
import os
import sys
import traceback

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "r1"))

import r5_common as C
import r1_utils

log = C.log

TRAIT_MAP = {
    "astle_eosinophil": {"component": "eosinophil_pct", "label": "Eosinophil %"},
    "astle_wbc": {"component": "wbc", "label": "WBC count"},
}
WC_CACHE = "r5_cache_panukbb_wc_corrected_grch37.csv"


def prethin_position(df: pd.DataFrame, window_bp: int = 2_000_000) -> pd.DataFrame:
    """Position-based pre-thinning so the downstream LD-clump matrix stays
    tractable (a 35k x 35k float32 r2 matrix OOM-kills the process).
    Keeps the lowest-p SNP per `window_bp` window per chromosome; the
    stringent LD clump (r2<0.001 / 10 Mb) is then applied on this thinned
    set. Conservative, acceptable for a sensitivity analysis."""
    if df is None or df.empty:
        return df
    kept = []
    for _, sub in df.groupby("chrom"):
        sub = sub.sort_values("p")
        kept_pos = []
        rows = []
        for _, r in sub.iterrows():
            pos = int(r["pos"])
            i = bisect.bisect_left(kept_pos, pos)
            clash = False
            if i < len(kept_pos) and kept_pos[i] - pos < window_bp:
                clash = True
            if not clash and i > 0 and pos - kept_pos[i - 1] < window_bp:
                clash = True
            if not clash:
                bisect.insort(kept_pos, pos)
                rows.append(r)
        if rows:
            kept.append(pd.DataFrame(rows))
    if not kept:
        return df.iloc[0:0].copy()
    return pd.concat(kept, ignore_index=True)


def select_and_clump(trait: str) -> tuple:
    df = C.select_astle_instruments(
        trait, p_thresh=5e-8, cache_name=f"r5_cache_astle_gws_{trait}.csv")
    df = prethin_position(df, window_bp=2_000_000)
    cl = C.clump_instruments(df, r2_thresh=0.001, kb=10000)
    return df, cl


def build_corrected(insts: dict) -> pd.DataFrame:
    """Astle instrument rsids -> 1KG-bim GRCh37 keys -> Pan-UKBB betas."""
    all_rsids = set()
    for df in insts.values():
        all_rsids.update(df["snp"].tolist())
    keys, k2r, n_bim = C.instrument_keys_grch37(
        pd.DataFrame({"snp": sorted(all_rsids)}))
    log(f"  {n_bim}/{len(all_rsids)} Astle instrument rsids have bim GRCh37 "
        f"coordinates -> querying Pan-UKBB ({len(keys)} keys)")
    comps = sorted({TRAIT_MAP[t]["component"] for t in insts})
    panu = C.scan_panukbb_components(keys, components=comps, workers=2,
                                     cache_name=WC_CACHE)
    panu = C.map_panukbb_to_rsid(panu, k2r)
    # guard: a stale/prebuilt cache must cover the current instrument set
    if not panu.empty:
        covered = set(panu["rsid"].unique())
        miss = [r for r in all_rsids if r not in covered]
        if len(miss) > 0.15 * len(all_rsids):
            log(f"  WARNING: cache covers {len(covered)}/{len(all_rsids)} rsids "
                f"-> forcing rescan")
            try:
                os.remove(os.path.join(C.CACHE_DIR, WC_CACHE))
            except OSError:
                pass
            panu = C.scan_panukbb_components(keys, components=comps, workers=2,
                                             cache_name=WC_CACHE)
            panu = C.map_panukbb_to_rsid(panu, k2r)
    log(f"  Pan-UKBB corrected betas recovered for "
        f"{panu['rsid'].nunique() if not panu.empty else 0} rsids")
    return panu


def ivw_fixed_row(har: pd.DataFrame, label: str) -> dict:
    bx = har["beta"].values.astype(float)
    by = har["beta_outcome"].values.astype(float)
    sy = har["se_outcome"].values.astype(float)
    k = len(har)
    if k == 0:
        return {f"{label}_n": 0}
    if k == 1:
        b, se, p = r1_utils.mr_wald_ratio(bx[0], by[0], sy[0])
    else:
        b, se, p = r1_utils.mr_ivw(bx, by, sy, random=False)
    return {f"{label}_n": k, f"{label}_beta": b, f"{label}_se": se,
            f"{label}_p": p, f"{label}_or": np.exp(b),
            f"{label}_or_lower": np.exp(b - 1.96 * se),
            f"{label}_or_upper": np.exp(b + 1.96 * se)}


def run_trait(trait: str, inst: pd.DataFrame, panu: pd.DataFrame,
              outcomes_fg: dict) -> list:
    """Naive (Astle) vs corrected (Pan-UKBB) IVW for one blood trait."""
    rows = []
    comp = TRAIT_MAP[trait]["component"]
    n_panu = C.PANUKBB_COMPONENTS[comp][1]
    n_astle = int(pd.to_numeric(inst["n"], errors="coerce").median()) \
        if "n" in inst.columns else C.N_ASTLE.get(trait, 172000)

    exp_naive = inst[["snp", "ea", "oa", "beta", "se"]].copy()
    sub = panu[panu["component"] == comp] if panu is not None and not panu.empty \
        else pd.DataFrame()
    exp_corr = pd.DataFrame()
    if not sub.empty:
        exp_corr = pd.DataFrame({
            "snp": sub["rsid"].astype(str),
            "ea": sub["alt"].astype(str),
            "oa": sub["ref"].astype(str),
            "beta": sub["beta"].astype(float),
            "se": sub["se"].astype(float),
        }).drop_duplicates("snp")

    for outcome, fg in outcomes_fg.items():
        row = {"trait": trait, "component": comp, "outcome": outcome}
        n_out = C.FINNGEN_N[outcome]
        out_df = fg[["snp", "ref", "alt", "beta", "se", "af_alt", "p"]].copy()
        har_n = C.harmonise_for_mr(exp_naive, out_df, n_astle, n_out,
                                   apply_steiger=False)[0]
        row.update(ivw_fixed_row(har_n, "naive"))
        if not exp_corr.empty:
            har_c = C.harmonise_for_mr(exp_corr, out_df, n_panu, n_out,
                                       apply_steiger=False)[0]
            row.update(ivw_fixed_row(har_c, "corrected"))
        else:
            row.update({"corrected_n": 0})
        rows.append(row)
    return rows


def add_r4_comparison(df: pd.DataFrame) -> pd.DataFrame:
    r4_path = os.path.join(C.R5_OUT, "..", "r4", "tables",
                           "r4_mr_diagnostics.csv")
    r4_path = os.path.abspath(r4_path)
    if not os.path.exists(r4_path):
        log(f"  R4 diagnostics not found: {r4_path}")
        return df
    r4 = pd.read_csv(r4_path)
    m = r4.drop_duplicates("outcome", keep="first").set_index("outcome")
    for i, row in df.iterrows():
        oc = row["outcome"]
        if oc not in m.index:
            continue
        r4r = m.loc[oc]
        r4_or = r4r.get("ivw_fixed_or")
        r4_beta = np.log(r4_or) if pd.notna(r4_or) and r4_or > 0 else np.nan
        df.at[i, "r4_index_or"] = r4_or
        df.at[i, "r4_index_beta"] = r4_beta
        df.at[i, "r4_index_p"] = r4r.get("ivw_fixed_p")
        for lab in ("naive", "corrected"):
            b = row.get(f"{lab}_beta")
            rb = r4_beta
            if pd.notna(b) and pd.notna(rb):
                df.at[i, f"{lab}_direction_concordant_with_r4"] = bool(
                    np.sign(b) == np.sign(rb))
        lo_n, hi_n = row.get("naive_or_lower"), row.get("naive_or_upper")
        lo_c, hi_c = row.get("corrected_or_lower"), row.get("corrected_or_upper")
        if all(pd.notna(x) for x in (lo_n, hi_n, lo_c, hi_c)):
            df.at[i, "naive_corrected_ci_overlap"] = bool(
                max(lo_n, lo_c) <= min(hi_n, hi_c))
    return df


def build_sensitivity(insts_raw: dict, insts: dict, panu: pd.DataFrame,
                      wc: pd.DataFrame) -> pd.DataFrame:
    """Winner's-curse shrinkage diagnostics per trait and per outcome."""
    rows = []
    for trait in TRAIT_MAP:
        comp = TRAIT_MAP[trait]["component"]
        raw, cl = insts_raw[trait], insts[trait]
        sub = panu[panu["component"] == comp] if panu is not None and not panu.empty \
            else pd.DataFrame()
        merged = cl.merge(sub[["rsid", "beta"]].rename(
            columns={"rsid": "snp", "beta": "beta_panu"}), on="snp", how="inner") \
            if not sub.empty else pd.DataFrame()
        ratio = np.nan
        if not merged.empty:
            # align corrected beta to the Astle effect allele (A1)
            a1 = merged["A1"].astype(str).str.upper() if "A1" in merged.columns else None
            alt = merged["alt"].astype(str).str.upper() if "alt" in merged.columns else None
            sign = np.where((a1 is not None) & (alt is not None) & (a1 == alt), 1.0, -1.0)
            aligned = merged["beta_panu"].values * sign
            ok = np.abs(aligned) > 1e-12
            if ok.any():
                ratio = float(np.median(np.abs(merged["beta"].values[ok]) /
                                        np.abs(aligned[ok])))
        rows.append({"trait": trait, "category": "instrument_selection",
                     "item": "counts",
                     "n_raw_gws": len(raw), "n_clumped": len(cl),
                     "median_abs_beta_naive": float(cl["beta"].abs().median()),
                     "median_abs_beta_corrected": float(np.abs(
                         merged["beta_panu"]).median()) if not merged.empty else np.nan,
                     "wc_inflation_ratio_naive_over_corrected": ratio})
    for _, r in wc.iterrows():
        nb, cb = r.get("naive_beta"), r.get("corrected_beta")
        rows.append({"trait": r["trait"], "category": "mr_effect",
                     "item": r["outcome"],
                     "naive_beta": nb, "corrected_beta": cb,
                     "naive_over_corrected_beta_ratio": (
                         nb / cb if pd.notna(nb) and pd.notna(cb) and abs(cb) > 1e-12
                         else np.nan),
                     "naive_p": r.get("naive_p"), "corrected_p": r.get("corrected_p"),
                     "naive_n": r.get("naive_n"), "corrected_n": r.get("corrected_n")})
    return pd.DataFrame(rows)


def main():
    log("=" * 72)
    log("R5-04 winner's curse correction (Astle-select / Pan-UKBB-estimate / "
        "FinnGen-outcome)")
    log("=" * 72)
    log("NOTE: Astle's UKBB subset partially overlaps Pan-UKBB, so this design\n"
        "      attenuates but does not fully eliminate selection-estimation\n"
        "      coupling.  With selection on the exposure and an independent\n"
        "      outcome, winner's curse inflates exposure betas and biases IVW\n"
        "      toward the null (conservative); it cannot manufacture R4's ORs.")

    log("\n[1/5] Selecting Astle instruments (p<5e-8, LD clump r2<0.001/10Mb)...")
    insts_raw, insts = {}, {}
    for trait in TRAIT_MAP:
        raw, cl = select_and_clump(trait)
        insts_raw[trait], insts[trait] = raw, cl
        log(f"  {trait}: {len(raw)} raw GWS SNPs -> {len(cl)} clumped instruments")

    log("\n[2/5] Building corrected (Pan-UKBB-estimated) exposures...")
    panu = build_corrected(insts)

    log("\n[3/5] Loading FinnGen outcome betas for instrument rsids...")
    all_rsids = set()
    for df in insts.values():
        all_rsids.update(df["snp"].tolist())
    outcomes_fg = {}
    for outcome in C.FINNGEN_OUTCOMES:
        fg = C.finngen_outcome_by_rsids(
            outcome, all_rsids,
            cache_name=f"r5_cache_finngen_outcome_wc_{outcome}.csv")
        outcomes_fg[outcome] = fg
        log(f"  {outcome}: {len(fg)} outcome SNPs")

    log("\n[4/5] Naive vs corrected IVW (trait x outcome)...")
    rows = []
    for trait in TRAIT_MAP:
        try:
            rows.extend(run_trait(trait, insts[trait], panu, outcomes_fg))
        except Exception as exc:
            log(f"  {trait} failed: {exc}")
            traceback.print_exc()
    wc = pd.DataFrame(rows)
    wc = add_r4_comparison(wc)
    p1 = os.path.join(C.R5_TABLES, "r5_winners_curse_corrected.csv")
    wc.to_csv(p1, index=False)
    log(f"Saved: {p1} ({len(wc)} rows)")
    for _, r in wc.iterrows():
        log(f"  {r['trait']} -> {r['outcome']}: "
            f"naive OR={r.get('naive_or', float('nan')):.3f} "
            f"(p={r.get('naive_p', float('nan')):.2g}, n={r.get('naive_n', 0)}), "
            f"corrected OR={r.get('corrected_or', float('nan')):.3f} "
            f"(p={r.get('corrected_p', float('nan')):.2g}, n={r.get('corrected_n', 0)}), "
            f"R4 index OR={r.get('r4_index_or', float('nan'))}")

    log("\n[5/5] Selection-sensitivity / shrinkage diagnostics...")
    sens = build_sensitivity(insts_raw, insts, panu, wc)
    p2 = os.path.join(C.R5_TABLES, "r5_astle_selection_sensitivity.csv")
    sens.to_csv(p2, index=False)
    log(f"Saved: {p2} ({len(sens)} rows)")

    log("\nR5-04 complete.")


if __name__ == "__main__":
    main()
