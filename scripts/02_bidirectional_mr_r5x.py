#!/usr/bin/env python3
"""
R5-02 Bidirectional MR  (r5x implementation -- standard LD-clump variant).

Authoritative implementation used for the final R5-02 outputs.  Kept under a
distinct filename because another parallel agent also maintains
src/r5/02_bidirectional_mr.py (see results/r5/logs/R5_02_03_04_notes.md).

Differences vs the other implementation:
  * Clumping = standard greedy-by-p LD clump at r2 < 0.001 within 10 Mb
    (1KG EUR, per-chromosome windowed), NO proximity pre-thinning.  This keeps
    multiple independent instruments per 10 Mb window (e.g. 185 asthma
    instruments vs 13 with pre-thinning), preserving power.
  * No Steiger filtering on the reverse direction (instruments are selected at
    p < 5e-8 for the allergy exposure; Steiger on the blood-component outcome
    would discard valid instruments and shrink n further).
  * Harmonisation by rsid via r1_utils.harmonise_pair (same primitive as R4).

Reverse direction: FinnGen allergy outcomes (asthma / rhinitis / atopic
dermatitis) as exposures -> 6 Pan-UKBB blood components (+ Astle eosinophil
count and WBC as GRCh37 replication).  Instruments: p < 5e-8 (fallback
p < 1e-6 if fewer than 10 clumped SNPs, logged).

BUILD NOTE: Pan-UKBB summary stats are GRCh37 (verified: rs429358 @
19:45411941), FinnGen R12 is GRCh38.  Raw coordinate matching is impossible;
the bridge is FinnGen rsid -> 1KG-bim GRCh37 coordinate -> Pan-UKBB
chr:pos:ref:alt key (allele letters are build-independent).  Astle matches by
rsid only.

Forward direction is read from results/r4/tables/r4_mr_diagnostics.csv
(index -> allergy IVW) and combined with the reverse IVW into a bidirectional
matrix + heatmap figure.

Run:  .venv/Scripts/python.exe src/r5/02_bidirectional_mr_r5x.py  (project root)
"""

import os
import sys
import traceback

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def _find_project_root():
    """Locate the project root robustly (works from any launch location)."""
    candidates = [os.getcwd()]
    try:
        candidates.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    for cand in candidates:
        cur = os.path.abspath(cand)
        for _ in range(6):
            if os.path.isdir(os.path.join(cur, "src", "r5")) and \
               os.path.isdir(os.path.join(cur, "data")):
                return cur
            cur = os.path.dirname(cur)
    return os.getcwd()


_ROOT_HINT = _find_project_root()
sys.path.insert(0, os.path.join(_ROOT_HINT, "src", "r5"))
import r5x_utils as C

OUTCOMES = list(C.FINNGEN_OUTCOMES)
COMPONENTS = list(C.PANUKBB_COMPONENTS)
ASTLE_TRAITS = list(C.ASTLE_FILES)
P_PRIMARY = 5e-8
P_FALLBACK = 1e-6
MIN_SNPS = 10
CACHE = os.path.join(C.LOGS_DIR, "cache")

log_lines = []


def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)


def select_and_clump(outcome):
    """Return (clumped df, p_used, n_sig, note), with disk caching."""
    cache = os.path.join(CACHE, f"inst02_{outcome}.csv")
    meta = os.path.join(CACHE, f"inst02_{outcome}.meta.txt")
    if os.path.exists(cache) and os.path.exists(meta):
        parts = open(meta, encoding="utf-8").read().split("\t")
        return pd.read_csv(cache), float(parts[0]), int(parts[1]), parts[2]
    df = C.select_finngen_instruments(outcome, p_thresh=P_PRIMARY)
    n_sig = len(df)
    note = ""
    p_used = P_PRIMARY
    cl = C.clump_instruments(df, r2=0.001, kb=10000)
    if len(cl) < MIN_SNPS:
        n_primary = len(cl)
        df2 = C.select_finngen_instruments(outcome, p_thresh=P_FALLBACK)
        cl2 = C.clump_instruments(df2, r2=0.001, kb=10000)
        if len(cl2) > len(cl):
            cl = cl2
            p_used = P_FALLBACK
            note = (f"p relaxed to {P_FALLBACK}: only {n_primary} clumped SNPs at "
                    f"{P_PRIMARY}")
    os.makedirs(CACHE, exist_ok=True)
    cl.to_csv(cache, index=False)
    with open(meta, "w", encoding="utf-8") as fh:
        fh.write(f"{p_used}\t{n_sig}\t{note}")
    return cl, p_used, n_sig, note


def scan_panukbb_cached(keys):
    """Scan 6 Pan-UKBB files, caching per component (incremental, OOM-safe)."""
    os.makedirs(CACHE, exist_ok=True)
    pan = {}
    missing = []
    for comp in COMPONENTS:
        cache = os.path.join(CACHE, f"pan02_{comp}.csv")
        if os.path.exists(cache):
            pan[comp] = pd.read_csv(cache, dtype={"chrom": str})
        else:
            missing.append(comp)
    if missing:
        log(f"  scanning Pan-UKBB for: {missing}")
        fresh = C.retry_call(
            lambda: C.scan_panukbb_components(keys, components=missing,
                                              n_jobs=min(2, len(missing))),
            label="panukbb_scan")
        for comp in missing:
            df = fresh.get(comp, pd.DataFrame())
            if df is None:
                df = pd.DataFrame()
            df.to_csv(os.path.join(CACHE, f"pan02_{comp}.csv"), index=False)
            pan[comp] = df
    return pan


def main():
    log("=" * 72)
    log("R5-02 Bidirectional MR (reverse: allergy -> blood components) [r5x]")
    log("=" * 72)

    # ---- 1. instruments per outcome ----
    instruments = {}
    for out in OUTCOMES:
        try:
            cl, p_used, n_sig, note = select_and_clump(out)
            instruments[out] = (cl, p_used)
            log(f"{out}: {n_sig} SNPs at p<{P_PRIMARY}; clumped={len(cl)} "
                f"(p_used={p_used})" + (f"  NOTE: {note}" if note else ""))
        except Exception:
            log(f"{out}: FAILED instrument selection\n{traceback.format_exc()}")
            instruments[out] = (pd.DataFrame(), P_PRIMARY)

    valid = {o: cl for o, (cl, _) in instruments.items() if not cl.empty}
    if not valid:
        log("No instruments for any outcome; aborting.")
        return

    # ---- 2. union scans ----
    all_inst = pd.concat(valid.values(), ignore_index=True).drop_duplicates("snp")
    keys, key2rsid = C.panukbb_keys_via_bim(all_inst)
    log(f"Union instrument SNPs: {len(all_inst)}; Pan-UKBB keys (GRCh37): {len(keys)}")

    try:
        log("Scanning Pan-UKBB component files ...")
        pan = scan_panukbb_cached(keys)
        for comp, df in pan.items():
            log(f"  {comp}: {len(df)} matched rows")
    except Exception:
        log(f"Pan-UKBB scan FAILED\n{traceback.format_exc()}")
        pan = {c: pd.DataFrame() for c in COMPONENTS}

    try:
        log("Scanning Astle eo/wbc files for union rsids ...")
        astle = C.scan_astle_rsids(all_inst["snp"].tolist())
        for t, df in astle.items():
            log(f"  {t}: {len(df)} matched rows")
    except Exception:
        log(f"Astle scan FAILED\n{traceback.format_exc()}")
        astle = {t: pd.DataFrame() for t in ASTLE_TRAITS}

    # ---- 3. reverse MR per outcome x trait ----
    rows = []
    for out in OUTCOMES:
        cl, p_used = instruments[out]
        if cl.empty:
            continue
        exp = cl[["snp", "ea", "oa", "beta", "se", "p"]].copy()
        for comp in COMPONENTS:
            try:
                df = pan.get(comp, pd.DataFrame())
                if df is None or df.empty:
                    log(f"{out} x {comp}: no Pan-UKBB matches")
                    continue
                mapped = C.attach_rsid_via_keys(df, key2rsid)
                if mapped.empty:
                    log(f"{out} x {comp}: 0 SNPs mapped via GRCh37 bridge")
                    continue
                out_frame = C.panukbb_outcome_frame(mapped)
                har = C.r1_utils.harmonise_pair(exp, out_frame)
                if har.empty or len(har) < 3:
                    log(f"{out} x {comp}: <3 harmonised SNPs (n={len(har)})")
                    continue
                for m in C.mr_battery(har):
                    m.update({"outcome": out, "trait": comp, "source": "panukbb",
                              "n_instruments": len(cl), "p_threshold": p_used})
                    rows.append(m)
                log(f"{out} x {comp}: harmonised n={len(har)}")
            except Exception:
                log(f"{out} x {comp}: FAILED\n{traceback.format_exc()}")
        for trait in ASTLE_TRAITS:
            try:
                df = astle.get(trait, pd.DataFrame())
                if df is None or df.empty:
                    log(f"{out} x {trait}: no Astle matches")
                    continue
                out_frame = df[df["snp"].isin(set(cl["snp"]))].copy()
                har = C.r1_utils.harmonise_pair(exp, out_frame)
                if har.empty or len(har) < 3:
                    log(f"{out} x {trait}: <3 harmonised SNPs (n={len(har)})")
                    continue
                for m in C.mr_battery(har):
                    m.update({"outcome": out, "trait": trait, "source": "astle",
                              "n_instruments": len(cl), "p_threshold": p_used})
                    rows.append(m)
                log(f"{out} x {trait}: harmonised n={len(har)}")
            except Exception:
                log(f"{out} x {trait}: FAILED\n{traceback.format_exc()}")

    res = pd.DataFrame(rows)
    if not res.empty:
        cols = ["outcome", "trait", "source", "method", "n", "n_instruments",
                "p_threshold", "beta", "se", "p", "or", "or_lower", "or_upper",
                "egger_intercept", "egger_intercept_p"]
        cols = [c for c in cols if c in res.columns]
        res = res[cols]
        path = os.path.join(C.TABLES_DIR, "r5_reverse_mr_components.csv")
        res.to_csv(path, index=False)
        log(f"Saved {path} ({len(res)} rows)")

    # ---- 4. bidirectional matrix (forward from R4, reverse from this task) ----
    matrix_rows = []
    try:
        r4 = pd.read_csv(os.path.join(C.ROOT, "results", "r4", "tables", "r4_mr_diagnostics.csv"))
        for _, r in r4.iterrows():
            matrix_rows.append({
                "outcome": r["outcome"], "component": "immunosenescence_index",
                "direction": "forward_index_to_allergy",
                "beta": np.log(r["ivw_fixed_or"]), "se": np.nan,
                "or": r["ivw_fixed_or"], "or_lower": np.nan, "or_upper": np.nan,
                "p": r["ivw_fixed_p"], "n": r.get("n_snps_harmonised", np.nan),
                "source": "R4_index_IVW",
            })
    except Exception:
        log(f"Forward (R4) read FAILED\n{traceback.format_exc()}")
    if not res.empty:
        ivw = res[res["method"] == "IVW_fixed"]
        for _, r in ivw.iterrows():
            matrix_rows.append({
                "outcome": r["outcome"], "component": r["trait"],
                "direction": "reverse_allergy_to_component",
                "beta": r["beta"], "se": r["se"], "or": r["or"],
                "or_lower": r["or_lower"], "or_upper": r["or_upper"],
                "p": r["p"], "n": r["n"], "source": r["source"],
            })
    mat = pd.DataFrame(matrix_rows)
    if not mat.empty:
        path = os.path.join(C.TABLES_DIR, "r5_bidirectional_matrix.csv")
        mat.to_csv(path, index=False)
        log(f"Saved {path} ({len(mat)} rows)")

    # ---- 5. figure ----
    try:
        make_figure(mat)
    except Exception:
        log(f"Figure FAILED\n{traceback.format_exc()}")

    with open(os.path.join(C.LOGS_DIR, "r5_02_run.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    log("R5-02 complete.")


def make_figure(mat):
    if mat.empty:
        return
    fwd = mat[mat["direction"] == "forward_index_to_allergy"].set_index("outcome")
    rev = mat[mat["direction"] == "reverse_allergy_to_component"]
    trait_order = COMPONENTS + ASTLE_TRAITS
    trait_order = [t for t in trait_order if t in rev["component"].unique()]

    beta_piv = rev.pivot_table(index="outcome", columns="component", values="beta")
    beta_piv = beta_piv.reindex(index=OUTCOMES, columns=trait_order)
    p_piv = rev.pivot_table(index="outcome", columns="component", values="p")
    p_piv = p_piv.reindex(index=OUTCOMES, columns=trait_order)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5),
                             gridspec_kw={"width_ratios": [1, 4]})

    lf = fwd.reindex(OUTCOMES)
    vals = lf["beta"].values.reshape(-1, 1)
    annot = []
    for o in OUTCOMES:
        p = lf.loc[o, "p"] if o in lf.index else np.nan
        star = "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else "ns"))
        b = lf.loc[o, "beta"] if o in lf.index else np.nan
        annot.append(f"{b:.2f}\n{star}" if np.isfinite(b) else "NA")
    sns.heatmap(vals, annot=np.array(annot).reshape(-1, 1), fmt="", cmap="RdBu_r",
                center=0, ax=axes[0], yticklabels=OUTCOMES,
                xticklabels=["index -> allergy\nlog(OR)"], cbar_kws={"label": "log(OR)"})
    axes[0].set_title("Forward (R4 IVW)")

    annot2 = []
    for o in OUTCOMES:
        line = []
        for t in trait_order:
            b = beta_piv.loc[o, t] if (o in beta_piv.index and t in beta_piv.columns) else np.nan
            p = p_piv.loc[o, t] if (o in p_piv.index and t in p_piv.columns) else np.nan
            star = "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else ""))
            line.append(f"{b:+.3f}{star}" if np.isfinite(b) else "NA")
        annot2.append(line)
    sns.heatmap(beta_piv.values, annot=np.array(annot2), fmt="", cmap="RdBu_r",
                center=0, ax=axes[1], yticklabels=OUTCOMES, xticklabels=trait_order,
                cbar_kws={"label": "beta (allergy -> component)"})
    axes[1].set_title("Reverse MR (R5, IVW fixed)")
    axes[1].tick_params(axis="x", rotation=35)
    fig.suptitle("Bidirectional MR: immunosenescence index <-> allergy", fontsize=12)
    fig.tight_layout()
    out = os.path.join(C.FIGURES_DIR, "r5_bidirectional_diagram.png")
    C.r1_utils.save_fig(fig, out)
    log(f"Saved {out}")


if __name__ == "__main__":
    main()
