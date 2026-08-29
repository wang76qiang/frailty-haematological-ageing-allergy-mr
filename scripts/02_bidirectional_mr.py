#!/usr/bin/env python3
"""
R5-02 Bidirectional MR: allergy <-> blood components.

Reverse direction: FinnGen allergy outcomes (asthma / rhinitis / atopic
dermatitis) as exposures -> 6 Pan-UKBB blood components, with Astle 2016
eosinophil-count and WBC GWAS as independent replication outcomes.
Instruments: p<5e-8 genome-wide significant allergy SNPs, LD-clumped
(r2<0.001, 10 Mb, 1KG EUR); falls back to p<1e-6 if <10 clumped SNPs (logged).

BUILD NOTE: Pan-UKBB summary stats are GRCh37 (verified: rs429358 @
19:45411941), FinnGen R12 is GRCh38.  Raw coordinate matching is impossible;
the bridge is FinnGen rsid -> 1KG-bim GRCh37 coordinate -> Pan-UKBB
chr:pos:ref:alt key (allele letters are build-independent).  Astle matches by
rsid directly.

Forward direction is read from results/r4/tables/r4_mr_diagnostics.csv
(index -> allergy IVW) and combined with the reverse IVW into a bidirectional
matrix + heatmap figure.
"""
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

log = C.log
OUTCOMES = list(C.FINNGEN_OUTCOMES)
COMPONENTS = list(C.PANUKBB_COMPONENTS)
ASTLE_TRAITS = list(C.ASTLE_TRAITS)
MIN_SNPS = 10


def build_instruments() -> dict:
    insts = {}
    for outcome in OUTCOMES:
        gwas = C.select_finngen_instruments(outcome, p_thresh=5e-8)
        inst = C.clump_instruments(gwas, r2_thresh=0.001, kb=10000)
        if len(inst) < MIN_SNPS:
            log(f"  {outcome}: only {len(inst)} SNPs at 5e-8 -> fallback p<1e-6")
            gwas = C.select_finngen_instruments(outcome, p_thresh=1e-6)
            inst = C.clump_instruments(gwas, r2_thresh=0.001, kb=10000)
        insts[outcome] = inst
        log(f"  {outcome}: {len(inst)} clumped instruments")
    return insts


def build_panukbb_outcomes(insts: dict):
    all_rsids = set()
    for df in insts.values():
        all_rsids.update(df["snp"].tolist())
    keys, k2r, n_bim = C.instrument_keys_grch37(
        pd.DataFrame({"snp": sorted(all_rsids)}))
    log(f"  {n_bim}/{len(all_rsids)} instrument rsids mapped to GRCh37 keys "
        f"({len(keys)} keys)")
    cname = "r5_cache_panukbb_reverse_mr_grch37.csv"
    panu = C.scan_panukbb_components(keys, components=COMPONENTS, workers=2,
                                     cache_name=cname)
    panu = C.map_panukbb_to_rsid(panu, k2r)
    n_rec = panu["rsid"].nunique() if not panu.empty else 0
    # coverage guard: a partial/stale cache must cover >=85% of mapped rsids
    if n_rec < 0.85 * n_bim:
        log(f"  coverage {n_rec}/{n_bim} < 85% -> forcing full rescan")
        try:
            os.remove(os.path.join(C.CACHE_DIR, cname))
        except OSError:
            pass
        panu = C.scan_panukbb_components(keys, components=COMPONENTS, workers=2,
                                         cache_name=cname)
        panu = C.map_panukbb_to_rsid(panu, k2r)
        n_rec = panu["rsid"].nunique() if not panu.empty else 0
    log(f"  Pan-UKBB betas recovered for {n_rec} rsids")
    return panu


def build_astle_outcomes(insts: dict):
    all_rsids = set()
    for df in insts.values():
        all_rsids.update(df["snp"].tolist())
    out = {}
    for trait in ASTLE_TRAITS:
        df = C.scan_astle_rsids(trait, all_rsids,
                                cache_name=f"r5_cache_astle_outcome_{trait}.csv")
        out[trait] = df
        log(f"  {trait}: {len(df)} rsids recovered")
    return out


def reverse_mr(insts: dict, panu: pd.DataFrame, astle: dict) -> list:
    rows = []
    for outcome in OUTCOMES:
        inst = insts[outcome]
        exp = inst[["snp", "ea", "oa", "beta", "se"]].copy()
        n_exp = C.FINNGEN_N[outcome]
        targets = [(c, "panukbb") for c in COMPONENTS] + [(t, "astle") for t in ASTLE_TRAITS]
        for trait, source in targets:
            if source == "panukbb":
                sub = panu[panu["component"] == trait] if not panu.empty else pd.DataFrame()
                if sub.empty:
                    continue
                out_df = pd.DataFrame({
                    "snp": sub["rsid"].astype(str),
                    "ref": sub["ref"].astype(str),
                    "alt": sub["alt"].astype(str),
                    "beta": sub["beta"].astype(float),
                    "se": sub["se"].astype(float),
                    "af_alt": sub["af_alt"].astype(float),
                    "p": sub["p"].astype(float),
                })
                n_out = C.PANUKBB_COMPONENTS[trait][1]
            else:
                sub = astle.get(trait, pd.DataFrame())
                if sub is None or sub.empty:
                    continue
                out_df = pd.DataFrame({
                    "snp": sub["snp"].astype(str),
                    "ref": sub["ref"].astype(str),
                    "alt": sub["alt"].astype(str),
                    "beta": sub["beta"].astype(float),
                    "se": sub["se"].astype(float),
                    "af_alt": sub["af_alt"].astype(float),
                    "p": sub["p"].astype(float),
                })
                n_out = int(pd.to_numeric(sub["n"], errors="coerce").median()) \
                    if "n" in sub.columns else 172000
            try:
                har, diag = C.harmonise_for_mr(exp, out_df, n_exp, n_out,
                                               apply_steiger=True)
                if har.empty:
                    log(f"  {outcome} -> {trait}: 0 harmonised")
                    continue
                battery = C.mr_battery(har, do_presso=True)
                for m in battery:
                    m.update({"outcome": outcome, "trait": trait, "source": source,
                              "n_harmonised": diag["n_harmonised"],
                              "steiger_applied": diag["steiger_applied"]})
                    rows.append(m)
                ivw = next((r for r in battery if r["method"] == "IVW_fixed"),
                           battery[0])
                log(f"  {outcome} -> {trait} ({source}): n={ivw['n_snps']} "
                    f"beta={ivw['beta']:.4f} p={ivw['p']:.3g}")
            except Exception as exc:
                log(f"  {outcome} -> {trait} failed: {exc}")
                traceback.print_exc()
    return rows


def build_matrix(comp_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # forward (R4 index -> allergy)
    r4_path = os.path.abspath(os.path.join(C.R5_OUT, "..", "r4", "tables",
                                           "r4_mr_diagnostics.csv"))
    if os.path.exists(r4_path):
        r4 = pd.read_csv(r4_path).drop_duplicates("outcome").set_index("outcome")
        for outcome in OUTCOMES:
            if outcome not in r4.index:
                continue
            r = r4.loc[outcome]
            orv = r.get("ivw_fixed_or")
            rows.append({"direction": "forward_index_to_allergy", "outcome": outcome,
                         "trait": "immune_index", "method": "IVW_fixed",
                         "n_snps": r.get("n_snps_harmonised"),
                         "beta": np.log(orv) if pd.notna(orv) and orv > 0 else np.nan,
                         "se": np.nan, "p": r.get("ivw_fixed_p"),
                         "or": orv, "or_lower": np.nan, "or_upper": np.nan})
    ivw = comp_df[comp_df["method"] == "IVW_fixed"]
    for _, r in ivw.iterrows():
        rows.append({"direction": "reverse_allergy_to_component",
                     "outcome": r["outcome"], "trait": r["trait"],
                     "method": "IVW_fixed", "n_snps": r["n_snps"],
                     "beta": r["beta"], "se": r["se"], "p": r["p"],
                     "or": r["or"], "or_lower": r["or_lower"], "or_upper": r["or_upper"]})
    return pd.DataFrame(rows)


def make_figure(matrix: pd.DataFrame):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    rev = matrix[matrix["direction"] == "reverse_allergy_to_component"].copy()
    piv = rev.pivot_table(index="trait", columns="outcome", values="beta")
    sig = rev.pivot_table(index="trait", columns="outcome", values="p")
    annot = piv.copy().astype(str)
    for t in piv.index:
        for o in piv.columns:
            b = piv.loc[t, o]
            p = sig.loc[t, o]
            if pd.isna(b):
                annot.loc[t, o] = ""
            else:
                star = "*" if p < 0.05 else ""
                annot.loc[t, o] = f"{b:+.3f}{star}"
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(piv, annot=annot, fmt="", cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "reverse-MR beta (allergy -> component)"})
    ax.set_title("R5-02 Reverse MR: genetically predicted allergy -> blood components\n"
                 "(* p<0.05, IVW fixed; Astle = independent replication)")
    fig.tight_layout()
    out = os.path.join(C.R5_FIGURES, "r5_bidirectional_diagram.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved figure: {out}")


def main():
    log("=" * 72)
    log("R5-02 bidirectional MR (allergy <-> blood components)")
    log("=" * 72)
    log("\n[1/5] Building allergy instruments (p<5e-8, clump r2<0.001/10Mb)...")
    insts = build_instruments()
    log("\n[2/5] Scanning Pan-UKBB components (GRCh37 bridge via 1KG bim)...")
    panu = build_panukbb_outcomes(insts)
    log("\n[3/5] Extracting Astle replication outcomes (rsid match)...")
    astle = build_astle_outcomes(insts)
    log("\n[4/5] Reverse MR battery (3 outcomes x 8 targets x 5 methods)...")
    rows = reverse_mr(insts, panu, astle)
    comp = pd.DataFrame(rows)
    p1 = os.path.join(C.R5_TABLES, "r5_reverse_mr_components.csv")
    comp.to_csv(p1, index=False)
    log(f"Saved: {p1} ({len(comp)} rows)")
    log("\n[5/5] Bidirectional matrix + heatmap...")
    matrix = build_matrix(comp)
    p2 = os.path.join(C.R5_TABLES, "r5_bidirectional_matrix.csv")
    matrix.to_csv(p2, index=False)
    log(f"Saved: {p2} ({len(matrix)} rows)")
    make_figure(matrix)
    log("\nR5-02 complete.")


if __name__ == "__main__":
    main()
