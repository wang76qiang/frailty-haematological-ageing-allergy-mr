#!/usr/bin/env python3
"""
R5-12: Cell-type-specific MR matrix (DICE cis-eQTL -> FinnGen allergy outcomes),
with OneK1K same-direction replication.

Outputs
-------
results/r5/tables/r5_celltype_mr_matrix.csv
results/r5/tables/r5_celltype_replication.csv
results/r5/figures/r5_celltype_heatmap.png
results/r5/logs/r5_12_summary.json
"""
import os
import re
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "r1"))
from utils import (load_dice, load_onek1k, harmonise_pair, mr_wald_ratio,
                   save_fig, FINNGEN_FULL_OUTCOMES)  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE, "data", "real")
OUT_TABLES = os.path.join(BASE, "results", "r5", "tables")
OUT_FIGS = os.path.join(BASE, "results", "r5", "figures")
OUT_LOGS = os.path.join(BASE, "results", "r5", "logs")
for d in (OUT_TABLES, OUT_FIGS, OUT_LOGS):
    os.makedirs(d, exist_ok=True)

SEED = 42
F_THRESH = 10.0

GENES = ["IL4", "IL5", "IL13", "IL33", "TSLP", "GATA3", "STAT6", "FOXP3",
         "RORC", "IL6", "TNF", "IL1R1", "STAT5A"]
NEG_CONTROLS = ["ARMS2", "OCA2", "COL1A1"]
ALL_GENES = GENES + NEG_CONTROLS

OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]

# DICE cell files present on disk (TREG not downloaded -> 15 cells)
DICE_CELLS = ["B_CELL_NAIVE", "CD4_NAIVE", "CD4_STIM", "CD8_NAIVE", "CD8_STIM",
              "M2", "MONOCYTES", "NK", "TFH", "TH1", "TH17", "TH2", "THSTAR",
              "TREG_MEM", "TREG_NAIVE"]

# DICE -> OneK1K cell-type mapping for replication
DICE2ONEK = {
    "TH1": "cd4et", "TH2": "cd4et", "TH17": "cd4et", "TFH": "cd4et",
    "THSTAR": "cd4et", "CD4_STIM": "cd4et", "TREG_MEM": "cd4et",
    "CD4_NAIVE": "cd4nc", "TREG_NAIVE": "cd4nc",
    "CD8_NAIVE": "cd8et", "CD8_STIM": "cd8et",
    "MONOCYTES": "monoc", "M2": "monoc",
    "NK": "nk", "B_CELL_NAIVE": "bin",
}

LOG = {"steps": []}


def log(msg):
    print(msg, flush=True)
    LOG["steps"].append(msg)


def collect_lead_snps():
    """For each (cell, gene) take the cis lead SNP (min Pvalue) with F>10."""
    rows = []
    for cell in DICE_CELLS:
        path = os.path.join(DATA_DIR, "dice_eqtl", f"{cell}.vcf")
        if not os.path.exists(path):
            log(f"WARN missing DICE file for {cell}")
            continue
        try:
            df = load_dice(cell, ALL_GENES, path=path)
        except Exception as e:
            log(f"WARN load_dice({cell}) failed: {e}")
            continue
        for gene in ALL_GENES:
            sub = df[df["gene"] == gene] if not df.empty else pd.DataFrame()
            if sub.empty:
                rows.append({"cell": cell, "gene": gene, "snp": None,
                             "instrument": "no instrument"})
                continue
            lead = sub.sort_values("p").iloc[0]
            fstat = float(lead["fstat"]) if np.isfinite(lead["fstat"]) else 0.0
            # keep utils-expected exposure names: beta / se / p
            rows.append({
                "cell": cell, "gene": gene, "snp": lead["snp"],
                "beta": lead["beta"], "se": lead["se"],
                "p": lead["p"], "ea": lead["ea"], "oa": lead["oa"],
                "fstat": fstat,
                "instrument": "pass" if fstat > F_THRESH else "weak (F<=10)",
            })
        n_pass = sum(1 for r in rows if r["cell"] == cell and r["instrument"] == "pass")
        log(f"  DICE {cell}: {n_pass}/{len(ALL_GENES)} genes with F>{F_THRESH}")
    return pd.DataFrame(rows)


CACHE_DIR = os.path.join(OUT_LOGS, "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def _atomic_write_csv(df, path):
    tmp = path + ".tmp"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def load_finngen_by_rsid(outcome, rsids):
    """Stream a FinnGen R12 file once, keeping rows whose rsid tokens match.
    Caches the matched subset so re-runs are fast (disk-safe atomic write)."""
    cache = os.path.join(CACHE_DIR, f"finngen_{outcome}.csv")
    rsids = set(rsids)
    if os.path.exists(cache):
        out = pd.read_csv(cache, dtype={"chrom": str})
        out = out[out["snp"].isin(rsids)]
        log(f"  FinnGen {outcome}: {len(out)} SNPs from cache")
        return out
    path = FINNGEN_FULL_OUTCOMES[outcome]
    # memory-safe pre-filter: stream the gz in 64MB text blocks, keep lines
    # containing any rsid (substring), then exact rsid-token match in pandas.
    import gzip as _gzip
    pat = re.compile("|".join(re.escape(s) for s in sorted(rsids)))
    hits = []
    with _gzip.open(path, "rt", errors="replace") as f:
        buf = ""
        while True:
            block = f.read(64 * 1024 * 1024)
            if not block:
                break
            buf += block
            lines = buf.split("\n")
            buf = lines.pop()
            hits.extend(ln for ln in lines if pat.search(ln))
        if buf and pat.search(buf):
            hits.append(buf)
    if not hits:
        return pd.DataFrame()
    from io import StringIO
    fg_cols = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
               "pval", "mlogp", "beta", "sebeta", "af_alt",
               "af_alt_cases", "af_alt_controls"]
    raw = pd.read_csv(StringIO("\n".join(hits)), sep="\t", low_memory=False,
                      header=None, names=fg_cols)
    raw = raw.rename(columns={"#chrom": "chrom", "rsids": "snp", "pval": "p",
                              "sebeta": "se"})
    tokens = raw["snp"].astype(str).str.split(",")
    exact = tokens.apply(lambda ts: any(t in rsids for t in ts))
    out = raw[exact].copy()
    out["snp"] = tokens[exact].apply(lambda ts: next(t for t in ts if t in rsids))
    out = out.dropna(subset=["snp", "beta", "se", "af_alt"]).drop_duplicates("snp")
    _atomic_write_csv(out, cache)
    log(f"  FinnGen {outcome}: matched {len(out)}/{len(rsids)} lead SNPs (cached)")
    return out


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted q-values."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    qv = np.empty(m)
    qv[order] = adj
    q[ok] = qv
    return q


def run_mr(leads):
    """Harmonise each lead SNP with the 3 FinnGen outcomes -> Wald ratio."""
    inst = leads[leads["instrument"] == "pass"].copy()
    rsids = inst["snp"].unique().tolist()
    log(f"Instruments passing F>{F_THRESH}: {len(inst)} SNP-gene-cell pairs, "
        f"{len(rsids)} unique SNPs")
    rows = []
    for outcome in OUTCOMES:
        try:
            out_df = load_finngen_by_rsid(outcome, rsids)
        except Exception as e:
            log(f"WARN FinnGen {outcome} load failed: {e}")
            out_df = pd.DataFrame()
        if out_df.empty:
            continue
        harm = harmonise_pair(inst, out_df)
        if harm.empty:
            log(f"  {outcome}: no harmonised SNPs")
            continue
        for _, r in harm.iterrows():
            b, se, p = mr_wald_ratio(r["beta"], r["beta_outcome"], r["se_outcome"])
            rows.append({
                "cell": r["cell"], "gene": r["gene"], "snp": r["snp"],
                "ea": r["ea"], "oa": r["oa"],
                "outcome": outcome,
                "beta_exp": r["beta"], "se_exp": r["se"], "p_exp": r["p"],
                "fstat": r["fstat"],
                "beta_outcome": r["beta_outcome"], "se_outcome": r["se_outcome"],
                "p_outcome": r["p_outcome"],
                "wald_beta": b, "wald_se": se, "wald_p": p,
                "or": np.exp(b),
                "or_lower": np.exp(b - 1.96 * se), "or_upper": np.exp(b + 1.96 * se),
            })
        log(f"  {outcome}: harmonised {len(harm)} pairs")
    res = pd.DataFrame(rows)
    if res.empty:
        return res, inst
    # BH-FDR across the full gene x cell family, per outcome
    res["wald_q"] = np.nan
    for outcome in OUTCOMES:
        m = res["outcome"] == outcome
        res.loc[m, "wald_q"] = bh_fdr(res.loc[m, "wald_p"].values)
    res["significant"] = res["wald_q"] < 0.05
    n_sig = int(res["significant"].sum())
    log(f"MR matrix: {len(res)} tests, {n_sig} FDR<0.05")
    return res, inst


def _onek1k_norm(df):
    out = pd.DataFrame({
        "gene": df["GENE"], "snp": df["RSID"], "a1": df["A1"], "a2": df["A2"],
        "rho": pd.to_numeric(df["SPEARMANS_RHO"], errors="coerce"),
        "p": pd.to_numeric(df["P_VALUE"], errors="coerce"),
    }).dropna(subset=["snp", "rho"])
    return out


def load_onek1k_full(cell, genes, rsids=None):
    """Exact-SNP rows from the full OneK1K eqtl table.

    Memory-light block streaming (~16MB blocks, keeps only matched lines) so it
    works under RAM pressure; result cached with atomic write.
    """
    cache = os.path.join(CACHE_DIR, f"onek1k_full_{cell}.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache)
        log(f"  OneK1K {cell}: {len(df)} rows from cache")
        return df
    import gzip as _gzip
    path = os.path.join(DATA_DIR, "onek1k", f"{cell}_eqtl_table.tsv.gz")
    if not os.path.exists(path):
        return pd.DataFrame()
    gene_set = set(genes)
    rsid_set = set(rsids) if rsids is not None and len(rsids) else None
    pat = None
    if rsid_set:
        pat = re.compile("|".join(re.escape(s) for s in sorted(rsid_set)))
    def _parse_line(ln):
        p = ln.split("\t")
        if len(p) < 15 or p[4] not in gene_set:
            return None
        if rsid_set is not None and p[2] not in rsid_set:
            return None
        try:
            return {"gene": p[4], "snp": p[2], "a1": p[8], "a2": p[9],
                    "rho": float(p[12]), "p": float(p[14])}
        except (ValueError, IndexError):
            return None

    rows = []
    try:
        with _gzip.open(path, "rt", errors="replace") as f:
            buf = ""
            while True:
                block = f.read(64 * 1024 * 1024)
                if not block:
                    break
                buf += block
                cut = buf.rfind("\n")
                if cut == -1:
                    continue
                segment, buf = buf[:cut], buf[cut + 1:]
                if pat is None:
                    for ln in segment.split("\n"):
                        r = _parse_line(ln)
                        if r:
                            rows.append(r)
                else:
                    for m in pat.finditer(segment):
                        s = segment.rfind("\n", 0, m.start()) + 1
                        e = segment.find("\n", m.end())
                        if e == -1:
                            e = len(segment)
                        r = _parse_line(segment[s:e])
                        if r:
                            rows.append(r)
            if buf:
                r = _parse_line(buf)
                if r:
                    rows.append(r)
    except Exception as e:
        log(f"WARN OneK1K stream {cell} failed: {e}")
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    _atomic_write_csv(df, cache)
    log(f"  OneK1K {cell}: streamed {len(df)} exact-SNP rows (cached)")
    return df


def load_onek1k_esnp(cell, genes):
    """Small per-gene lead-eSNP table (FDR-significant eQTLs only)."""
    esnp = os.path.join(DATA_DIR, "onek1k", f"{cell}_esnp_table.tsv.gz")
    try:
        if not os.path.exists(esnp):
            return pd.DataFrame()
        df = pd.read_csv(esnp, sep="\t", low_memory=False)
        df = df[df["GENE"].isin(set(genes))]
        return _onek1k_norm(df) if not df.empty else pd.DataFrame()
    except Exception as e:
        log(f"WARN OneK1K esnp load {cell} failed: {e}")
        return pd.DataFrame()


def replicate_onek1k(mr_res):
    """Same-direction replication of FDR-significant pairs in OneK1K.

    Primary: exact same SNP (from the per-gene lead-eSNP table; full table if
    needed).  Secondary: the gene's lead eSNP in OneK1K as a gene-level
    direction proxy (flagged in `note`).  Signs are allele-aligned to the
    DICE effect allele (ALT)."""
    sig = mr_res[mr_res["significant"]] if not mr_res.empty else pd.DataFrame()
    if sig.empty:
        log("No significant pairs -> OneK1K replication skipped")
        return pd.DataFrame()
    pairs = sig[["cell", "gene", "snp", "ea", "oa", "wald_beta",
                 "outcome"]].drop_duplicates()
    need_cells = sorted({DICE2ONEK[c] for c in pairs["cell"] if c in DICE2ONEK})
    need_genes = sorted(pairs["gene"].unique())
    log(f"OneK1K replication: {len(pairs)} pairs, cells={need_cells}, "
        f"genes={need_genes}")
    rows = []
    for ocell in need_cells:
        full = load_onek1k_full(ocell, need_genes,
                                rsids=pairs["snp"].unique())
        esnp = load_onek1k_esnp(ocell, need_genes)
        for _, pr in pairs[pairs["cell"].map(DICE2ONEK) == ocell].iterrows():
            base = {**pr.to_dict(), "onek1k_cell": ocell}
            exact = full[(full["gene"] == pr["gene"]) & (full["snp"] == pr["snp"])] \
                if not full.empty else pd.DataFrame()
            if not exact.empty:
                h, note, exact_match = exact.iloc[0], "", True
            else:
                glead = esnp[esnp["gene"] == pr["gene"]] if not esnp.empty else pd.DataFrame()
                if glead.empty:
                    rows.append({**base, "onek1k_snp": None,
                                 "onek1k_beta": np.nan, "onek1k_p": np.nan,
                                 "allele_aligned": None, "same_direction": None,
                                 "exact_snp": False,
                                 "note": "SNP absent in OneK1K full table; "
                                         "gene not a significant eQTL"})
                    continue
                h, note, exact_match = glead.sort_values("p").iloc[0], \
                    "lead-eSNP proxy (gene-level)", False
            ea = str(pr.get("ea", "")).upper()
            a1, a2 = str(h["a1"]).upper(), str(h["a2"]).upper()
            if ea and ea == a2:
                ok_beta, aligned = h["rho"], True
            elif ea and ea == a1:
                ok_beta, aligned = -h["rho"], True
            else:
                ok_beta, aligned = h["rho"], False
            same = (bool(np.sign(ok_beta) == np.sign(pr["wald_beta"]))
                    if aligned else None)
            rows.append({**base, "onek1k_snp": h["snp"],
                         "onek1k_beta": h["rho"], "onek1k_p": h["p"],
                         "onek1k_beta_aligned": ok_beta,
                         "allele_aligned": aligned, "same_direction": same,
                         "exact_snp": exact_match,
                         "note": note if aligned else note + " allele mismatch"})
    rep = pd.DataFrame(rows)
    valid = rep[(rep["same_direction"].notna()) & (rep["exact_snp"])] \
        if not rep.empty else rep
    if len(valid):
        n_same = int(valid["same_direction"].sum())
        p_sign = stats.binomtest(n_same, len(valid), 0.5).pvalue
        log(f"OneK1K exact-SNP same-direction: {n_same}/{len(valid)} "
            f"(binomial p={p_sign:.3e})")
    valid2 = rep[rep["same_direction"].notna()] if not rep.empty else rep
    if len(valid2):
        log(f"OneK1K same-direction incl. proxies: "
            f"{int(valid2['same_direction'].sum())}/{len(valid2)}")
    return rep


def plot_heatmap(mr_res, leads):
    import matplotlib.pyplot as plt
    import seaborn as sns
    cells = DICE_CELLS
    genes = GENES + NEG_CONTROLS
    fig, axes = plt.subplots(1, 3, figsize=(24, 7), sharey=True)
    for ax, outcome in zip(axes, OUTCOMES):
        sub = mr_res[mr_res["outcome"] == outcome] if not mr_res.empty else pd.DataFrame()
        piv = pd.DataFrame(np.nan, index=cells, columns=genes)
        annot = pd.DataFrame("", index=cells, columns=genes)
        for _, r in sub.iterrows():
            piv.loc[r["cell"], r["gene"]] = r["or"]
            annot.loc[r["cell"], r["gene"]] = (
                f"{r['or']:.2f}" + ("*" if r["significant"] else ""))
        # grey = no instrument
        for c in cells:
            for g in genes:
                inst = leads[(leads["cell"] == c) & (leads["gene"] == g)]
                if not inst.empty and inst.iloc[0]["instrument"] != "pass":
                    if pd.isna(piv.loc[c, g]):
                        annot.loc[c, g] = "—"
        cmap = sns.color_palette("RdBu_r", as_cmap=True)
        cmap.set_bad("#d9d9d9")
        sns.heatmap(piv, cmap=cmap, center=1.0, annot=annot, fmt="",
                    linewidths=0.4, ax=ax, cbar=outcome == OUTCOMES[-1],
                    cbar_kws={"label": "OR (Wald ratio)"}, vmin=0.25, vmax=4.0)
        ax.set_title(outcome.replace("_", " "))
        ax.set_xlabel("")
        ax.set_ylabel("DICE cell type" if outcome == OUTCOMES[0] else "")
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle("R5-12: Cell-type-specific MR (DICE eQTL → FinnGen). "
                 "* = BH-FDR<0.05 (gene×cell family); grey = no usable instrument",
                 fontsize=11)
    fig.tight_layout()
    save_fig(fig, os.path.join(OUT_FIGS, "r5_celltype_heatmap.png"))


def main():
    try:
        log("Collecting DICE lead SNPs ...")
        leads = collect_lead_snps()
        log("Running MR against FinnGen outcomes ...")
        mr_res, inst = run_mr(leads)

        # full matrix (incl. no-instrument rows) for the record
        if not mr_res.empty:
            full = leads.merge(
                mr_res.drop(columns=["beta_exp", "se_exp", "p_exp", "fstat",
                                     "ea", "oa"], errors="ignore"),
                on=["cell", "gene"], how="left",
                suffixes=("_lead", "_mr"))
        else:
            full = leads.copy()
        full.to_csv(os.path.join(OUT_TABLES, "r5_celltype_mr_matrix.csv"),
                    index=False)

        log("OneK1K replication ...")
        rep = replicate_onek1k(mr_res)
        if not rep.empty:
            rep.to_csv(os.path.join(OUT_TABLES, "r5_celltype_replication.csv"),
                       index=False)
        else:
            pd.DataFrame(columns=["cell", "gene", "snp", "outcome", "wald_beta",
                                  "onek1k_cell", "onek1k_beta", "same_direction",
                                  "note"]).to_csv(
                os.path.join(OUT_TABLES, "r5_celltype_replication.csv"), index=False)

        plot_heatmap(mr_res, leads)

        sig_pairs = []
        if not mr_res.empty:
            for _, r in mr_res[mr_res["significant"]].iterrows():
                sig_pairs.append(f"{r['cell']}×{r['gene']} ({r['outcome']}, "
                                 f"OR={r['or']:.2f}, q={r['wald_q']:.2e})")
        LOG["summary"] = {
            "n_instruments_pass": int((leads["instrument"] == "pass").sum()),
            "n_tests": int(len(mr_res)) if not mr_res.empty else 0,
            "n_significant": len(sig_pairs),
            "significant_pairs": sig_pairs,
            "n_replication_pairs": int(len(rep)) if not rep.empty else 0,
        }
        with open(os.path.join(OUT_LOGS, "r5_12_summary.json"), "w") as f:
            json.dump(LOG, f, indent=2, default=str)
        log("DONE R5-12")
    except Exception as e:
        import traceback
        log(f"FAILED R5-12: {e}")
        with open(os.path.join(OUT_LOGS, "r5_12_summary.json"), "w") as f:
            json.dump({"failed": str(e), "trace": traceback.format_exc(),
                       "steps": LOG["steps"]}, f, indent=2)
        raise


if __name__ == "__main__":
    main()
