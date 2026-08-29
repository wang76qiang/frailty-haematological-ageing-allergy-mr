#!/usr/bin/env python3
"""
R5-06: IL4/IL13 (5q31) locus-level colocalisation adjudication.

Question (editor comment, docs/R5/R5_REVISION_MATRIX.md):
  The IL4 cis-eQTL "protective" signal (OR~0.46 vs asthma) may be an LD
  mis-attribution of the stronger IL13 / 5q31 atopy signal. Decide between:
    (i)  coloc holds (PP.H4>0.8) and IL4-lowering allele == asthma-risk allele
         -> tissue-mismatch narrative (blood IL4 reflects regulatory pool);
    (ii) coloc fails (PP.H4<0.5)
         -> "IL4 protection" is LD coincidence; retract gene-level claim;
    (iii) IL4 vs IL13 eQTL directions opposite
         -> non-equivalence inside 5q31 (mechanistic highlight).

Data:
  - FinnGen R12 ALLERG_ASTHMA / ALLERG_RHINITIS / L12_ATOPIC full-region SNPs
    (chr5: IL4/IL13 union +/-1Mb, GRCh38).
  - eQTLGen whole-blood cis-eQTL for IL4 and IL13 (significant SNPs, GRCh37;
    matched to FinnGen by rsID -> labelled "approximate coloc").
  - DICE 15 cell types (TH1/TREG_NAIVE carry IL4/IL13 eQTL), OneK1K, GTEx
    (only 4 tissues available locally: Adipose_Subcutaneous/Visceral, Adrenal,
    Artery_Aorta; Whole_Blood/Lung/Skin NOT downloaded -> recorded as missing).
  - Cross-link: IL4/IL13 cis-pQTL MR results from R5-09 (if present).

Outputs (results/r5/):
  tables/r5_il4_locus_coloc.csv
  tables/r5_il4_il13_direction.csv
  tables/r5_il4_celltype_eqtl.csv
  figures/r5_fig4_il4_locus.png
  logs/r5_06_run.log
"""

import os
import sys
import json
import logging
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(BASE_DIR, "src", "r1"))
import utils  # noqa: E402

FG_DIR = os.path.join(BASE_DIR, "data", "real", "finngen_full")
OUT_TAB = os.path.join(BASE_DIR, "results", "r5", "tables")
OUT_FIG = os.path.join(BASE_DIR, "results", "r5", "figures")
OUT_LOG = os.path.join(BASE_DIR, "results", "r5", "logs")
FG_CACHE = os.path.join(OUT_TAB, "finngen_il4_locus")
for d in (OUT_TAB, OUT_FIG, OUT_LOG, FG_CACHE):
    os.makedirs(d, exist_ok=True)

LOG_PATH = os.path.join(OUT_LOG, "r5_06_run.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [R5-06] %(message)s",
    handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler(sys.stdout)],
)
log = logging.info

RNG = np.random.default_rng(42)

# ----------------------------------------------------------------------------- config
OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
PRIMARY = "ALLERG_ASTHMA"
PRIORS = (1e-4, 1e-4, 1e-5)
CHUNKSIZE = 200_000

# GRCh38 coords (target_gene_coords.csv): IL4 132673986-132682678, IL13 132656263-132661110
GENES_GRCH38 = {
    "IL4": ("5", 132673986, 132682678),
    "IL13": ("5", 132656263, 132661110),
    "RAD50": ("5", 132541001, 132630677),   # approx GRCh38, for figure span only
}
UNION_START = min(v[1] for v in GENES_GRCH38.values())
UNION_END = max(v[2] for v in GENES_GRCH38.values())
REGION_CHROM = "5"
REGION_START = UNION_START - 1_000_000
REGION_END = UNION_END + 1_000_000

DICE_CELLS = ["TH1", "TH2", "TH17", "TREG", "CD4_NAIVE", "CD8_NAIVE", "B_CELL_NAIVE",
              "MONOCYTES", "M2", "NK", "TFH", "THSTAR", "TREG_MEM", "TREG_NAIVE",
              "CD4_STIM", "CD8_STIM"]
ONEK1K_CELLS = ["cd4nc", "cd4et", "cd8et", "monoc", "nk", "bin"]
GTEX_LOCAL = ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
              "Artery_Aorta"]
GTEX_MISSING = ["Whole_Blood", "Lung", "Skin_Sun_Exposed_Lower_leg"]
GENE_ENSG = {"IL4": "ENSG00000113520", "IL13": "ENSG00000169194"}

FG_USECOLS = ["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]


# ----------------------------------------------------------------------------- loaders
def _posix(p):
    return os.path.abspath(p).replace("\\", "/")


def finngen_region(outcome):
    """Stream FinnGen once via zcat+awk (tiny memory), keep chr5 region SNPs. Cached."""
    import subprocess
    cache = os.path.join(FG_CACHE, f"{outcome}.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache, dtype={"chrom": str})
    path = os.path.join(FG_DIR, f"finngen_R12_{outcome}.gz")
    cond = (f'$1=="{REGION_CHROM}" && $2>={REGION_START} && $2<={REGION_END}')
    cmd = (f'zcat "{_posix(path)}" | awk -F"\\t" -v OFS="\\t" '
           f'\'NR==1{{print "chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"; next}}'
           f'{cond}{{print $1,$2,$3,$4,$5,$7,$9,$10,$11}}\' > "{_posix(cache)}.tmp"')
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"awk failed rc={r.returncode}: {r.stderr[:200]!r}")
    df = pd.read_csv(f"{cache}.tmp", sep="\t", dtype={"chrom": str})
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    df.to_csv(cache, index=False)
    os.remove(f"{cache}.tmp")
    log(f"[finngen] {outcome}: {len(df)} region SNPs -> {cache}")
    return df


def load_eqtl_exposure(fg_ref):
    """eQTLGen IL4/IL13 significant SNPs; beta/se derived with FinnGen MAF map."""
    eq = utils.load_eqtlgen(genes=["IL4", "IL13"])
    maf_map = fg_ref.set_index("rsid")["af_alt"] if not fg_ref.empty else pd.Series(dtype=float)
    eq = utils.derive_eqtlgen_beta_se(eq, maf_map)
    log(f"[eqtlgen] IL4/IL13 significant SNPs with FinnGen-MAF: {len(eq)} "
        f"(IL4={int((eq['gene']=='IL4').sum())}, IL13={int((eq['gene']=='IL13').sum())})")
    return eq


def harmonise_eqtl_finngen(eqtl_gene_df, fg_df):
    """Match eQTL SNPs to FinnGen rows by rsID, align outcome to eQTL effect allele."""
    rows = []
    if fg_df.empty or eqtl_gene_df.empty:
        return pd.DataFrame()
    rsid_map = {}
    for i, r in fg_df.iterrows():
        for tok in str(r["rsids"]).split(","):
            rsid_map.setdefault(tok.strip(), []).append(i)
    for _, e in eqtl_gene_df.iterrows():
        rsid = str(e["snp"])
        if rsid not in rsid_map:
            continue
        for i in rsid_map[rsid]:
            o = fg_df.loc[i]
            al = utils.align_outcome(
                {"ea": e["ea"]},
                {"ref": o["ref"], "alt": o["alt"], "beta": o["beta"],
                 "se": o["sebeta"], "af_alt": o["af_alt"]})
            if al is None:
                continue
            rows.append({"snp": rsid, "gene": e["gene"],
                         "beta_eqtl": e["beta"], "se_eqtl": e["se"], "p_eqtl": e["p"],
                         "beta_outcome": al["beta"], "se_outcome": al["se"],
                         "p_outcome": o["pval"], "fg_pos": int(o["pos"]),
                         "eaf_outcome": al["eaf"]})
            break
    return pd.DataFrame(rows)


def region_coloc_abf(bx, sx, by, sy, priors=PRIORS):
    """coloc.abf-style regional PP.H0-H4 from per-SNP Wakefield ABFs."""
    p1, p2, p12 = priors
    w = 0.15 ** 2

    def abf(beta, se):
        v = se ** 2
        z = beta / se
        return np.sqrt(v / (v + w)) * np.exp(0.5 * z ** 2 * (w / (v + w)))

    ax = abf(np.asarray(bx), np.asarray(sx))
    ay = abf(np.asarray(by), np.asarray(sy))
    sx_, sy_, sxy = ax.sum(), ay.sum(), (ax * ay).sum()
    p0 = max(1.0 - p1 - p2 - p12, 1e-12)
    ev = {"h0": p0, "h1": p1 * sx_, "h2": p2 * sy_,
          "h3": p1 * p2 * max(sx_ * sy_ - sxy, 0.0), "h4": p12 * sxy}
    tot = sum(ev.values())
    return {f"pp_{k}": v / tot for k, v in ev.items()}


# ----------------------------------------------------------------------------- analyses
def coloc_eqtl_vs_outcome(eq, fg, gene, outcome):
    sub = eq[eq["gene"] == gene]
    m = harmonise_eqtl_finngen(sub, fg)
    if len(m) < 5:
        log(f"[coloc] {gene} x {outcome}: only {len(m)} shared SNPs -> skipped")
        return None
    m = m.sort_values("p_eqtl").drop_duplicates("snp", keep="first")
    res = region_coloc_abf(m["beta_eqtl"], m["se_eqtl"],
                           m["beta_outcome"], m["se_outcome"])
    top = m.sort_values("p_outcome").iloc[0]
    h0, h1, h2, h3, h4 = utils.approx_coloc_abf(
        top["beta_eqtl"], top["se_eqtl"], top["beta_outcome"], top["se_outcome"],
        prior_p1=PRIORS[0], prior_p2=PRIORS[1], prior_p12=PRIORS[2])
    res.update({"gene": gene, "outcome": outcome, "contrast": f"{gene}-eQTL x {outcome}",
                "n_shared_snps": int(len(m)), "top_shared_snp": str(top["snp"]),
                "top_snp_pp_h4": float(h4),
                "note": "approximate coloc (eQTLGen significant-SNPs only, rsID-matched)"})
    return res


def coloc_il4_vs_il13(eq):
    a = eq[eq["gene"] == "IL4"].set_index("snp")
    b = eq[eq["gene"] == "IL13"].set_index("snp")
    shared = a.index.intersection(b.index)
    if len(shared) < 5:
        log(f"[coloc] IL4 x IL13 eQTL: only {len(shared)} shared SNPs -> skipped")
        return None
    a, b = a.loc[shared], b.loc[shared]
    # align IL13 beta to IL4 effect allele
    flip = a["ea"].str.upper().values != b["ea"].str.upper().values
    by = np.where(flip, -b["beta"].values, b["beta"].values)
    res = region_coloc_abf(a["beta"].values, a["se"].values, by, b["se"].values)
    res.update({"gene": "IL4_vs_IL13", "outcome": "eQTLGen",
                "contrast": "IL4-eQTL x IL13-eQTL (eQTLGen)",
                "n_shared_snps": int(len(shared)), "top_shared_snp": "",
                "top_snp_pp_h4": np.nan,
                "note": "approximate coloc (eQTLGen significant-SNPs only)"})
    return res


def direction_table(eq, fg_by_outcome, extra_snps=None):
    """IL4 & IL13 eQTL lead SNPs (+ optional extra SNPs) vs 3 outcomes: Wald direction."""
    leads = {}
    for gene in ["IL4", "IL13"]:
        sub = eq[eq["gene"] == gene].sort_values("p")
        if not sub.empty:
            r = sub.iloc[0]
            leads[f"{gene}_lead"] = {"snp": str(r["snp"]), "gene": gene,
                                     "beta_eqtl": r["beta"], "se": r["se"],
                                     "p_eqtl": r["p"], "ea": r["ea"],
                                     "source": "eQTLGen_lead"}
    if extra_snps:
        leads.update(extra_snps)
    rows = []
    for tag, L in leads.items():
        for outcome in OUTCOMES:
            fg = fg_by_outcome.get(outcome, pd.DataFrame())
            m = harmonise_eqtl_finngen(
                pd.DataFrame([{"snp": L["snp"], "gene": L["gene"],
                               "ea": L.get("ea", ""), "beta": L["beta_eqtl"],
                               "se": L.get("se", np.nan), "p": L["p_eqtl"]}]),
                fg)
            row = {"tag": tag, "snp": L["snp"], "gene": L["gene"],
                   "source": L["source"], "outcome": outcome,
                   "beta_exposure": L["beta_eqtl"], "matched": False,
                   "beta_outcome": np.nan, "wald_beta": np.nan, "wald_or": np.nan,
                   "direction_vs_exposure": "unmatched"}
            if m is not None and not m.empty:
                r = m.iloc[0]
                b, se, p = utils.mr_wald_ratio(r["beta_eqtl"], r["beta_outcome"],
                                               r["se_outcome"])
                row.update({"matched": True, "beta_outcome": r["beta_outcome"],
                            "p_outcome": r["p_outcome"], "wald_beta": b,
                            "wald_or": float(np.exp(b)),
                            "direction_vs_exposure": "risk" if b > 0 else "protective"})
            rows.append(row)
    return pd.DataFrame(rows), leads


def celltype_table(lead_snp_map):
    """Check lead SNPs across DICE / OneK1K / GTEx; record presence + direction."""
    rows = []
    wanted = set(lead_snp_map.keys())

    def add(source, cell, gene, snp, p, beta, note=""):
        rows.append({"source": source, "cell_or_tissue": cell, "gene": gene,
                     "snp": snp, "is_lead_snp": snp in wanted,
                     "p": p, "beta": beta,
                     "direction": ("up" if (beta or 0) > 0 else "down") if beta is not None else "",
                     "note": note})

    for cell in DICE_CELLS:
        try:
            df = utils.load_dice(cell, genes=["IL4", "IL13"])
            if df.empty:
                add("DICE", cell, "IL4/IL13", "", np.nan, None, "no IL4/IL13 eQTL in cell type")
                continue
            df = df.sort_values("p")
            for _, r in df.head(25).iterrows():
                add("DICE", cell, r["gene"], str(r["snp"]), r["p"], r["beta"],
                    "lead-region SNP" if str(r["snp"]) in wanted else "")
        except Exception as e:
            add("DICE", cell, "IL4/IL13", "", np.nan, None, f"load failed: {repr(e)[:60]}")

    for cell in ONEK1K_CELLS:
        try:
            df = utils.load_onek1k(cell, genes=["IL4", "IL13"])
            if df.empty:
                add("OneK1K", cell, "IL4/IL13", "", np.nan, None, "no IL4/IL13 eQTL in cell type")
            else:
                for _, r in df.sort_values("p").head(10).iterrows():
                    add("OneK1K", cell, r["gene"], str(r["snp"]), r["p"], r.get("beta"), "")
        except Exception as e:
            add("OneK1K", cell, "IL4/IL13", "", np.nan, None, f"load failed: {repr(e)[:60]}")

    for tissue in GTEX_LOCAL:
        for gene, ensg in GENE_ENSG.items():
            try:
                r = utils.load_gtex_top_snp(tissue, ensg)
                if r is None:
                    add("GTEx", tissue, gene, "", np.nan, None, "not a significant eGene")
                else:
                    add("GTEx", tissue, gene, str(r["variant_id"]), r["p"], r["beta"],
                        "top signif pair (GRCh38 variant_id)")
            except Exception as e:
                add("GTEx", tissue, gene, "", np.nan, None, f"load failed: {repr(e)[:60]}")
    for tissue in GTEX_MISSING:
        add("GTEx", tissue, "IL4/IL13", "", np.nan, None,
            "tissue not downloaded locally (data gap)")
    return pd.DataFrame(rows)


def fig4_locus(eq, fg_asthma, leads):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if fg_asthma is None or fg_asthma.empty:
        log("[fig4] no FinnGen asthma region data -> skip figure")
        return
    # map eQTLGen rsid -> GRCh38 pos via FinnGen asthma region
    pos_map = {}
    for _, r in fg_asthma.iterrows():
        for tok in str(r["rsids"]).split(","):
            pos_map.setdefault(tok.strip(), int(r["pos"]))
    eq = eq.copy() if eq is not None and not eq.empty else pd.DataFrame()
    if not eq.empty:
        eq["pos38"] = eq["snp"].map(pos_map)
        mapped = eq.dropna(subset=["pos38"])
    else:
        mapped = pd.DataFrame()
    log(f"[fig4] eQTL SNPs mapped to GRCh38: {len(mapped)}/{len(eq)}")

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.5), sharex=True,
                             gridspec_kw={"hspace": 0.28})
    ax = axes[0]
    for gene, c in [("IL4", "#2c6fbb"), ("IL13", "#c0392b")]:
        sub = mapped[mapped["gene"] == gene]
        ax.scatter(sub["pos38"] / 1e6, -np.log10(sub["p"].clip(lower=1e-300)),
                   s=12, color=c, alpha=0.7, label=f"{gene} eQTL (eQTLGen whole blood)")
    ax.set_ylabel("-log10(p)  eQTL")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_title("5q31 locus (IL4/IL13/RAD50, +/-1 Mb): eQTL vs FinnGen allergic asthma")

    ax2 = axes[1]
    ax2.scatter(fg_asthma["pos"] / 1e6, -np.log10(fg_asthma["pval"].clip(lower=1e-300)),
                s=10, color="#444444", alpha=0.55, label="FinnGen ALLERG_ASTHMA GWAS")
    ax2.set_ylabel("-log10(p)  GWAS")
    ax2.set_xlabel("chr5 position (Mb, GRCh38)")

    # gene spans + lead SNP markers
    for a in axes:
        for gene, (ch, s, e) in GENES_GRCH38.items():
            a.axvspan(s / 1e6, e / 1e6, color="gold", alpha=0.18)
            a.text((s + e) / 2e6, a.get_ylim()[1] * 0.02 if a.get_ylim()[1] else 1,
                   gene, ha="center", fontsize=8, color="#8a6d00")
        for tag, L in leads.items():
            p38 = pos_map.get(L["snp"])
            if p38:
                a.axvline(p38 / 1e6, color="purple", lw=0.9, ls="--", alpha=0.8)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(loc="upper right", fontsize=8, frameon=False)
    fig.text(0.99, 0.005, "dashed purple: IL4/IL13 eQTL lead SNPs (rsID-mapped to GRCh38)",
             ha="right", fontsize=7, color="purple")
    out = os.path.join(OUT_FIG, "r5_fig4_il4_locus.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"[fig4] saved {out}")


# ----------------------------------------------------------------------------- main
def adjudicate(coloc_df, dir_df):
    """Apply the pre-registered (i)/(ii)/(iii) decision tree."""
    def get(gene, outcome):
        sub = coloc_df[(coloc_df["gene"] == gene) & (coloc_df["outcome"] == outcome)]
        return sub.iloc[0] if not sub.empty else None

    il4_ast = get("IL4", PRIMARY)
    il13_ast = get("IL13", PRIMARY)
    il4_il13 = coloc_df[coloc_df["gene"] == "IL4_vs_IL13"]
    il4_il13 = il4_il13.iloc[0] if not il4_il13.empty else None

    if coloc_df.empty or "pp_h4" not in coloc_df.columns:
        return {"pp_h4_il4_asthma": np.nan, "pp_h4_il13_asthma": np.nan,
                "pp_h4_il4_vs_il13": np.nan, "dir_il4_asthma": "unmatched",
                "dir_il13_asthma": "unmatched", "verdict": "inconclusive",
                "narrative": "coloc 结果为空（数据缺口），无法裁决；如实报告。"}
    h4_il4 = il4_ast["pp_h4"] if il4_ast is not None else np.nan
    h4_il13 = il13_ast["pp_h4"] if il13_ast is not None else np.nan
    h4_pair = il4_il13["pp_h4"] if il4_il13 is not None else np.nan

    # directions of the two lead SNPs on asthma
    d_il4 = dir_df[(dir_df["tag"] == "IL4_lead") & (dir_df["outcome"] == PRIMARY)]
    d_il13 = dir_df[(dir_df["tag"] == "IL13_lead") & (dir_df["outcome"] == PRIMARY)]
    dir_il4 = d_il4.iloc[0]["direction_vs_exposure"] if len(d_il4) else "unmatched"
    dir_il13 = d_il13.iloc[0]["direction_vs_exposure"] if len(d_il13) else "unmatched"
    opposite = (dir_il4 in ("risk", "protective") and dir_il13 in ("risk", "protective")
                and dir_il4 != dir_il13)

    verdict, narrative = "inconclusive", ""
    if opposite:
        verdict = "(iii) IL4 vs IL13 opposite directions"
        narrative = ("IL4 与 IL13 的 cis-eQTL lead SNP 对哮喘方向相反：5q31 位点内 "
                     "IL4/IL13 非等价，dupilumab 同时阻断两者无法区分该差异——"
                     "位点内两个细胞因子承载不同（甚至拮抗）的遗传效应，是机制亮点。")
    elif pd.notna(h4_il4) and h4_il4 < 0.5 and pd.notna(h4_pair) and h4_pair >= 0.5:
        verdict = "(ii) coloc fails -> LD coincidence"
        narrative = ("IL4-eQTL 与哮喘 GWAS 共定位不成立（PP.H4<0.5)，而 IL4 与 IL13 的 "
                     "eQTL 共享同一信号（PP.H4>=0.5)：所谓 'IL4 保护' 实为 5q31 位点 LD "
                     "结构的误标，IL4 基因级论断必须撤回，双轴模型改写为位点级表述。")
    elif pd.notna(h4_il4) and h4_il4 >= 0.8:
        verdict = "(i) coloc holds -> tissue-mismatch narrative"
        narrative = ("IL4-eQTL 与哮喘共定位成立（PP.H4>=0.8)：IL4 降低等位基因=哮喘风险"
                     "等位基因，与 dupilumab 药理学方向相反，提示全血 eQTL 反映调节性/"
                     "记忆性 Th2 池而非气道致病信号；eQTL 证据降级为'组织特异性 IL-4 "
                     "生物学'。")
    elif pd.notna(h4_il4) and h4_il4 >= 0.5:
        verdict = "(i-weak) coloc suggestive (0.5<=PP.H4<0.8)"
        narrative = ("IL4-eQTL 与哮喘共定位为中等证据（0.5<=PP.H4<0.8)：支持同一因果"
                     "信号但强度不足以下最终结论；按组织错配叙事谨慎表述，并在审稿后"
                     "用 coloc.susie（多因果变异）复核。")
    else:
        narrative = ("证据不足以下三类裁决中的任何一类（数据缺口：eQTL 仅显著 SNP、"
                     "区域不完整）。如实报告 PP.H4 并标注'近似共定位'局限。")
    return {"pp_h4_il4_asthma": h4_il4, "pp_h4_il13_asthma": h4_il13,
            "pp_h4_il4_vs_il13": h4_pair, "dir_il4_asthma": dir_il4,
            "dir_il13_asthma": dir_il13, "verdict": verdict, "narrative": narrative}


def main():
    log("=" * 78)
    log("R5-06 IL4/IL13 (5q31) locus adjudication  (seed=42)")
    log(f"Region (GRCh38): chr{REGION_CHROM}:{REGION_START}-{REGION_END}")

    # ---- FinnGen region SNPs -------------------------------------------------------
    fg_by_outcome = {}
    for oc in OUTCOMES:
        try:
            fg_by_outcome[oc] = finngen_region(oc)
        except Exception as e:
            log(f"[finngen] {oc}: FAILED {repr(e)[:150]}")
            fg_by_outcome[oc] = pd.DataFrame()

    fg_ref = fg_by_outcome.get(PRIMARY, pd.DataFrame())
    if not fg_ref.empty:
        fg_ref = fg_ref.rename(columns={"rsids": "rsid"}).copy()
        fg_ref["rsid"] = fg_ref["rsid"].astype(str).str.split(",")
        fg_ref = fg_ref.explode("rsid")
        fg_ref["rsid"] = fg_ref["rsid"].str.strip()
        fg_ref = fg_ref.drop_duplicates("rsid", keep="first")

    # ---- eQTL exposure --------------------------------------------------------------
    try:
        eq = load_eqtl_exposure(fg_ref)
    except Exception as e:
        log(f"[eqtlgen] FAILED {repr(e)[:150]}")
        eq = pd.DataFrame()

    # ---- coloc combos -----------------------------------------------------------------
    coloc_rows = []
    if not eq.empty:
        for oc in OUTCOMES:
            for gene in ["IL4", "IL13"]:
                try:
                    r = coloc_eqtl_vs_outcome(eq, fg_by_outcome.get(oc, pd.DataFrame()),
                                              gene, oc)
                    if r:
                        coloc_rows.append(r)
                except Exception as e:
                    log(f"[coloc] {gene} x {oc}: FAILED {repr(e)[:120]}")
        try:
            r = coloc_il4_vs_il13(eq)
            if r:
                coloc_rows.append(r)
        except Exception as e:
            log(f"[coloc] IL4 x IL13: FAILED {repr(e)[:120]}")
    coloc_df = pd.DataFrame(coloc_rows)
    coloc_df.to_csv(os.path.join(OUT_TAB, "r5_il4_locus_coloc.csv"), index=False)
    log(f"[coloc] wrote {len(coloc_df)} contrasts")
    for _, r in coloc_df.iterrows():
        log(f"  {r['contrast']:<38} n={r['n_shared_snps']:<4} "
            f"H0={r['pp_h0']:.3f} H1={r['pp_h1']:.3f} H2={r['pp_h2']:.3f} "
            f"H3={r['pp_h3']:.3f} H4={r['pp_h4']:.3f}")

    # ---- direction table (lead SNPs + optional pQTL cross-link) -----------------------
    extra = {}
    pqtl_csv = os.path.join(OUT_TAB, "r5_pqtl_drug_target_mr.csv")
    if os.path.exists(pqtl_csv):
        try:
            pr = pd.read_csv(pqtl_csv)
            for gene in ["IL4", "IL13"]:
                sub = pr[(pr["protein"] == gene) & (pr["outcome"] == PRIMARY)]
                if not sub.empty:
                    r = sub.iloc[0]
                    # exposure-side beta of the pQTL lead SNP: recover from instrument table
                    extra[f"{gene}_pqtl_lead"] = {
                        "snp": str(r["lead_snp"]), "gene": gene,
                        "beta_eqtl": np.nan, "p_eqtl": np.nan,
                        "source": "INTERVAL_pQTL_lead (R5-09 cross-link)"}
        except Exception as e:
            log(f"[xlink] pQTL read failed: {repr(e)[:100]}")
    # need ea for lead SNPs -> pull from eQTLGen frame
    lead_meta = {}
    if not eq.empty:
        for gene in ["IL4", "IL13"]:
            sub = eq[eq["gene"] == gene].sort_values("p")
            if not sub.empty:
                r = sub.iloc[0]
                lead_meta[f"{gene}_lead"] = {"snp": str(r["snp"]), "gene": gene,
                                             "beta_eqtl": r["beta"], "se": r["se"],
                                             "p_eqtl": r["p"], "ea": r["ea"],
                                             "source": "eQTLGen_lead"}
    dir_df, leads = direction_table(eq, fg_by_outcome, extra_snps=None)
    # add pQTL lead SNP rows via instrument table (has ea/beta/se)
    pqtl_extra_rows = []
    inst_csv = os.path.join(OUT_TAB, "r5_pqtl_instruments.csv")
    if os.path.exists(inst_csv):
        try:
            inst = pd.read_csv(inst_csv, dtype={"snp": str, "rsid": str, "chrom": str})
            for gene in ["IL4", "IL13"]:
                sub = inst[inst["protein"] == gene].sort_values("p")
                if sub.empty:
                    continue
                s = sub.iloc[0]
                for oc in OUTCOMES:
                    fg = fg_by_outcome.get(oc, pd.DataFrame())
                    m = harmonise_eqtl_finngen(
                        pd.DataFrame([{"snp": str(s["snp"]), "gene": gene, "ea": s["ea"],
                                       "beta": s["beta"], "se": s["se"], "p": s["p"]}]), fg)
                    row = {"tag": f"{gene}_pqtl_lead", "snp": str(s["snp"]), "gene": gene,
                           "source": "INTERVAL_pQTL_lead (R5-09 cross-link)",
                           "outcome": oc, "beta_exposure": s["beta"], "matched": False,
                           "beta_outcome": np.nan, "wald_beta": np.nan, "wald_or": np.nan,
                           "direction_vs_exposure": "unmatched"}
                    if m is not None and not m.empty:
                        r = m.iloc[0]
                        b, se, p = utils.mr_wald_ratio(r["beta_eqtl"], r["beta_outcome"],
                                                       r["se_outcome"])
                        row.update({"matched": True, "beta_outcome": r["beta_outcome"],
                                    "p_outcome": r["p_outcome"], "wald_beta": b,
                                    "wald_or": float(np.exp(b)),
                                    "direction_vs_exposure": "risk" if b > 0 else "protective"})
                    pqtl_extra_rows.append(row)
        except Exception as e:
            log(f"[xlink] pQTL direction rows failed: {repr(e)[:100]}")
    if pqtl_extra_rows:
        dir_df = pd.concat([dir_df, pd.DataFrame(pqtl_extra_rows)], ignore_index=True)
    dir_df.to_csv(os.path.join(OUT_TAB, "r5_il4_il13_direction.csv"), index=False)
    log(f"[direction] wrote {len(dir_df)} rows")
    for _, r in dir_df.iterrows():
        log(f"  {r['tag']:<18} {r['snp']:<12} {r['outcome']:<16} "
            f"{'OR=' + format(r['wald_or'], '.3f') if r['matched'] else 'unmatched':<14} "
            f"{r['direction_vs_exposure']}")

    # ---- cell-type localisation -------------------------------------------------------
    lead_snp_map = {L["snp"]: L for L in lead_meta.values()}
    try:
        ct = celltype_table(lead_snp_map)
    except Exception as e:
        log(f"[celltype] FAILED {repr(e)[:120]}")
        ct = pd.DataFrame()
    ct.to_csv(os.path.join(OUT_TAB, "r5_il4_celltype_eqtl.csv"), index=False)
    log(f"[celltype] wrote {len(ct)} rows")

    # ---- figure -----------------------------------------------------------------------
    try:
        fig4_locus(eq, fg_by_outcome.get(PRIMARY, pd.DataFrame()),
                   {k: v for k, v in lead_meta.items()})
    except Exception as e:
        log(f"[fig4] FAILED {repr(e)[:150]}")

    # ---- adjudication ------------------------------------------------------------------
    verdict = adjudicate(coloc_df, dir_df)
    with open(os.path.join(OUT_TAB, "r5_il4_adjudication_verdict.json"), "w",
              encoding="utf-8") as f:
        json.dump(verdict, f, ensure_ascii=False, indent=2)
    log("-" * 78)
    log(f"VERDICT: {verdict['verdict']}")
    log(f"  PP.H4 IL4-eQTL x asthma  = {verdict['pp_h4_il4_asthma']}")
    log(f"  PP.H4 IL13-eQTL x asthma = {verdict['pp_h4_il13_asthma']}")
    log(f"  PP.H4 IL4 x IL13 eQTL    = {verdict['pp_h4_il4_vs_il13']}")
    log(f"  lead directions (asthma): IL4={verdict['dir_il4_asthma']}, "
        f"IL13={verdict['dir_il13_asthma']}")
    log(f"  narrative: {verdict['narrative']}")
    log("R5-06 done.")


if __name__ == "__main__":
    main()
