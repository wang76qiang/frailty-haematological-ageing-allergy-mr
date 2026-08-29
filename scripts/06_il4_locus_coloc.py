#!/usr/bin/env python3
"""
R5-06: IL4/IL13 locus colocalisation adjudication.

Region: IL4/IL13/RAD50 +/- 1 Mb on chr5 (GRCh38; union of IL4 and IL13 gene
bodies from target_gene_coords.csv).

Data reality (verified before writing this script):
  * eQTLGen whole blood: 204 FDR-significant cis SNPs for IL4, ZERO for IL13
    (IL13 is barely expressed in whole blood).
  * eQTL Catalogue target file: no IL4/IL13 rows.
  * DICE: IL4 eQTLs in TREG_NAIVE (46 SNPs); IL13 eQTLs in TH1 (6 SNPs).
  * OneK1K: no IL4/IL13 eQTL in any of the 6 cell types (full eqtl tables
    verified empty; re-checked here via top-eSNP tables).
  * GTEx v8 EUR on disk: only Adipose_Subcutaneous, Adipose_Visceral_Omentum,
    Adrenal_Gland, Artery_Aorta (Whole_Blood/Lung/Skin NOT downloaded).
    signif_pairs files are content-truncated (no IL13 pairs despite IL13 being
    a significant eGene), so gene-level *egenes* results are used instead:
    IL4 eGene in Artery_Aorta; IL13 eGene in Adipose_Subcutaneous /
    Adipose_Visceral_Omentum / Artery_Aorta.

Primary coloc: IL4 (eQTLGen) and IL13 (DICE TH1, best available) against
FinnGen R12 ALLERG_ASTHMA / ALLERG_RHINITIS / L12_ATOPIC, plus the
IL4-eQTL x IL13-eQTL self-colocalisation.  Priors p1=p2=1e-4, p12=1e-5;
per-SNP kernel = utils.approx_coloc_abf; regional combination = coloc.abf-style
sum of Wakefield ABFs.  Wald-direction concordance of the two lead SNPs across
the three outcomes; cell-type localisation table; dual-track locus figure.

Honest limitations (also written to the notes file):
  - Only FDR-significant cis SNPs are available for every eQTL source, so the
    regional coloc runs on the available-SNP intersection (trait-1 association
    pre-filtered); PP.H4 is biased upward relative to full-summary-stat coloc.
  - IL13 coloc uses DICE TH1 (n=91 donors, 6 SNPs): low power.

Seed = 42.  No pip installs.
"""

import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "r1"))
import utils  # noqa: E402

BASE = utils.BASE_DIR
DATA = utils.DATA_DIR
FG_DIR = os.path.join(DATA, "finngen_full")
OUT_T = os.path.join(BASE, "results", "r5", "tables")
OUT_F = os.path.join(BASE, "results", "r5", "figures")
OUT_L = os.path.join(BASE, "results", "r5", "logs")
FG_CACHE = os.path.join(OUT_T, "r506_fg_cache")
for d in (OUT_T, OUT_F, OUT_L, FG_CACHE):
    os.makedirs(d, exist_ok=True)

SEED = 42
np.random.seed(SEED)

CHROM = "5"
OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
FG_FILES = {o: os.path.join(FG_DIR, f"finngen_R12_{o}.gz") for o in OUTCOMES}
FG_N = {"ALLERG_ASTHMA": 13450 + 270290, "ALLERG_RHINITIS": 15569 + 474650,
        "L12_ATOPIC": 31245 + 432874}

GENES = ["IL4", "IL13"]
ENSG = {"IL4": "ENSG00000113520", "IL13": "ENSG00000169194"}
GRCH37_TO_38_OFFSET = 659963  # chr5 IL4 start 132014023 (b37) -> 132673986 (b38)

DICE_CELLS = ["TH1", "TH2", "TH17", "TREG", "CD4_NAIVE", "CD8_NAIVE",
              "B_CELL_NAIVE", "MONOCYTES", "M2", "NK", "TFH", "THSTAR",
              "TREG_MEM", "TREG_NAIVE", "CD4_STIM", "CD8_STIM"]
DICE_FOCUS = ["TH2", "TREG_NAIVE", "MONOCYTES"]
ONEK1K_CELLS = ["cd4nc", "cd4et", "cd8et", "monoc", "nk", "bin"]
GTEX_WANTED = ["Whole_Blood", "Lung", "Skin_Sun_Exposed_Lower_leg"]
GTEX_DIR = os.path.join(DATA, "gtex_eqtl", "GTEx_Analysis_v8_eQTL_EUR")

LOG_LINES = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


# --------------------------------------------------------------------------- #
# FinnGen region extraction
# --------------------------------------------------------------------------- #
FG_COLS = ["#chrom", "pos", "rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"]


def extract_fg_region(outcome, chrom, start, end):
    cache = os.path.join(FG_CACHE, f"{outcome}_region_{chrom}_{start}_{end}.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache, dtype={"#chrom": str})
    import gc
    rows = []
    reader = pd.read_csv(FG_FILES[outcome], sep="\t",
                         usecols=lambda c: c in FG_COLS, chunksize=250_000,
                         low_memory=False, dtype={"#chrom": str})
    for chunk in reader:
        chunk["#chrom"] = chunk["#chrom"].str.replace("chr", "", regex=False)
        keep = chunk[(chunk["#chrom"] == str(chrom)) &
                     (chunk["pos"] >= start) & (chunk["pos"] <= end)]
        if not keep.empty:
            rows.append(keep.copy())
        del chunk
        gc.collect()
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=FG_COLS)
    df.to_csv(cache, index=False)
    return df


def build_outcome_lookup(fg_df):
    by_rsid, by_coord = {}, {}
    for _, r in fg_df.iterrows():
        entry = {"chrom": str(r["#chrom"]), "pos": int(r["pos"]), "ref": r["ref"],
                 "alt": r["alt"], "beta": r["beta"], "se": r["sebeta"],
                 "af_alt": r["af_alt"], "p": r["pval"], "rsids": r["rsids"]}
        for rs in str(r["rsids"]).split(","):
            rs = rs.strip()
            if rs and rs.lower() != "nan":
                by_rsid.setdefault(rs, entry)
        by_coord.setdefault((str(r["#chrom"]), int(r["pos"])), entry)
    return by_rsid, by_coord


# --------------------------------------------------------------------------- #
# GTEx egenes (gene-level top eQTL; signif_pairs on disk is content-truncated)
# --------------------------------------------------------------------------- #
def gtex_egene_top(tissue, gene_id):
    path = os.path.join(GTEX_DIR, f"{tissue}.v8.EUR.egenes.txt.gz")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:  # noqa: BLE001
        log(f"  GTEx {tissue} egenes read error: {repr(e)[:60]}")
        return None
    sub = df[df["phenotype_id"].str.startswith(gene_id + ".")]
    if sub.empty:
        return None
    r = sub.iloc[0]
    parts = str(r["variant_id"]).split("_")
    return {"snp": str(r["variant_id"]), "chrom": parts[0].replace("chr", ""),
            "pos": int(parts[1]), "ref": parts[2], "alt": parts[3],
            "p": float(r["pval_nominal"]), "beta": float(r["slope"]),
            "se": float(r["slope_se"]), "qval": float(r["qval"])}


# --------------------------------------------------------------------------- #
# Regional coloc (coloc.abf-style sum of Wakefield ABFs over harmonised SNPs)
# --------------------------------------------------------------------------- #
def wakefield_abf(beta, se, w=0.15 ** 2):
    z = beta / se
    v = se ** 2
    return np.sqrt(v / (v + w)) * np.exp(0.5 * z ** 2 * (w / (v + w)))


def regional_coloc_from_pairs(pairs, p1=1e-4, p2=1e-4, p12=1e-5):
    if pairs.empty:
        return pd.DataFrame(), None
    per = pairs.drop_duplicates("snp").copy()
    h = [utils.approx_coloc_abf(r.beta_x, r.se_x, r.beta_y, r.se_y,
                                prior_p1=p1, prior_p2=p2, prior_p12=p12)
         for r in per.itertuples()]
    per[["pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4"]] = pd.DataFrame(h, index=per.index)
    abf1 = wakefield_abf(per["beta_x"].values, per["se_x"].values)
    abf2 = wakefield_abf(per["beta_y"].values, per["se_y"].values)
    s1, s2 = abf1.sum(), abf2.sum()
    s12 = float((abf1 * abf2).sum())
    w0, w1, w2 = 1.0, p1 * s1, p2 * s2
    w3 = p1 * p2 * max(s1 * s2 - s12, 0.0)
    w4 = p12 * s12
    tot = w0 + w1 + w2 + w3 + w4
    summary = {"n_snps": len(per), "pp_h0": w0 / tot, "pp_h1": w1 / tot,
               "pp_h2": w2 / tot, "pp_h3": w3 / tot, "pp_h4": w4 / tot,
               "max_snp_pph4": float(per["pp_h4"].max()),
               "lead_coloc_snp": str(per.loc[per["pp_h4"].idxmax(), "snp"])}
    return per, summary


def harmonise_eqtl_vs_fg(eqtl_df, fg_lookup):
    by_rsid, _ = fg_lookup
    rows = []
    for _, e in eqtl_df.iterrows():
        o = by_rsid.get(str(e["snp"]))
        if o is None:
            continue
        ea = str(e["ea"]).upper()
        al = utils.align_outcome({"ea": ea}, o)
        if al is None:
            continue
        other = o["alt"].upper() if al["flipped"] else o["ref"].upper()
        if other != str(e["oa"]).upper():
            continue
        rows.append({"snp": str(e["snp"]), "pos37": int(e["pos"]),
                     "pos38": o["pos"], "beta_x": e["beta"], "se_x": e["se"],
                     "p_x": e["p"], "beta_y": al["beta"], "se_y": al["se"],
                     "p_y": o["p"]})
    return pd.DataFrame(rows)


def harmonise_eqtl_vs_eqtl(a_df, b_df):
    b = b_df.set_index("snp")
    rows = []
    for _, e in a_df.iterrows():
        if e["snp"] not in b.index:
            continue
        o = b.loc[e["snp"]]
        ea, oa = str(e["ea"]).upper(), str(e["oa"]).upper()
        if {ea, oa} != {str(o["ea"]).upper(), str(o["oa"]).upper()}:
            continue
        beta2 = o["beta"] if ea == str(o["ea"]).upper() else -o["beta"]
        rows.append({"snp": str(e["snp"]), "pos37": int(e["pos"]),
                     "beta_x": e["beta"], "se_x": e["se"], "p_x": e["p"],
                     "beta_y": beta2, "se_y": o["se"], "p_y": o["p"]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    t0 = time.time()
    log("=== R5-06 IL4/IL13 locus colocalisation ===")

    coords = utils.load_gene_coords().set_index("gene")
    lo = int(min(coords.loc["IL4", "start"], coords.loc["IL13", "start"])) - 1_000_000
    hi = int(max(coords.loc["IL4", "end"], coords.loc["IL13", "end"])) + 1_000_000
    log(f"region: chr{CHROM}:{lo}-{hi} (GRCh38, IL4/IL13 union +/-1Mb)")

    fg_lookup = {}
    fg_region = {}
    for o in OUTCOMES:
        df = extract_fg_region(o, CHROM, lo, hi)
        fg_region[o] = df
        fg_lookup[o] = build_outcome_lookup(df)
        log(f"  {o}: {len(df)} region SNPs")

    # ---- eQTL frames
    eqg_raw = utils.load_eqtlgen(genes=GENES)
    maf_map = {}
    for o in OUTCOMES:
        maf_map.update({rs: e["af_alt"] for rs, e in fg_lookup[o][0].items()})
    eqg = utils.derive_eqtlgen_beta_se(eqg_raw.copy(), pd.Series(maf_map))
    eqg = eqg[eqg["gene"].isin(GENES)].reset_index(drop=True)
    il4_eqtlgen = eqg[eqg["gene"] == "IL4"].sort_values("p").reset_index(drop=True)
    log(f"eQTLGen: IL4={len(il4_eqtlgen)} SNPs, "
        f"IL13={int((eqg['gene']=='IL13').sum())} SNPs (significant-only file)")

    # DICE frames (GRCh37 positions, rsid-matched to FinnGen downstream)
    dice_th1_il13 = utils.load_dice("TH1", genes=["IL13"])
    dice_th1_il13 = dice_th1_il13.sort_values("p").reset_index(drop=True)
    dice_treg_il4 = utils.load_dice("TREG_NAIVE", genes=["IL4"])
    dice_treg_il4 = dice_treg_il4.sort_values("p").reset_index(drop=True)
    log(f"DICE: TH1 IL13={len(dice_th1_il13)} SNPs; TREG_NAIVE IL4={len(dice_treg_il4)} SNPs")

    eqtl_sets = []  # (label, gene, df)
    if len(il4_eqtlgen):
        eqtl_sets.append(("IL4_eQTLGen", "IL4", il4_eqtlgen))
    if len(dice_th1_il13):
        eqtl_sets.append(("IL13_DICE_TH1", "IL13", dice_th1_il13))
    if len(dice_treg_il4):
        eqtl_sets.append(("IL4_DICE_TREG_NAIVE", "IL4", dice_treg_il4))

    lead_il4 = il4_eqtlgen.iloc[0] if len(il4_eqtlgen) else dice_treg_il4.iloc[0]
    lead_il13 = dice_th1_il13.iloc[0]
    log(f"  IL4 lead: {lead_il4['snp']} (p={lead_il4['p']:.2e}); "
        f"IL13 lead: {lead_il13['snp']} (p={lead_il13['p']:.2e}, DICE TH1)")

    # ---- coloc combos x outcomes
    coloc_rows = []
    for label, gene, gdf in eqtl_sets:
        for o in OUTCOMES:
            try:
                pairs = harmonise_eqtl_vs_fg(gdf, fg_lookup[o])
                per, summ = regional_coloc_from_pairs(pairs)
                if summ is None:
                    coloc_rows.append({"trait1": label, "trait2": o, "n_snps": 0})
                    continue
                per.insert(0, "trait1", label)
                per.insert(1, "trait2", o)
                per.to_csv(os.path.join(OUT_T, f"r5_il4_locus_coloc_{label}_{o}_persnp.csv"),
                           index=False)
                coloc_rows.append({"trait1": label, "trait2": o, **summ})
                log(f"  coloc {label} x {o}: PP.H4={summ['pp_h4']:.3f} (n={summ['n_snps']})")
            except Exception as e:  # noqa: BLE001
                log(f"  coloc {label} x {o} failed: {repr(e)[:120]}")
    # IL4 x IL13 self-coloc (eQTLGen IL4 vs DICE TH1 IL13, rsid-matched)
    try:
        pairs = harmonise_eqtl_vs_eqtl(il4_eqtlgen, dice_th1_il13)
        per, summ = regional_coloc_from_pairs(pairs)
        if summ is not None:
            per.insert(0, "trait1", "IL4_eQTLGen")
            per.insert(1, "trait2", "IL13_DICE_TH1")
            per.to_csv(os.path.join(OUT_T, "r5_il4_locus_coloc_IL4xIL13_persnp.csv"),
                       index=False)
            coloc_rows.append({"trait1": "IL4_eQTLGen", "trait2": "IL13_DICE_TH1", **summ})
            log(f"  coloc IL4-eQTL x IL13-eQTL: PP.H4={summ['pp_h4']:.3f} (n={summ['n_snps']})")
        else:
            coloc_rows.append({"trait1": "IL4_eQTLGen", "trait2": "IL13_DICE_TH1", "n_snps": 0})
            log("  coloc IL4-eQTL x IL13-eQTL: no shared rsids")
    except Exception as e:  # noqa: BLE001
        log(f"  IL4xIL13 coloc failed: {repr(e)[:120]}")
    coloc_df = pd.DataFrame(coloc_rows)
    coloc_df.to_csv(os.path.join(OUT_T, "r5_il4_locus_coloc.csv"), index=False)

    # ---- Wald direction table: lead SNPs x 3 outcomes
    dir_rows = []
    for gene, lead in (("IL4", lead_il4), ("IL13", lead_il13)):
        for o in OUTCOMES:
            by_rsid, _ = fg_lookup[o]
            fg_hit = by_rsid.get(str(lead["snp"]))
            row = {"gene": gene, "lead_snp": str(lead["snp"]),
                   "eqtl_source": "eQTLGen" if gene == "IL4" else "DICE_TH1",
                   "eqtl_beta": lead["beta"], "eqtl_se": lead["se"],
                   "eqtl_p": lead["p"], "outcome": o}
            if fg_hit is None:
                row.update({"gwas_beta": np.nan, "gwas_se": np.nan, "gwas_p": np.nan,
                            "wald_beta": np.nan, "wald_p": np.nan, "direction": "unmatched"})
            else:
                al = utils.align_outcome({"ea": str(lead["ea"]).upper()}, fg_hit)
                if al is None:
                    row.update({"gwas_beta": np.nan, "gwas_se": np.nan,
                                "gwas_p": fg_hit["p"], "wald_beta": np.nan,
                                "wald_p": np.nan, "direction": "allele_fail"})
                else:
                    b, se, p = utils.mr_wald_ratio(lead["beta"], al["beta"], al["se"])
                    row.update({"gwas_beta": al["beta"], "gwas_se": al["se"],
                                "gwas_p": fg_hit["p"], "wald_beta": b, "wald_se": se,
                                "wald_p": p,
                                "direction": "risk" if b > 0 else "protective"})
            dir_rows.append(row)
    dir_df = pd.DataFrame(dir_rows)
    dir_df.to_csv(os.path.join(OUT_T, "r5_il4_il13_direction.csv"), index=False)
    for d in dir_df.itertuples():
        log(f"  direction {d.gene} {d.lead_snp} x {d.outcome}: {d.direction}")

    # ---- cell-type localisation of lead SNPs
    lead_rsids = {str(lead_il4["snp"]): "IL4", str(lead_il13["snp"]): "IL13"}
    ct_rows = []
    for cell in DICE_CELLS:
        path = os.path.join(DATA, "dice_eqtl", f"{cell}.vcf")
        if not os.path.exists(path):
            ct_rows.append({"gene": "IL4/IL13", "source": "DICE", "cell_or_tissue": cell,
                            "is_focus": cell in DICE_FOCUS, "any_eqtl_for_gene": False,
                            "lead_snp_present": False, "note": "vcf not on disk"})
            continue
        try:
            d = utils.load_dice(cell, genes=GENES)
        except Exception as e:  # noqa: BLE001
            log(f"  DICE {cell} load failed: {repr(e)[:80]}")
            d = pd.DataFrame()
        if d.empty:
            ct_rows.append({"gene": "IL4/IL13", "source": "DICE", "cell_or_tissue": cell,
                            "is_focus": cell in DICE_FOCUS, "any_eqtl_for_gene": False,
                            "lead_snp_present": False})
            continue
        for g in GENES:
            gsub = d[d["gene"] == g]
            if gsub.empty:
                continue
            top = gsub.sort_values("p").iloc[0]
            ct_rows.append({"gene": g, "source": "DICE", "cell_or_tissue": cell,
                            "is_focus": cell in DICE_FOCUS, "any_eqtl_for_gene": True,
                            "lead_snp_present": bool(top["snp"] in lead_rsids),
                            "top_snp": str(top["snp"]), "top_p": top["p"],
                            "top_beta": top["beta"],
                            "direction": "up" if top["beta"] > 0 else "down"})
    # OneK1K via top-eSNP tables (full eqtl tables verified empty separately)
    for cell in ONEK1K_CELLS:
        esnp = os.path.join(DATA, "onek1k", f"{cell}_esnp_table.tsv.gz")
        try:
            k = utils.load_onek1k(cell, genes=GENES, path=esnp)
        except Exception as e:  # noqa: BLE001
            log(f"  OneK1K {cell} load failed: {repr(e)[:80]}")
            k = pd.DataFrame()
        if k.empty:
            ct_rows.append({"gene": "IL4/IL13", "source": "OneK1K",
                            "cell_or_tissue": cell, "is_focus": False,
                            "any_eqtl_for_gene": False, "lead_snp_present": False,
                            "note": "absent in top-eSNP table; full eqtl table verified empty"})
        else:
            for g in GENES:
                gsub = k[k["gene"] == g]
                if not gsub.empty:
                    top = gsub.sort_values("p").iloc[0]
                    ct_rows.append({"gene": g, "source": "OneK1K", "cell_or_tissue": cell,
                                    "is_focus": False, "any_eqtl_for_gene": True,
                                    "lead_snp_present": bool(top["snp"] in lead_rsids),
                                    "top_snp": str(top["snp"]), "top_p": top["p"],
                                    "top_beta": top["beta"],
                                    "direction": "up" if top["beta"] > 0 else "down"})
    # GTEx egenes
    gtex_avail = sorted({f.split(".v8.EUR")[0] for f in os.listdir(GTEX_DIR)
                         if f.endswith("egenes.txt.gz")})
    wanted_present = [t for t in GTEX_WANTED if t in gtex_avail]
    log(f"  GTEx wanted {GTEX_WANTED}; present: {wanted_present}; using: {gtex_avail}")
    lead_pos38 = {rs: fg_lookup["ALLERG_ASTHMA"][0].get(rs, {}).get("pos")
                  for rs in lead_rsids}
    for tissue in gtex_avail:
        for g in GENES:
            top = gtex_egene_top(tissue, ENSG[g])
            if top is None:
                ct_rows.append({"gene": g, "source": "GTEx", "cell_or_tissue": tissue,
                                "is_focus": tissue in GTEX_WANTED,
                                "any_eqtl_for_gene": False, "lead_snp_present": False})
            else:
                ct_rows.append({"gene": g, "source": "GTEx", "cell_or_tissue": tissue,
                                "is_focus": tissue in GTEX_WANTED,
                                "any_eqtl_for_gene": True,
                                "lead_snp_present": bool(top["pos"] in
                                                         [p for p in lead_pos38.values() if p]),
                                "top_snp": top["snp"], "top_p": top["p"],
                                "top_beta": top["beta"],
                                "direction": "up" if top["beta"] > 0 else "down"})
    ct_df = pd.DataFrame(ct_rows)
    ct_df.to_csv(os.path.join(OUT_T, "r5_il4_celltype_eqtl.csv"), index=False)

    # ---- figure: dual-track locus plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rsid2pos38 = fg_lookup["ALLERG_ASTHMA"][0]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True,
                             gridspec_kw={"hspace": 0.28})
    ax = axes[0]
    if len(il4_eqtlgen):
        xs, ys = [], []
        for r in il4_eqtlgen.itertuples():
            hit = rsid2pos38.get(str(r.snp))
            xs.append(hit["pos"] if hit else r.pos + GRCH37_TO_38_OFFSET)
            ys.append(-np.log10(max(r.p, 1e-300)))
        ax.scatter(xs, ys, c="#2c6fbb", s=14, alpha=0.75, label="IL4 eQTL (eQTLGen blood)")
    if len(dice_th1_il13):
        xs, ys = [], []
        for r in dice_th1_il13.itertuples():
            hit = rsid2pos38.get(str(r.snp))
            xs.append(hit["pos"] if hit else r.pos + GRCH37_TO_38_OFFSET)
            ys.append(-np.log10(max(r.p, 1e-300)))
        ax.scatter(xs, ys, c="#27ae60", s=26, marker="D", alpha=0.9,
                   label="IL13 eQTL (DICE TH1)")
    for tissue in gtex_avail:
        for g, col in (("IL4", "#8e44ad"), ("IL13", "#16a085")):
            top = gtex_egene_top(tissue, ENSG[g])
            if top and lo <= top["pos"] <= hi:
                ax.scatter([top["pos"]], [-np.log10(top["p"])], c=col, s=60,
                           marker="*", edgecolor="k", zorder=5)
    ax.scatter([], [], c="#8e44ad", marker="*", edgecolor="k", s=60,
               label="GTEx eGene top (IL4/IL13)")
    ax.set_ylabel("-log10(p) eQTL")
    ax.legend(loc="upper left", fontsize=7.5, frameon=False)
    ax.set_title(f"IL4/IL13 locus  chr{CHROM}:{lo:,}-{hi:,} (GRCh38)", fontsize=10)
    ax2 = axes[1]
    gwas = fg_region["ALLERG_ASTHMA"]
    ax2.scatter(gwas["pos"], -np.log10(gwas["pval"].clip(lower=1e-300)),
                c="#c0392b", s=10, alpha=0.6)
    ax2.set_ylabel("-log10(p) FinnGen\nALLERG_ASTHMA")
    ax2.set_xlabel(f"chr{CHROM} position (bp)")
    for rs, gname in lead_rsids.items():
        hit = rsid2pos38.get(rs)
        if hit:
            for a in axes:
                a.axvline(hit["pos"], color="black", ls="--", lw=0.7, alpha=0.6)
            axes[0].annotate(f"{rs}\n({gname} lead)", (hit["pos"], axes[0].get_ylim()[1]),
                             fontsize=7, ha="right", va="top")
    for a in axes:
        for gname in GENES:
            a.axvspan(int(coords.loc[gname, "start"]), int(coords.loc[gname, "end"]),
                      color="gray", alpha=0.12)
    utils.save_fig(fig, os.path.join(OUT_F, "r5_fig4_il4_locus.png"))
    log("locus figure saved")

    # ---- adjudication + notes
    def verdict(h4):
        if h4 is None or pd.isna(h4):
            return "NA"
        return "COLOCALISED" if h4 > 0.7 else ("suggestive" if h4 > 0.5 else "NOT colocalised")

    notes = ["## R5-06 IL4/IL13 locus colocalisation adjudication"]
    notes.append(f"- Region: chr{CHROM}:{lo}-{hi} (GRCh38; IL4+IL13 bodies +/-1Mb).")
    notes.append(f"- eQTL availability: eQTLGen IL4={len(il4_eqtlgen)} sig cis SNPs, "
                 f"IL13=0 (IL13 not detectably expressed in whole blood); eQTL Catalogue "
                 f"has no IL4/IL13; DICE IL4 in TREG_NAIVE ({len(dice_treg_il4)} SNPs), "
                 f"IL13 in TH1 ({len(dice_th1_il13)} SNPs); OneK1K: none in 6 cell types; "
                 f"GTEx on disk: {gtex_avail} (Whole_Blood/Lung/Skin NOT downloaded).")
    notes.append("- LIMITATION: all eQTL sources provide FDR-significant SNPs only, so the "
                 "regional coloc runs on the available-SNP intersection (trait-1 "
                 "pre-filtered on significance), not a full summary-statistic region; "
                 "PP.H4 is biased upward. IL13 coloc uses DICE TH1 (n=91, 6 SNPs): low power.")
    for c in coloc_df.itertuples():
        h4 = getattr(c, "pp_h4", np.nan)
        if pd.notna(h4):
            notes.append(f"- Coloc {c.trait1} x {c.trait2}: PP.H0={c.pp_h0:.3f} "
                         f"PP.H1={c.pp_h1:.3f} PP.H2={c.pp_h2:.3f} PP.H3={c.pp_h3:.3f} "
                         f"PP.H4={c.pp_h4:.3f} (n={c.n_snps}) -> {verdict(h4)}")
        else:
            notes.append(f"- Coloc {c.trait1} x {c.trait2}: no overlapping SNPs")
    for o in OUTCOMES:
        sub = dir_df[dir_df["outcome"] == o].set_index("gene")
        d4, d13 = sub.loc["IL4", "direction"], sub.loc["IL13", "direction"]
        same = (d4 == d13) and d4 in ("risk", "protective")
        notes.append(f"- Direction {o}: IL4({sub.loc['IL4','lead_snp']})={d4}, "
                     f"IL13({sub.loc['IL13','lead_snp']})={d13} -> "
                     f"{'CONCORDANT' if same else 'discordant/unmatched'}")
    dice_hits = ct_df[(ct_df["source"] == "DICE") & (ct_df["any_eqtl_for_gene"] == True)]
    notes.append("- DICE cells with IL4/IL13 eQTL: " +
                 (", ".join(f"{r.cell_or_tissue}({r.gene},top={r.top_snp})"
                            for r in dice_hits.itertuples()) or "none") +
                 "; TH2 and MONOCYTES: none.")
    onek_hits = ct_df[(ct_df["source"] == "OneK1K") & (ct_df["any_eqtl_for_gene"] == True)]
    notes.append("- OneK1K: " + ("present in " + ", ".join(onek_hits["cell_or_tissue"])
                                  if len(onek_hits)
                                  else "no IL4/IL13 eQTL in any of 6 cell types "
                                       "(top-eSNP tables; full eqtl tables verified empty)"))
    gtex_hits = ct_df[(ct_df["source"] == "GTEx") & (ct_df["any_eqtl_for_gene"] == True)]
    notes.append("- GTEx (egenes; signif_pairs on disk content-truncated): " +
                 (", ".join(f"{r.cell_or_tissue}({r.gene},{r.top_snp},"
                             f"beta={r.top_beta:.2f})" for r in gtex_hits.itertuples())
                  or "none"))
    def _h4(t1, t2):
        sub = coloc_df[(coloc_df["trait1"] == t1) & (coloc_df["trait2"] == t2)]["pp_h4"]
        return float(sub.iloc[0]) if len(sub) and pd.notna(sub.iloc[0]) else float("nan")
    notes.append(f"- ADJUDICATION: IL4(eQTLGen)xASTHMA PP.H4={_h4('IL4_eQTLGen','ALLERG_ASTHMA'):.3f}; "
                 f"IL13(DICE_TH1)xASTHMA PP.H4={_h4('IL13_DICE_TH1','ALLERG_ASTHMA'):.3f}; "
                 f"IL4xIL13 eQTL PP.H4={_h4('IL4_eQTLGen','IL13_DICE_TH1'):.3f}")
    notes.append(f"- Runtime: {(time.time() - t0) / 60:.1f} min")
    with open(os.path.join(OUT_L, "r5_06_notes.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(notes) + "\n")
    with open(os.path.join(OUT_L, "r5_06_run.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG_LINES) + "\n")
    parts = []
    for name in ("r5_09_notes.md", "r5_06_notes.md"):
        p = os.path.join(OUT_L, name)
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                parts.append(f.read().strip())
    if parts:
        with open(os.path.join(OUT_L, "R5_06_09_notes.md"), "w", encoding="utf-8") as f:
            f.write("# R5 revision notes: R5-09 + R5-06\n\n" + "\n\n".join(parts) + "\n")
    log("=== R5-06 done ===")


if __name__ == "__main__":
    main()
