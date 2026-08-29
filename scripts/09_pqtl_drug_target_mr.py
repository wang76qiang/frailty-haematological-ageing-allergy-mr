#!/usr/bin/env python3
"""
R5-09: cis-pQTL drug-target Mendelian randomisation.

For each plasma pQTL summary-statistic file in data/real/interval_pqtl
(45+ proteins), take the lead cis SNP (gene TSS +/- 500 kb, GRCh38) after a
build-version self-check, harmonise against FinnGen R12 allergic outcomes
(asthma, allergic rhinitis, atopic dermatitis) plus two negative-control
outcomes (MI, RA), Steiger-filter, and estimate the causal effect with a
single-SNP Wald ratio.  BH-FDR across the protein x allergic-outcome grid;
regional approximate colocalisation (approx_coloc_abf kernel, coloc.abf-style
regional sum) for FDR < 0.1 signals; static drug-target map; faceted OR forest
plot.  Reuses src/r1/utils.py.  Seed = 42.  No pip installs.
"""

import os
import sys
import json
import time
import socket
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "r1"))
sys.path.insert(0, _HERE)
import utils  # noqa: E402
from _r509_worker import process_pqtl_file  # noqa: E402

BASE = utils.BASE_DIR
DATA = utils.DATA_DIR
PQTL_DIR = os.path.join(DATA, "interval_pqtl")
FG_DIR = os.path.join(DATA, "finngen_full")

OUT_T = os.path.join(BASE, "results", "r5", "tables")
OUT_F = os.path.join(BASE, "results", "r5", "figures")
OUT_L = os.path.join(BASE, "results", "r5", "logs")
PQTL_CACHE = os.path.join(OUT_T, "r509_pqtl_cache")
FG_CACHE = os.path.join(OUT_T, "r509_fg_cache")
for d in (OUT_T, OUT_F, OUT_L, PQTL_CACHE, FG_CACHE):
    os.makedirs(d, exist_ok=True)

SEED = 42
np.random.seed(SEED)

CIS_WINDOW = 500_000
BUILD_CHECK_WINDOW = 5_000_000
REGION_CACHE_WINDOW = 1_100_000
COLOC_WINDOW = 250_000
MIN_SNPS_BUILD_OK = 20
F_THRESH = 10.0

ALLERGIC = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]
NEG_CTRL = ["I9_MI_STRICT", "M13_RHEUMA"]
FG_FILES = {o: os.path.join(FG_DIR, f"finngen_R12_{o}.gz") for o in ALLERGIC + NEG_CTRL}
FG_N = {"ALLERG_ASTHMA": 13450 + 270290, "ALLERG_RHINITIS": 15569 + 474650,
        "L12_ATOPIC": 31245 + 432874, "I9_MI_STRICT": 31666 + 416171,
        "M13_RHEUMA": 16314 + 315115}

POSITIVE_CONTROL = ["IL4", "IL5", "IL13", "IL33", "TSLP"]
HYPOTHESIS = ["IL6", "TNF", "IFNG", "IL1A", "IL18", "IL8", "OSM", "S100A12", "CCL2"]

PROTEIN_TO_GENE = {
    "CCL2_MCP1": "CCL2", "CCL3_MIP1A": "CCL3", "IL8": "CXCL8",
    "LAP_TGFB1": "TGFB1", "OPG": "TNFRSF11B", "PDL1": "CD274",
    "SCF": "KITLG", "TNFB": "LTA", "TRAIL": "TNFSF10",
}

DRUG_MAP = {
    "TSLP": ("tezepelumab", "approved anti-TSLP mAb"),
    "IL5": ("mepolizumab; reslizumab", "approved anti-IL5 mAbs"),
    "IL13": ("lebrikizumab; tralokinumab", "approved anti-IL13 mAbs"),
    "IL4": ("dupilumab", "acts via IL4RA blockade; ligand (IL4/IL13) vs receptor distinction"),
    "IL6": ("tocilizumab", "acts via IL6R; direction depends on ligand-vs-receptor instrument"),
    "IL1A": ("anakinra", "IL1R1 antagonist; indirect for IL1A ligand"),
    "TNF": ("etanercept; infliximab; adalimumab", "approved anti-TNF biologics"),
    "IL33": ("itepekimab", "investigational anti-IL33 mAb"),
    "IFNG": ("emapalumab", "approved anti-IFN-gamma mAb"),
    "IL18": ("tadekinig alfa", "investigational IL-18BP"),
    "OSM": ("vixarelimab", "investigational anti-OSMR mAb"),
    "IL10": ("none marketed", ""),
    "IL17A": ("secukinumab; ixekizumab", "approved anti-IL17A mAbs"),
    "IL8": ("none marketed", "anti-CXCL8 programmes discontinued pre-approval"),
    "CCL2": ("carlumab", "anti-CCL2 mAb; development discontinued"),
    "S100A12": ("none marketed", ""),
    "IL2": ("aldesleukin", "recombinant IL-2; approved (oncology)"),
    "IL7": ("none", ""),
    "TNFB": ("none", "lymphotoxin-alpha; no approved drug"),
    "TRAIL": ("none", "TNFSF10 agonists trialled in oncology; none approved"),
    "VEGFA": ("bevacizumab", "approved anti-VEGFA mAb (oncology/ophthalmology)"),
}

PQTL_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
             "beta", "standard_error", "effect_allele_frequency", "p_value",
             "variant_id", "rsid", "n"]
FG_COLS = ["#chrom", "pos", "rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"]

LOG_LINES = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG_LINES.append(line)


# --------------------------------------------------------------------------- #
# Gene coordinates (target file + Ensembl REST cache, GRCh38)
# --------------------------------------------------------------------------- #
def ensembl_lookup(symbol, retries=3):
    url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
           f"?content-type=application/json")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            return {"gene": symbol, "ensembl_id": d.get("id", ""),
                    "chr": str(d.get("seq_region_name", "")),
                    "start": int(d["start"]), "end": int(d["end"]),
                    "strand": int(d.get("strand", 1))}
        except Exception as e:  # noqa: BLE001
            log(f"  Ensembl lookup failed for {symbol} (try {attempt + 1}): {repr(e)[:80]}")
            time.sleep(2 * (attempt + 1))
    return None


def get_gene_coords(genes):
    base = utils.load_gene_coords()
    coords = {}
    for _, r in base.iterrows():
        tss = int(r["start"]) if int(r["strand"]) == 1 else int(r["end"])
        coords[r["gene"]] = {"chr": str(r["chr"]), "start": int(r["start"]),
                             "end": int(r["end"]), "strand": int(r["strand"]), "tss": tss}
    cache_path = os.path.join(OUT_T, "r5_gene_coords_extended.csv")
    cache = pd.DataFrame()
    if os.path.exists(cache_path):
        cache = pd.read_csv(cache_path)
        for _, r in cache.iterrows():
            if r["gene"] not in coords:
                tss = int(r["start"]) if int(r["strand"]) == 1 else int(r["end"])
                coords[r["gene"]] = {"chr": str(r["chr"]), "start": int(r["start"]),
                                     "end": int(r["end"]), "strand": int(r["strand"]), "tss": tss}
    missing = [g for g in genes if g not in coords]
    new_rows = []
    for g in missing:
        socket.setdefaulttimeout(30)
        res = ensembl_lookup(g)
        if res is None:
            log(f"  !! no coordinates for {g}; protein will be skipped")
            continue
        tss = res["start"] if res["strand"] == 1 else res["end"]
        coords[g] = {**res, "tss": tss}
        new_rows.append(res)
    if new_rows:
        cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
        cache.to_csv(cache_path, index=False)
        log(f"  cached {len(new_rows)} new Ensembl coordinates -> {cache_path}")
    return coords, [g for g in missing if g not in coords]


# --------------------------------------------------------------------------- #
# FinnGen extraction
# --------------------------------------------------------------------------- #
def _read_fg_filtered(path, mask_fn):
    import gc
    rows = []
    reader = pd.read_csv(path, sep="\t", usecols=lambda c: c in FG_COLS,
                         chunksize=250_000, low_memory=False, dtype={"#chrom": str})
    for chunk in reader:
        chunk["#chrom"] = chunk["#chrom"].str.replace("chr", "", regex=False)
        keep = chunk[mask_fn(chunk)]
        if not keep.empty:
            rows.append(keep.copy())
        del chunk
        gc.collect()
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=FG_COLS)


def extract_fg_lead_hits(outcome, rsid_set, coord_set):
    cache = os.path.join(FG_CACHE, f"{outcome}_leadhits.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache, dtype={"#chrom": str})
    pattern = "|".join(sorted(rsid_set)) if rsid_set else None
    pos_by_chrom = {}
    for c, p in coord_set:
        pos_by_chrom.setdefault(str(c), set()).add(int(p))

    def mask(chunk):
        m = pd.Series(False, index=chunk.index)
        if pattern:
            m |= chunk["rsids"].astype(str).str.contains(pattern, regex=True, na=False)
        for c, posset in pos_by_chrom.items():
            m |= (chunk["#chrom"] == c) & chunk["pos"].isin(posset)
        return m

    df = _read_fg_filtered(FG_FILES[outcome], mask)
    df.to_csv(cache, index=False)
    return df


def extract_fg_region(outcome, chrom, start, end):
    cache = os.path.join(FG_CACHE, f"{outcome}_region_{chrom}_{start}_{end}.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache, dtype={"#chrom": str})

    def mask(chunk):
        return ((chunk["#chrom"] == str(chrom)) &
                (chunk["pos"] >= start) & (chunk["pos"] <= end))

    df = _read_fg_filtered(FG_FILES[outcome], mask)
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


def match_outcome(lead, by_rsid, by_coord):
    rsid = str(lead["rsid"])
    if rsid and rsid.lower() not in ("nan", "na") and rsid in by_rsid:
        return by_rsid[rsid], "rsid"
    o = by_coord.get((str(lead["chromosome"]), int(lead["base_pair_location"])))
    if o is not None:
        return o, "coord"
    return None, "unmatched"


def harmonise_lead(lead, o):
    ea = str(lead["effect_allele"]).upper()
    oa = str(lead["other_allele"]).upper()
    aligned = utils.align_outcome({"ea": ea}, o)
    if aligned is None:
        return None
    other_outcome = o["alt"].upper() if aligned["flipped"] else o["ref"].upper()
    if other_outcome != oa:
        return None
    return aligned


# --------------------------------------------------------------------------- #
# Regional approximate coloc
# --------------------------------------------------------------------------- #
def wakefield_abf(beta, se, w=0.15 ** 2):
    z = beta / se
    v = se ** 2
    return np.sqrt(v / (v + w)) * np.exp(0.5 * z ** 2 * (w / (v + w)))


def regional_coloc(exp_region, fg_region, p1=1e-4, p2=1e-4, p12=1e-5):
    by_rsid, by_coord = build_outcome_lookup(fg_region)
    pairs = []
    for _, e in exp_region.iterrows():
        o, how = match_outcome(e, by_rsid, by_coord)
        if o is None:
            continue
        al = harmonise_lead(e, o)
        if al is None:
            continue
        pairs.append({"snp": str(e["rsid"]), "pos": int(e["base_pair_location"]),
                      "beta_x": e["beta"], "se_x": e["standard_error"],
                      "p_x": e["p_value"], "beta_y": al["beta"], "se_y": al["se"],
                      "p_y": o["p"], "match": how})
    if not pairs:
        return pd.DataFrame(), None
    per = pd.DataFrame(pairs).drop_duplicates("snp")
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


def extract_pqtl_window(path, chrom, start, end, cache_name):
    cache = os.path.join(PQTL_CACHE, cache_name)
    if os.path.exists(cache):
        return pd.read_csv(cache, dtype={"chromosome": str})
    rows = []
    reader = pd.read_csv(path, sep="\t", usecols=PQTL_COLS, chunksize=500_000,
                         dtype={"chromosome": str})
    for chunk in reader:
        chunk["chromosome"] = chunk["chromosome"].str.replace("chr", "", regex=False)
        keep = chunk[(chunk["chromosome"] == str(chrom)) &
                     (chunk["base_pair_location"] >= start) &
                     (chunk["base_pair_location"] <= end)]
        if not keep.empty:
            rows.append(keep)
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=PQTL_COLS)
    df.to_csv(cache, index=False)
    return df


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    t0 = time.time()
    log("=== R5-09 cis-pQTL drug-target MR (single lead-SNP Wald ratio) ===")

    files = sorted(f for f in os.listdir(PQTL_DIR) if f.endswith(".tsv.gz"))
    proteins = []
    for f in files:
        prot = f.split("_GCST")[0]
        gene = PROTEIN_TO_GENE.get(prot, prot)
        proteins.append((prot, gene, os.path.join(PQTL_DIR, f)))
    prot2path = {p: path for p, _, path in proteins}
    log(f"found {len(proteins)} pQTL files")

    genes_needed = sorted({g for _, g, _ in proteins})
    coords, failed_genes = get_gene_coords(genes_needed)
    log(f"coordinates: {len(genes_needed) - len(failed_genes)}/{len(genes_needed)} ok; "
        f"failed: {failed_genes if failed_genes else 'none'}")

    tasks = [(prot, path, coords[gene]["chr"], coords[gene]["tss"], PQTL_CACHE,
              REGION_CACHE_WINDOW, BUILD_CHECK_WINDOW)
             for prot, gene, path in proteins if gene in coords]
    scan = {}
    use_pool = os.environ.get("R509_POOL", "0") == "1"
    if use_pool:
        try:
            with ProcessPoolExecutor(max_workers=3) as ex:
                futs = {ex.submit(process_pqtl_file, t): t[0] for t in tasks}
                for fut in as_completed(futs):
                    res = fut.result()
                    scan[res["protein"]] = res
            if not any(r.get("ok") for r in scan.values()):
                raise RuntimeError("all pool workers failed")
        except Exception as e:  # noqa: BLE001
            log(f"  pool failed ({repr(e)[:100]}); falling back to serial")
            scan = {}
            use_pool = False
    if not use_pool:
        for i, t in enumerate(tasks, 1):
            res = process_pqtl_file(t)
            scan[res["protein"]] = res
            tag = "cached" if res.get("cached") else "scanned"
            log(f"  [{i}/{len(tasks)}] {tag} {t[0]}: ok={res.get('ok')} "
                f"n5mb={res.get('n5mb')} err={res.get('error','')[:60]}")

    for prot in sorted(scan):
        r = scan[prot]
        if not r.get("ok"):
            log(f"  !! {prot}: scan error: {r.get('error')}")

    # ---- instrument table
    inst_rows = []
    for prot, gene, path in proteins:
        res = scan.get(prot, {"ok": False})
        if gene not in coords or not res.get("ok"):
            inst_rows.append({"protein": prot, "gene": gene, "status": "failed_scan",
                              "error": res.get("error", "no_coords")})
            continue
        c = coords[gene]
        region = pd.read_csv(os.path.join(PQTL_CACHE, f"{prot}_region.csv"),
                             dtype={"chromosome": str})
        cis = region[(region["base_pair_location"] >= c["tss"] - CIS_WINDOW) &
                     (region["base_pair_location"] <= c["tss"] + CIS_WINDOW)]
        if res["n5mb"] >= MIN_SNPS_BUILD_OK and not cis.empty:
            build_flag, lead = "grch38_cis", cis.loc[cis["p_value"].idxmin()]
        elif res["n5mb"] >= MIN_SNPS_BUILD_OK:
            build_flag = "grch38_wide_window"
            lead = region.loc[region["p_value"].idxmin()]
        else:
            build_flag = "grch37_fallback"
            lead = pd.read_csv(os.path.join(PQTL_CACHE, f"{prot}_globalmin.csv")).iloc[0]
            log(f"  !! {prot}: only {res['n5mb']} SNPs within 5Mb of {gene} TSS -> "
                f"GRCh37/mismatch fallback (genome-wide min-p SNP)")
        fstat = float((lead["beta"] / lead["standard_error"]) ** 2)
        inst_rows.append({
            "protein": prot, "gene": gene, "status": "ok",
            "chromosome": str(lead["chromosome"]),
            "base_pair_location": int(lead["base_pair_location"]),
            "variant_id": lead["variant_id"], "rsid": str(lead["rsid"]),
            "effect_allele": lead["effect_allele"], "other_allele": lead["other_allele"],
            "beta": lead["beta"], "standard_error": lead["standard_error"],
            "effect_allele_frequency": lead["effect_allele_frequency"],
            "p_value": lead["p_value"], "n": lead["n"],
            "fstat": round(fstat, 2), "build_flag": build_flag,
            "n_snps_cis": len(cis), "n_snps_5mb": res["n5mb"],
            "instrument_pass": bool(fstat > F_THRESH),
            "protein_category": (
                "positive_control" if (prot in POSITIVE_CONTROL or gene in POSITIVE_CONTROL)
                else "hypothesis" if (prot in HYPOTHESIS or gene in HYPOTHESIS)
                else "exploratory"),
        })
    inst = pd.DataFrame(inst_rows)
    inst.to_csv(os.path.join(OUT_T, "r5_pqtl_instrument_qc.csv"), index=False)
    ok = inst[inst["status"] == "ok"].copy()
    if len(ok) == 0:
        log("FATAL: no instruments built; writing empty outputs")
        pd.DataFrame().to_csv(os.path.join(OUT_T, "r5_pqtl_mr_results.csv"), index=False)
        pd.DataFrame().to_csv(os.path.join(OUT_T, "r5_pqtl_coloc.csv"), index=False)
        return
    ok["instrument_pass"] = ok["instrument_pass"].astype(bool)
    n_pass = int(ok["instrument_pass"].sum())
    log(f"instruments: {len(ok)} proteins; F>10: {n_pass}; "
        f"F<=10 (insufficient evidence, retained): {len(ok) - n_pass}")
    for r in ok[~ok["instrument_pass"]].itertuples():
        log(f"  insufficient evidence: {r.protein} (F={r.fstat:.2f}, lead={r.rsid})")

    # ---- FinnGen extraction for lead SNPs
    rsid_set = {r for r in ok["rsid"] if r and r.lower() not in ("nan", "na")}
    coord_set = set(zip(ok["chromosome"], ok["base_pair_location"]))
    log(f"extracting FinnGen rows for {len(rsid_set)} lead rsids across 5 outcomes ...")
    fg_look = {}
    for outcome in ALLERGIC + NEG_CTRL:
        df = extract_fg_lead_hits(outcome, rsid_set, coord_set)
        fg_look[outcome] = build_outcome_lookup(df)
        log(f"  {outcome}: {len(df)} candidate rows kept")

    # ---- harmonise + Steiger + Wald ratio
    mr_rows = []
    for r in ok.itertuples():
        lead = {"chromosome": r.chromosome, "base_pair_location": r.base_pair_location,
                "rsid": r.rsid, "effect_allele": r.effect_allele,
                "other_allele": r.other_allele, "beta": r.beta,
                "standard_error": r.standard_error}
        for outcome in ALLERGIC + NEG_CTRL:
            by_rsid, by_coord = fg_look[outcome]
            o, how = match_outcome(lead, by_rsid, by_coord)
            base_row = {"protein": r.protein, "gene": r.gene, "lead_snp": r.rsid,
                        "variant_id": r.variant_id, "chromosome": r.chromosome,
                        "pos": r.base_pair_location, "fstat": r.fstat,
                        "n_exp": int(r.n), "outcome": outcome, "n_outcome": FG_N[outcome],
                        "outcome_class": "allergic" if outcome in ALLERGIC else "negative_control",
                        "protein_category": r.protein_category,
                        "build_flag": r.build_flag, "match_method": how,
                        "instrument_pass": r.instrument_pass}
            if o is None:
                mr_rows.append({**base_row, "beta": np.nan, "se": np.nan, "p": np.nan,
                                "raw_p": np.nan, "or": np.nan, "or_lower": np.nan,
                                "or_upper": np.nan, "steiger_pass": False,
                                "p_outcome": np.nan})
                continue
            al = harmonise_lead(lead, o)
            if al is None:
                mr_rows.append({**base_row, "beta": np.nan, "se": np.nan, "p": np.nan,
                                "raw_p": np.nan, "or": np.nan, "or_lower": np.nan,
                                "or_upper": np.nan, "steiger_pass": False,
                                "p_outcome": o["p"], "match_method": "allele_mismatch"})
                continue
            sdf = pd.DataFrame([{"beta": r.beta, "se": r.standard_error,
                                 "beta_outcome": al["beta"], "se_outcome": al["se"]}])
            sdf = utils.steiger_filter(sdf, n_exp=int(r.n), n_out=FG_N[outcome])
            steiger_pass = bool(sdf["steiger_pass"].iloc[0])
            b, se, p = utils.mr_wald_ratio(r.beta, al["beta"], al["se"])
            mr_rows.append({**base_row, "beta": b, "se": se,
                            "p": p if steiger_pass else np.nan, "raw_p": p,
                            "or": np.exp(b), "or_lower": np.exp(b - 1.96 * se),
                            "or_upper": np.exp(b + 1.96 * se),
                            "steiger_pass": steiger_pass,
                            "steiger_p": float(sdf["steiger_p"].iloc[0]),
                            "beta_outcome_aligned": al["beta"], "p_outcome": o["p"]})
    mr = pd.DataFrame(mr_rows)
    mr["fdr"] = np.nan
    mask = mr["outcome"].isin(ALLERGIC) & mr["p"].notna()
    if mask.sum() > 0:
        mr.loc[mask, "fdr"] = multipletests(mr.loc[mask, "p"], method="fdr_bh")[1]
    mr.to_csv(os.path.join(OUT_T, "r5_pqtl_mr_results.csv"), index=False)
    sig = mr[mr["fdr"] < 0.1]
    log(f"FDR<0.1 signals (allergic, Steiger-passed): {len(sig)}")
    for _, s in sig.iterrows():
        log(f"  * {s['protein']} x {s['outcome']}: OR={s['or']:.3f} [{s['or_lower']:.3f},"
            f"{s['or_upper']:.3f}] FDR={s['fdr']:.3g} F={s['fstat']}")

    # ---- coloc for FDR<0.1 signals (+/-250 kb)
    coloc_rows = []
    for s in sig.itertuples():
        try:
            lo, hi = s.pos - COLOC_WINDOW, s.pos + COLOC_WINDOW
            if s.build_flag == "grch37_fallback":
                exp_reg = extract_pqtl_window(prot2path[s.protein], s.chromosome, lo, hi,
                                              f"{s.protein}_coloc_{s.chromosome}_{lo}_{hi}.csv")
            else:
                region = pd.read_csv(os.path.join(PQTL_CACHE, f"{s.protein}_region.csv"),
                                     dtype={"chromosome": str})
                exp_reg = region[(region["base_pair_location"] >= lo) &
                                 (region["base_pair_location"] <= hi)]
            fg_reg = extract_fg_region(s.outcome, s.chromosome, lo, hi)
            per, summ = regional_coloc(exp_reg, fg_reg)
            if summ is None:
                coloc_rows.append({"protein": s.protein, "outcome": s.outcome,
                                   "lead_snp": s.lead_snp, "n_snps": 0})
                continue
            per.to_csv(os.path.join(OUT_T, f"r5_pqtl_coloc_{s.protein}_{s.outcome}_persnp.csv"),
                       index=False)
            coloc_rows.append({"protein": s.protein, "outcome": s.outcome,
                               "lead_snp": s.lead_snp, **summ})
            log(f"  coloc {s.protein} x {s.outcome}: PP.H4={summ['pp_h4']:.3f} "
                f"(n={summ['n_snps']})")
        except Exception as e:  # noqa: BLE001
            log(f"  coloc failed {s.protein} x {s.outcome}: {repr(e)[:120]}")
            coloc_rows.append({"protein": s.protein, "outcome": s.outcome,
                               "lead_snp": s.lead_snp, "error": repr(e)[:120]})
    pd.DataFrame(coloc_rows).to_csv(os.path.join(OUT_T, "r5_pqtl_coloc.csv"), index=False)

    # ---- drug map
    best = (mr[mr["outcome"].isin(ALLERGIC) & mr["or"].notna()]
            .sort_values("fdr").groupby("protein", as_index=False).first())
    best_idx = best.set_index("protein") if len(best) else pd.DataFrame()
    best_gene_idx = best.set_index("gene") if len(best) else pd.DataFrame()
    dm_rows = []
    for target, (drugs, note) in DRUG_MAP.items():
        row = {"target": target, "drugs": drugs, "note": note}
        hit = None
        if len(best):
            if target in best_idx.index:
                hit = best_idx.loc[target]
            elif target in best_gene_idx.index:
                hit = best_gene_idx.loc[target]
        if hit is not None:
            b = hit
            row.update({"best_outcome": b["outcome"], "best_or": b["or"],
                        "best_fdr": b["fdr"]})
        dm_rows.append(row)
    pd.DataFrame(dm_rows).to_csv(os.path.join(OUT_T, "r5_pqtl_drug_map.csv"), index=False)

    # ---- forest plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plot_rows = []
    for prot in ok["protein"]:
        sub = mr[(mr["protein"] == prot) & mr["outcome"].isin(ALLERGIC) & mr["or"].notna()]
        if sub.empty:
            continue
        sub = sub.sort_values(["fdr", "raw_p"], na_position="last")
        b = sub.iloc[0]
        plot_rows.append({"protein": prot, "outcome": b["outcome"], "or": b["or"],
                          "lo": b["or_lower"], "hi": b["or_upper"], "fdr": b["fdr"],
                          "cat": b["protein_category"], "f_ok": bool(b["instrument_pass"])})
    pdf = pd.DataFrame(plot_rows)
    cats = [("positive_control", "Positive control (Th2 axis)"),
            ("hypothesis", "Hypothesis (inflammatory)"),
            ("exploratory", "Exploratory")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 8), sharex=True)
    for ax, (cat, title) in zip(axes, cats):
        sub = pdf[pdf["cat"] == cat].sort_values("or") if len(pdf) else pd.DataFrame()
        if sub.empty:
            ax.axis("off")
            continue
        y = np.arange(len(sub))
        colors = np.where(~sub["f_ok"], "#bbbbbb",
                          np.where(sub["fdr"] < 0.1, "#c0392b", "#2c6fbb"))
        ax.errorbar(sub["or"], y,
                    xerr=[(sub["or"] - sub["lo"]).clip(lower=1e-6),
                          (sub["hi"] - sub["or"]).clip(lower=1e-6)],
                    fmt="none", ecolor="#999999", elinewidth=1, capsize=2)
        ax.scatter(sub["or"], y, c=colors, s=45, zorder=3)
        ax.axvline(1, color="black", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{p}  ({o.replace('ALLERG_', '').replace('L12_', '')})"
                            for p, o in zip(sub["protein"], sub["outcome"])], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("OR per 1-SD higher protein level (log scale)")
        ax.grid(axis="x", alpha=0.3)
    legend = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b",
                     label="FDR<0.1 (F>10)"),
              Line2D([0], [0], marker="o", color="w", markerfacecolor="#2c6fbb",
                     label="NS (F>10)"),
              Line2D([0], [0], marker="o", color="w", markerfacecolor="#bbbbbb",
                     label="Weak instrument (F<=10)")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("R5-09 cis-pQTL drug-target MR vs FinnGen allergic outcomes", y=1.00)
    fig.tight_layout()
    utils.save_fig(fig, os.path.join(OUT_F, "r5_fig5_pqtl_forest.png"))
    log("forest plot saved")

    # ---- notes
    pos_hits = sig[sig["protein"].isin(POSITIVE_CONTROL)]
    neg = mr[mr["outcome_class"] == "negative_control"]
    neg_nom = neg[(neg["raw_p"].notna()) & (neg["raw_p"] < 0.05)]
    fb = ok[ok["build_flag"] != "grch38_cis"] if len(ok) else ok
    notes = ["## R5-09 cis-pQTL drug-target MR"]
    notes.append(f"- pQTL files processed: {len(proteins)}; instruments built: {len(ok)}; "
                 f"failed scans/coords: {len(inst) - len(ok)}")
    notes.append(f"- Instruments F>10 (main analysis): {n_pass}; F<=10 (insufficient-evidence "
                 f"table, retained): {len(ok) - n_pass}")
    for r in ok[~ok["instrument_pass"]].itertuples():
        notes.append(f"  - insufficient: {r.protein} F={r.fstat:.2f} lead={r.rsid}")
    notes.append(f"- Build self-check: {int((ok['build_flag']=='grch38_cis').sum())} files "
                 f"consistent with GRCh38 cis window; fallbacks: "
                 f"{', '.join(f'{r.protein}({r.build_flag})' for r in fb.itertuples()) or 'none'}")
    notes.append(f"- Positive-control hits (FDR<0.1, any allergic outcome): "
                 f"{len(pos_hits['protein'].unique())}/5 -> {sorted(pos_hits['protein'].unique())}")
    for _, s in sig.iterrows():
        notes.append(f"  - {s['protein']} x {s['outcome']}: OR={s['or']:.3f} "
                     f"({s['or_lower']:.3f}-{s['or_upper']:.3f}), FDR={s['fdr']:.3g}, F={s['fstat']}")
    notes.append(f"- Negative controls (MI/RA): {len(neg_nom)} nominal raw-p<0.05 of "
                 f"{int(neg['raw_p'].notna().sum())} tested")
    coloc_df = pd.DataFrame(coloc_rows)
    if not coloc_df.empty and "pp_h4" in coloc_df:
        for c in coloc_df.dropna(subset=["pp_h4"]).itertuples():
            verdict = ("COLOCALISED" if c.pp_h4 > 0.7 else
                       "suggestive" if c.pp_h4 > 0.5 else "not colocalised")
            notes.append(f"- Coloc {c.protein} x {c.outcome}: PP.H4={c.pp_h4:.3f} "
                         f"(n_snps={c.n_snps}) -> {verdict}")
    notes.append(f"- Runtime: {(time.time() - t0) / 60:.1f} min")
    with open(os.path.join(OUT_L, "r5_09_notes.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(notes) + "\n")
    with open(os.path.join(OUT_L, "r5_09_run.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG_LINES) + "\n")
    _update_combined_notes()
    log("=== R5-09 done ===")


def _update_combined_notes():
    parts = []
    for name in ("r5_09_notes.md", "r5_06_notes.md"):
        p = os.path.join(OUT_L, name)
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                parts.append(f.read().strip())
    if parts:
        with open(os.path.join(OUT_L, "R5_06_09_notes.md"), "w", encoding="utf-8") as f:
            f.write("# R5 revision notes: R5-09 + R5-06\n\n" + "\n\n".join(parts) + "\n")


if __name__ == "__main__":
    main()
