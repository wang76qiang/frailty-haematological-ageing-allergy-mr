#!/usr/bin/env python3
"""
R5-10: EWAS Catalog CpG convergence between aging and allergy EWAS.

Inputs
------
data/real/ewas_catalog/ewascatalog-studies.txt.gz   (study-level metadata)
data/real/ewas_catalog/ewascatalog-results.txt.gz   (CpG-level associations)

Outputs
-------
results/r5/tables/r5_ewas_overlap_enrichment.csv
results/r5/tables/r5_ewas_overlap_genes.csv
results/r5/figures/r5_ewas_convergence.png
results/r5/logs/r5_10_summary.json
"""
import os
import re
import sys
import json
import gzip

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "r1"))
from utils import save_fig  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE, "data", "real", "ewas_catalog")
OUT_TABLES = os.path.join(BASE, "results", "r5", "tables")
OUT_FIGS = os.path.join(BASE, "results", "r5", "figures")
OUT_LOGS = os.path.join(BASE, "results", "r5", "logs")
for d in (OUT_TABLES, OUT_FIGS, OUT_LOGS):
    os.makedirs(d, exist_ok=True)

STUDIES_GZ = os.path.join(DATA_DIR, "ewascatalog-studies.txt.gz")
RESULTS_GZ = os.path.join(DATA_DIR, "ewascatalog-results.txt.gz")

P_THRESH = 1e-7
N_PERM = 1000
SEED = 42
CHUNK = 500_000

# trait classifiers -------------------------------------------------------
AGE_RE = re.compile(r"\bage\b|aging|ageing|biological age|epigenetic age", re.I)
GEST_RE = re.compile(r"gestational", re.I)                      # excluded
ALLERGY_RE = re.compile(
    r"asthma|allerg|atopic|atopy|rhinitis|hay ?fever|\bIgE\b|eczema|urticaria|wheeze",
    re.I,
)
BLOOD_RE = re.compile(r"whole blood|peripheral blood|leukocy|leucocy|pbmc|\bblood\b", re.I)

# type-2 / alarmin gene panel ---------------------------------------------
TYPE2_GENES = ["IL33", "TSLP", "IL4R", "IL13", "GATA3", "IL5RA", "IL4", "IL5",
               "IL1RL1", "CRLF2", "IL17RB", "STAT6", "FOXP3", "RORC"]
TYPE2_SET = set(TYPE2_GENES)
TYPE2_RE = re.compile(r"(^|;)(" + "|".join(TYPE2_GENES) + r")(;|$)")

# 450k manifest for gene-mapping fallback (CpG -> UCSC_RefGene_Name)
MANIFEST = os.path.join(BASE, "data", "real", "geo_methylation",
                        "450k_manifest_compact.csv")

LOG = {"steps": []}


def log(msg):
    print(msg, flush=True)
    LOG["steps"].append(msg)


def classify_studies():
    """Map study_id -> (trait, tissue); derive age / allergy / blood study sets."""
    log("Loading studies table ...")
    st = pd.read_csv(STUDIES_GZ, sep="\t", compression="gzip",
                     usecols=["study_id", "trait", "tissue"], dtype=str,
                     low_memory=False)
    st["trait"] = st["trait"].fillna("")
    st["tissue"] = st["tissue"].fillna("")

    def is_age(t):
        return bool(AGE_RE.search(t)) and not bool(GEST_RE.search(t))

    st["is_age"] = st["trait"].apply(is_age)
    st["is_allergy"] = st["trait"].apply(lambda t: bool(ALLERGY_RE.search(t)))
    st["is_blood"] = st["tissue"].apply(lambda t: bool(BLOOD_RE.search(t)))

    age_ids = set(st.loc[st["is_age"], "study_id"])
    allergy_ids = set(st.loc[st["is_allergy"], "study_id"])
    log(f"studies total={len(st)}  age={st['is_age'].sum()}  "
        f"allergy={st['is_allergy'].sum()}  blood={st['is_blood'].sum()}")
    overlap_ids = age_ids & allergy_ids
    if overlap_ids:
        log(f"WARN: {len(overlap_ids)} studies match both age and allergy "
            f"(kept in both classes)")
    st.to_csv(os.path.join(OUT_TABLES, "r5_ewas_study_classification.csv"), index=False)
    return st, age_ids, allergy_ids


def stream_results(st, age_ids, allergy_ids):
    """Single streaming pass over the results table.

    Collects:
      all_cpgs      - every CpG in the catalog (background)
      blood_cpgs    - CpGs reported in blood-tissue studies
      age_cpg/gene  - CpGs with p<1e-7 in age studies (all / blood-only)
      allergy_cpg   - CpGs with p<1e-7 in allergy studies (all / blood-only)
      bg_type2      - background CpGs whose gene annotation hits the type-2 panel
    """
    blood_ids = set(st.loc[st["is_blood"], "study_id"])
    age_blood_ids = age_ids & blood_ids
    allergy_blood_ids = allergy_ids & blood_ids
    target_ids = age_ids | allergy_ids

    all_cpgs, blood_cpgs, bg_type2 = set(), set(), set()
    age_cpg, age_cpg_b = {}, {}
    al_cpg, al_cpg_b = {}, {}

    log("Streaming results table (this takes a few minutes) ...")
    n_rows = 0
    reader = pd.read_csv(
        RESULTS_GZ, sep="\t", compression="gzip", chunksize=CHUNK,
        usecols=[0, 3, 5, 10], names=["cpg", "p", "study_id", "gene"],
        header=0, dtype={"cpg": str, "study_id": str, "gene": str},
        low_memory=False,
    )
    for chunk in reader:
        n_rows += len(chunk)
        chunk["p"] = pd.to_numeric(chunk["p"], errors="coerce")
        all_cpgs.update(chunk["cpg"].dropna().tolist())

        is_blood = chunk["study_id"].isin(blood_ids)
        if is_blood.any():
            blood_cpgs.update(chunk.loc[is_blood, "cpg"].dropna().tolist())

        # type-2 membership for the whole background (vectorised)
        genes = chunk["gene"].fillna("")
        t2_mask = genes.str.contains(TYPE2_RE, na=False)
        if t2_mask.any():
            bg_type2.update(chunk.loc[t2_mask, "cpg"].tolist())

        sub = chunk[chunk["study_id"].isin(target_ids) & (chunk["p"] < P_THRESH)]
        for cpg, sid, gene in zip(sub["cpg"], sub["study_id"], sub["gene"]):
            if sid in age_ids:
                age_cpg.setdefault(cpg, gene)
                if sid in age_blood_ids:
                    age_cpg_b.setdefault(cpg, gene)
            if sid in allergy_ids:
                al_cpg.setdefault(cpg, gene)
                if sid in allergy_blood_ids:
                    al_cpg_b.setdefault(cpg, gene)

        if n_rows % 5_000_000 < CHUNK:
            log(f"  ... {n_rows/1e6:.1f}M rows | bg={len(all_cpgs)} "
                f"age={len(age_cpg)} allergy={len(al_cpg)}")

    log(f"results rows={n_rows}")
    log(f"|S_age|={len(age_cpg)}  |S_allergy|={len(al_cpg)}  |bg|={len(all_cpgs)}")
    log(f"|S_age_blood|={len(age_cpg_b)}  |S_allergy_blood|={len(al_cpg_b)}  "
        f"|bg_blood|={len(blood_cpgs)}")
    return {
        "all_cpgs": all_cpgs, "blood_cpgs": blood_cpgs, "bg_type2": bg_type2,
        "age_cpg": age_cpg, "age_cpg_b": age_cpg_b,
        "al_cpg": al_cpg, "al_cpg_b": al_cpg_b,
        "n_rows": n_rows,
    }


def overlap_test(label, set_a, set_b, background, rng):
    """Permutation + Fisher test of overlap between set_a and set_b."""
    set_a, set_b = set(set_a), set(set_b)
    background = set(background)
    # restrict sets to background (safety)
    set_a &= background
    set_b &= background
    n_a, n_b, n_bg = len(set_a), len(set_b), len(background)
    obs = len(set_a & set_b)

    null = np.zeros(N_PERM, dtype=int)
    if n_a > 0 and n_b > 0 and n_bg > 0:
        # vectorised: boolean mask of set_b over the background array
        bg_arr = np.array(sorted(background))
        b_mask = np.isin(bg_arr, list(set_b))
        for i in range(N_PERM):
            idx = rng.choice(n_bg, size=n_a, replace=False)
            null[i] = int(b_mask[idx].sum())
        perm_p = (1 + int((null >= obs).sum())) / (1 + N_PERM)
    else:
        perm_p = np.nan

    a, b = obs, n_b - obs
    c, d = n_a - obs, n_bg - n_a - n_b + obs
    if min(a, b, c, d) >= 0 and n_bg > 0:
        or_, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    else:
        or_, fisher_p = np.nan, np.nan
    expected = n_a * n_b / n_bg if n_bg else np.nan
    log(f"[{label}] n_age={n_a} n_allergy={n_b} bg={n_bg} overlap={obs} "
        f"expected={expected:.1f} perm_p={perm_p} fisher_p={fisher_p:.3e}")
    return {
        "analysis": label, "n_age": n_a, "n_allergy": n_b, "n_background": n_bg,
        "overlap": obs, "expected_overlap": round(expected, 3),
        "perm_p": perm_p, "perm_n": N_PERM,
        "fisher_p": fisher_p, "fisher_or": or_,
        "null": null,
    }


def type2_enrichment(label, overlap_cpgs, cpg2gene, background, bg_type2):
    """Fisher test: are overlap CpGs enriched for type-2 / alarmin genes?"""
    overlap_cpgs = set(overlap_cpgs)
    t2_in_overlap = {c for c in overlap_cpgs
                     if any(x.strip().upper() in TYPE2_SET
                            for x in str(cpg2gene.get(c, "")).split(";"))}
    a = len(t2_in_overlap)
    b = len(overlap_cpgs) - a
    c = len(bg_type2) - a
    d = len(background) - len(bg_type2) - b
    if min(a, b, c, d) >= 0 and (a + b) > 0:
        or_, p = stats.fisher_exact([[a, b], [c, d]])
    else:
        or_, p = np.nan, np.nan
    genes_hit = sorted({x.strip().upper() for c in t2_in_overlap
                        for x in str(cpg2gene.get(c, "")).split(";")
                        if x.strip().upper() in TYPE2_SET})
    log(f"[{label}] type2 CpGs in overlap={a}/{len(overlap_cpgs)} "
        f"bg_type2={len(bg_type2)} fisher_p={p:.3e} genes={genes_hit}")
    return {
        "analysis": label + " | type2 enrichment",
        "n_age": np.nan, "n_allergy": np.nan, "n_background": len(background),
        "overlap": len(overlap_cpgs),
        "expected_overlap": np.nan,
        "perm_p": np.nan, "perm_n": 0,
        "fisher_p": p, "fisher_or": or_,
        "type2_cpgs_in_overlap": a,
        "type2_bg_cpgs": len(bg_type2),
        "type2_genes_hit": ";".join(genes_hit),
    }


def plot_null(res_all, res_blood, out_path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, res, ttl in zip(axes, [res_all, res_blood],
                            ["All tissues", "Whole blood / leukocyte studies"]):
        null = res["null"]
        obs = res["overlap"]
        ax.hist(null, bins=min(40, max(10, len(np.unique(null)))),
                color="#9ecae1", edgecolor="#3182bd", alpha=0.85,
                label="Permutation null (n=1000)")
        ax.axvline(obs, color="#de2d26", lw=2.2, label=f"Observed overlap = {obs}")
        ax.set_xlabel("Overlap with allergy CpG set (permuted age sets)")
        ax.set_ylabel("Frequency")
        p_txt = (f"perm p = {res['perm_p']:.3g}"
                 if not np.isnan(res["perm_p"]) else "perm p = NA")
        ax.set_title(f"{ttl}\n{p_txt}; Fisher p = {res['fisher_p']:.2e}")
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("R5-10: EWAS convergence — aging vs allergy CpG sets (EWAS Catalog, p<1e-7)",
                 fontsize=11)
    fig.tight_layout()
    save_fig(fig, out_path)


def load_manifest_map():
    """CpG -> gene symbol map from the 450k manifest (fallback annotation)."""
    try:
        m = pd.read_csv(MANIFEST, usecols=["IlmnID", "UCSC_RefGene_Name"],
                        dtype=str)
        return dict(zip(m["IlmnID"], m["UCSC_RefGene_Name"].fillna("")))
    except Exception as e:
        log(f"WARN manifest load failed ({e}); gene fallback disabled")
        return {}


def main():
    rng = np.random.default_rng(SEED)
    try:
        manifest_map = load_manifest_map()
        st, age_ids, allergy_ids = classify_studies()
        data = stream_results(st, age_ids, allergy_ids)

        rows = []

        # ---- primary analysis: all tissues ------------------------------
        res_all = overlap_test("all_tissues", data["age_cpg"], data["al_cpg"],
                               data["all_cpgs"], rng)
        rows.append({k: v for k, v in res_all.items() if k != "null"})

        overlap_all = set(data["age_cpg"]) & set(data["al_cpg"])
        cpg2gene = {}
        cpg2gene.update({c: g for c, g in data["al_cpg"].items()})
        cpg2gene.update({c: g for c, g in data["age_cpg"].items()})
        rows.append(type2_enrichment("all_tissues", overlap_all, cpg2gene,
                                     data["all_cpgs"], data["bg_type2"]))

        # ---- stratified: whole blood / leukocyte studies ----------------
        res_blood = overlap_test("blood_stratified", data["age_cpg_b"], data["al_cpg_b"],
                                 data["blood_cpgs"], rng)
        rows.append({k: v for k, v in res_blood.items() if k != "null"})

        overlap_blood = set(data["age_cpg_b"]) & set(data["al_cpg_b"])
        cpg2gene_b = {}
        cpg2gene_b.update({c: g for c, g in data["al_cpg_b"].items()})
        cpg2gene_b.update({c: g for c, g in data["age_cpg_b"].items()})
        bg_type2_blood = data["bg_type2"] & data["blood_cpgs"]
        rows.append(type2_enrichment("blood_stratified", overlap_blood, cpg2gene_b,
                                     data["blood_cpgs"], bg_type2_blood))

        # ---- write tables ------------------------------------------------
        enrich = pd.DataFrame(rows)
        enrich.to_csv(os.path.join(OUT_TABLES, "r5_ewas_overlap_enrichment.csv"),
                      index=False)

        gene_rows = []
        for scope, ov, gmap in [("all_tissues", overlap_all, cpg2gene),
                                ("blood_stratified", overlap_blood, cpg2gene_b)]:
            for cpg in sorted(ov):
                gene = str(gmap.get(cpg, ""))
                if not gene or gene == "nan":
                    gene = str(manifest_map.get(cpg, ""))  # 450k fallback
                genes_split = [x.strip().upper() for x in gene.split(";") if x.strip()]
                gene_rows.append({
                    "scope": scope, "cpg": cpg, "gene": gene,
                    "in_type2_panel": any(g in TYPE2_SET for g in genes_split),
                    "type2_match": ";".join(g for g in genes_split if g in TYPE2_SET),
                })
        gene_df = pd.DataFrame(gene_rows)
        gene_df.to_csv(os.path.join(OUT_TABLES, "r5_ewas_overlap_genes.csv"),
                       index=False)
        # primary deliverable: all-tissues overlap CpG x gene table
        conv = gene_df[gene_df["scope"] == "all_tissues"].drop(columns="scope")
        conv.to_csv(os.path.join(OUT_TABLES, "r5_ewas_convergence.csv"),
                    index=False)
        log(f"r5_ewas_convergence.csv: {len(conv)} overlap CpGs, "
            f"{int(conv['in_type2_panel'].sum())} in type-2/alarmin panel")

        # ---- figure ------------------------------------------------------
        plot_null(res_all, res_blood,
                  os.path.join(OUT_FIGS, "r5_ewas_convergence.png"))

        LOG["summary"] = {
            "n_age_studies": int(st["is_age"].sum()),
            "n_allergy_studies": int(st["is_allergy"].sum()),
            "S_age": len(data["age_cpg"]), "S_allergy": len(data["al_cpg"]),
            "S_age_blood": len(data["age_cpg_b"]),
            "S_allergy_blood": len(data["al_cpg_b"]),
            "overlap_all": len(overlap_all), "overlap_blood": len(overlap_blood),
            "perm_p_all": res_all["perm_p"], "perm_p_blood": res_blood["perm_p"],
        }
        with open(os.path.join(OUT_LOGS, "r5_10_summary.json"), "w") as f:
            json.dump(LOG, f, indent=2, default=str)
        log("DONE R5-10")
    except Exception as e:  # never crash the pipeline silently
        import traceback
        log(f"FAILED R5-10: {e}")
        with open(os.path.join(OUT_LOGS, "r5_10_summary.json"), "w") as f:
            json.dump({"failed": str(e),
                       "trace": traceback.format_exc(), "steps": LOG["steps"]},
                      f, indent=2)
        raise


if __name__ == "__main__":
    main()
