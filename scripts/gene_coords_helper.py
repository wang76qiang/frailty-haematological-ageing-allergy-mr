#!/usr/bin/env python3
"""Shared helper: resolve GRCh38 gene coordinates for the INTERVAL pQTL panel.

Priority:
  1. data/real/finngen/target_gene_coords.csv (local, 28 genes)
  2. results/r5/tables/r5_gene_coords_extended.csv (cache from earlier runs)
  3. Ensembl REST lookup (one quick pass, short timeout)
  4. curated static fallback (marked source='static_fallback')
Genes that cannot be resolved are recorded with source='failed' and skipped.
"""

import json
import os
import socket
import time
import urllib.request

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LOCAL_COORDS = os.path.join(BASE_DIR, "data", "real", "finngen", "target_gene_coords.csv")
CACHE_PATH = os.path.join(BASE_DIR, "results", "r5", "tables", "r5_gene_coords_extended.csv")

# protein file prefix -> gene symbol used for coordinate lookup
PROTEIN_TO_GENE = {
    "ADA": "ADA", "ARTN": "ARTN", "CASP8": "CASP8", "CCL2_MCP1": "CCL2",
    "CCL3_MIP1A": "CCL3", "CST5": "CST5", "CX3CL1": "CX3CL1", "GDNF": "GDNF",
    "HGF": "HGF", "IFNG": "IFNG", "IL10": "IL10", "IL13": "IL13",
    "IL17A": "IL17A", "IL17C": "IL17C", "IL18": "IL18", "IL1A": "IL1A",
    "IL20": "IL20", "IL22RA1": "IL22RA1", "IL24": "IL24", "IL2": "IL2",
    "IL33": "IL33", "IL4": "IL4", "IL5": "IL5", "IL6": "IL6", "IL7": "IL7",
    "IL8": "CXCL8", "LAP_TGFB1": "TGFB1", "MMP10": "MMP10", "MMP1": "MMP1",
    "NGF": "NGF", "NRTN": "NRTN", "NTF3": "NTF3", "OPG": "TNFRSF11B",
    "OSM": "OSM", "PDL1": "CD274", "PLAU": "PLAU", "S100A12": "S100A12",
    "SCF": "KITLG", "SLAMF1": "SLAMF1", "TGFA": "TGFA", "TNFB": "LTA",
    "TNFSF12": "TNFSF12", "TNFSF14": "TNFSF14", "TNF": "TNF",
    "TRAIL": "TNFSF10", "TSLP": "TSLP", "VEGFA": "VEGFA",
}

# Curated GRCh38 coordinates (NCBI RefSeq, approximate), used only if Ensembl fails.
STATIC_FALLBACK = {
    "ADA": ("20", 44619532, 44652806, -1),
    "ARTN": ("1", 43937541, 43981898, 1),
    "CASP8": ("2", 201233528, 201287611, -1),
    "CCL2": ("17", 34255274, 34257208, -1),
    "CCL3": ("17", 36085514, 36087568, -1),
    "CST5": ("20", 23804778, 23809899, -1),
    "CX3CL1": ("16", 57373331, 57386769, -1),
    "GDNF": ("5", 37812446, 37838849, -1),
    "HGF": ("7", 81328124, 81399714, -1),
    "IL10": ("1", 206767602, 206772494, -1),
    "IL17A": ("6", 52119849, 52124145, -1),
    "IL17C": ("16", 88663965, 88665215, -1),
    "IL18": ("11", 112143825, 112164501, -1),
    "IL1A": ("2", 112773826, 112784764, -1),
    "IL20": ("1", 105668190, 105672252, 1),
    "IL22RA1": ("1", 243648260, 243664650, -1),
    "IL24": ("1", 105695053, 105702568, -1),
    "IL7": ("8", 78736963, 78808745, 1),
    "CXCL8": ("4", 73740519, 73747379, 1),
    "TGFB1": ("19", 41288203, 41353961, -1),
    "MMP10": ("11", 102640573, 102650645, -1),
    "MMP1": ("11", 102789032, 102797171, 1),
    "NGF": ("1", 115330575, 115382880, -1),
    "NRTN": ("19", 1433883, 1445984, 1),
    "NTF3": ("12", 5454436, 5836893, -1),
    "TNFRSF11B": ("8", 118923557, 118952251, -1),
    "OSM": ("22", 30345485, 30350636, -1),
    "CD274": ("9", 5450503, 5470566, 1),
    "PLAU": ("10", 73936997, 73943440, -1),
    "S100A12": ("1", 153352687, 153354751, 1),
    "KITLG": ("12", 88492789, 88580851, -1),
    "SLAMF1": ("1", 159897860, 159914653, 1),
    "TGFA": ("2", 70442353, 70543378, -1),
    "LTA": ("6", 31572054, 31574324, 1),
    "TNFSF12": ("17", 7483542, 7487853, -1),
    "TNFSF14": ("19", 6661728, 6669075, -1),
    "TNFSF10": ("3", 172505503, 172523529, -1),
    "VEGFA": ("6", 43737918, 43754224, -1),
}


def ensembl_lookup(gene, timeout=12):
    url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
           f"?content-type=application/json")
    socket.setdefaulttimeout(timeout)
    with urllib.request.urlopen(url) as r:
        d = json.load(r)
    return (str(d["seq_region_name"]), int(d["start"]), int(d["end"]),
            int(d["strand"]), d.get("id", ""))


def load_local_coords():
    import pandas as pd
    df = pd.read_csv(LOCAL_COORDS)
    out = {}
    for _, r in df.iterrows():
        out[r["gene_symbol"]] = {
            "gene": r["gene_symbol"], "ensembl_id": r["ensembl_id"],
            "chr": str(r["chr"]), "start": int(r["start"]), "end": int(r["end"]),
            "strand": int(r["strand"]), "source": "local_target_gene_coords",
        }
    return out


def load_cache():
    import pandas as pd
    out = {}
    if os.path.exists(CACHE_PATH):
        df = pd.read_csv(CACHE_PATH)
        for _, r in df.iterrows():
            if r.get("source") == "failed":
                continue
            out[r["gene"]] = {
                "gene": r["gene"], "ensembl_id": r.get("ensembl_id", ""),
                "chr": str(r["chr"]), "start": int(r["start"]), "end": int(r["end"]),
                "strand": int(r["strand"]), "source": r.get("source", "cache"),
            }
    return out


def resolve_panel_coords(log=print):
    """Return dict gene -> coord record for all proteins in PROTEIN_TO_GENE."""
    import pandas as pd
    os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)
    coords = load_local_coords()
    cache = load_cache()
    coords.update({g: c for g, c in cache.items() if g not in coords})

    needed = sorted(set(PROTEIN_TO_GENE.values()))
    missing = [g for g in needed if g not in coords]
    log(f"[coords] local+cache cover {len(needed)-len(missing)}/{len(needed)} genes; "
        f"missing: {missing}")

    fetched, failed = {}, []
    if missing:
        log("[coords] one quick Ensembl pass (short timeout, no long retries) ...")
        for g in missing:
            ok = False
            for attempt in range(2):
                try:
                    chrom, start, end, strand, eid = ensembl_lookup(g)
                    fetched[g] = {"gene": g, "ensembl_id": eid, "chr": chrom,
                                  "start": start, "end": end, "strand": strand,
                                  "source": "ensembl_rest"}
                    ok = True
                    break
                except Exception as e:
                    if attempt == 1:
                        log(f"[coords] Ensembl failed for {g}: {repr(e)[:80]}")
                    time.sleep(0.5)
            if not ok:
                failed.append(g)
    # static fallback for failures
    for g in failed:
        if g in STATIC_FALLBACK:
            chrom, start, end, strand = STATIC_FALLBACK[g]
            fetched[g] = {"gene": g, "ensembl_id": "", "chr": chrom, "start": start,
                          "end": end, "strand": strand, "source": "static_fallback"}
            log(f"[coords] {g}: using static fallback")
        else:
            log(f"[coords] {g}: UNRESOLVED -> protein will be skipped")

    # update cache with everything new (including failures for transparency)
    rows = list(cache.values())
    rows += [v for v in fetched.values()]
    for g in failed:
        if g not in STATIC_FALLBACK:
            rows.append({"gene": g, "ensembl_id": "", "chr": "", "start": -1,
                         "end": -1, "strand": 0, "source": "failed"})
    if rows:
        pd.DataFrame(rows).drop_duplicates("gene", keep="last").to_csv(
            CACHE_PATH, index=False)
        log(f"[coords] cache updated -> {CACHE_PATH}")

    coords.update(fetched)
    return coords


if __name__ == "__main__":
    c = resolve_panel_coords()
    print(f"resolved {len(c)} genes")
