# -*- coding: utf-8 -*-
"""Extract canonical LDSC (bulik/ldsc v1.0.1) results from logs and build
Supplementary Tables ST137 (canonical rg) and ST138 (canonical rg contrast).

Contrast SE conventions replicate ST134:
  rho = 0   : se_diff = sqrt(se_f^2 + se_c^2)
  rho = 0.5 : se_diff = sqrt(se_f^2 + se_c^2 - se_f*se_c)
q values: Benjamini-Hochberg within the 12-test canonical family.
"""
import csv
import math
import os
import re
from statistics import median

BASE = r"D:\衰老研究\v3_pipeline\data\r6\canonical_ldsc"
LOGS = os.path.join(BASE, "logs")
OUTPKG = r"D:\衰老研究\v3_pipeline\Genome Medicine\04_补充数据_supplementary_data"

EXPOSURES = ["GrimAge", "Hannum", "IEAA", "PhenoAge", "Frailty"]
OUTCOMES = ["ALLERG_ASTHMA", "ALLERG_RHINITIS", "L12_ATOPIC"]


def parse_rg_log(path):
    """Return (h2_trait1, pairs) where pairs = list of dicts for trait1 vs each subsequent trait."""
    txt = open(path, encoding="utf-8", errors="replace").read()
    h2_1 = re.search(r"Heritability of phenotype 1.*?Total Observed scale h2:\s*([-\d.e]+)\s*\(([\d.e]+)\)", txt, re.S)
    h2_trait1 = (float(h2_1.group(1)), float(h2_1.group(2))) if h2_1 else (None, None)
    # per phenotype-2 blocks
    blocks = re.split(r"Computing rg for phenotype (\d+)/\d+", txt)[1:]
    pairs = []
    for i in range(0, len(blocks), 2):
        body = blocks[i + 1]
        # NB: the phenotype-1 h2 section is restated only inside the FIRST block;
        # robust rule: h2/int lines before the gencov section belong to phenotype 2.
        pre, _, post = body.partition("Total Observed scale gencov")
        g = re.findall(r"Total Observed scale h2:\s*([-\d.e]+)\s*\(([\d.e]+)\)", pre)
        gc = re.search(r"Total Observed scale gencov:\s*([-\d.e]+)\s*\(([\d.e]+)\)", body)
        gr = re.search(r"Genetic Correlation:\s*([-\d.e]+)\s*\(([\d.e]+)\)", body)
        z = re.search(r"Z-score:\s*([-\d.e]+)", body)
        p = re.search(r"^P:\s*([-\d.e]+)", body, re.M)
        nvalid = re.search(r"(\d+) SNPs with valid alleles", body)
        h2t2 = g[-1] if g else None
        ints_pre = re.findall(r"Intercept:\s*([-\d.e]+)\s*\(([\d.e]+)\)", pre)
        ints_post = re.findall(r"Intercept:\s*([-\d.e]+)\s*\(([\d.e]+)\)", post)
        pairs.append(dict(
            h2_trait2=float(h2t2[0]) if h2t2 else None,
            h2_trait2_se=float(h2t2[1]) if h2t2 else None,
            gencov=float(gc.group(1)) if gc else None,
            gencov_se=float(gc.group(2)) if gc else None,
            rg=float(gr.group(1)) if gr else None,
            rg_se=float(gr.group(2)) if gr else None,
            z=float(z.group(1)) if z else None,
            p=float(p.group(1)) if p else None,
            n_valid=int(nvalid.group(1)) if nvalid else None,
            h2_int=float(ints_pre[-1][0]) if ints_pre else None,
            gcov_int=None,
            gcov_int_se=float(ints_post[0][1]) if ints_post else None,
        ))
    # gcov intercept from the summary table at the end
    summ = txt.split("Summary of Genetic Correlation Results")[-1]
    for j, line in enumerate(summ.splitlines()):
        m = re.match(r"\s*\S+\s+(\S+\.sumstats\.gz)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.e]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)", line)
        if m and j < len(pairs):
            pairs[j]["gcov_int"] = float(m.group(10))
            pairs[j]["gcov_int_se"] = float(m.group(11))
    return h2_trait1, pairs


def norm_p(p):
    return p


def main():
    rows = []          # ST137 rows
    h2_trait1 = {}
    rg = {}            # (exp, out) -> row dict (mhc_free)
    for exp in EXPOSURES:
        h1, pairs = parse_rg_log(os.path.join(LOGS, "ldsc_rg_%s.log" % exp))
        h2_trait1[exp] = h1
        for out, pr in zip(OUTCOMES, pairs):
            pr["trait1"], pr["trait2"] = exp, out
            rg[(exp, out)] = pr
            rows.append(dict(trait1=exp, trait2=out, ld_reference="canonical_bulik_ldsc_mhc_free",
                             n_valid=pr["n_valid"], rg=pr["rg"], rg_se=pr["rg_se"], z=pr["z"], p=pr["p"],
                             h2_trait2=pr["h2_trait2"], h2_trait2_se=pr["h2_trait2_se"],
                             h2_int=pr["h2_int"], gcov_int=pr["gcov_int"]))

    # MHC sensitivity (Frailty)
    _, pairs_mhc = parse_rg_log(os.path.join(LOGS, "ldsc_rg_mhc_Frailty.log"))
    for out, pr in zip(OUTCOMES, pairs_mhc):
        rows.append(dict(trait1="Frailty", trait2=out, ld_reference="canonical_bulik_ldsc_mhc_included",
                         n_valid=pr["n_valid"], rg=pr["rg"], rg_se=pr["rg_se"], z=pr["z"], p=pr["p"],
                         h2_trait2=pr["h2_trait2"], h2_trait2_se=pr["h2_trait2_se"],
                         h2_int=pr["h2_int"], gcov_int=pr["gcov_int"]))
        rg[("Frailty", out, "mhc")] = pr

    # FDR within the 15 mhc_free pairs (BH), to mirror ST60's q_bh convention
    ps = sorted((r["p"], i) for i, r in enumerate(rows) if r["ld_reference"].endswith("mhc_free"))
    m = len(ps)
    qvals = {}
    prev = 1.0
    for rank, (p, i) in reversed(list(enumerate(ps, start=1))):
        q = min(prev, p * m / rank)
        qvals[i] = q
        prev = q
    for i, r in enumerate(rows):
        r["q_bh_mhc_free_family"] = qvals.get(i, "")

    # ST137
    st137 = os.path.join(OUTPKG, "ST137_canonical_ldsc_rg.csv")
    with open(st137, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["trait1", "trait2", "ld_reference", "n_valid", "rg", "rg_se", "z", "p",
                                           "q_bh_mhc_free_family", "h2_trait2", "h2_trait2_se", "h2_int", "gcov_int"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote", st137)

    # ST138 contrasts: frailty vs each clock per outcome (canonical numbers)
    def bh(ps):
        m = len(ps)
        order = sorted(range(m), key=lambda i: ps[i])
        q = [0.0] * m
        prev = 1.0
        for rank in range(m, 0, -1):
            i = order[rank - 1]
            prev = min(prev, ps[i] * m / rank)
            q[i] = prev
        return q

    crows = []
    for out in OUTCOMES:
        f = rg[("Frailty", out)]
        for c in ["GrimAge", "Hannum", "IEAA", "PhenoAge"]:
            ck = rg[(c, out)]
            diff = f["rg"] - ck["rg"]
            se0 = math.sqrt(f["rg_se"] ** 2 + ck["rg_se"] ** 2)
            se05 = math.sqrt(f["rg_se"] ** 2 + ck["rg_se"] ** 2 - f["rg_se"] * ck["rg_se"])
            z0, z05 = diff / se0, diff / se05
            p0 = math.erfc(abs(z0) / math.sqrt(2))
            p05 = math.erfc(abs(z05) / math.sqrt(2))
            crows.append(dict(outcome=out, clock=c, source="canonical_bulik_ldsc",
                              rg_frailty=f["rg"], se_frailty=f["rg_se"],
                              rg_clock=ck["rg"], se_clock=ck["rg_se"],
                              clock_rg_hi95=ck["rg"] + 1.96 * ck["rg_se"],
                              frailty_rg_lo95=f["rg"] - 1.96 * f["rg_se"],
                              ci_nonoverlap=(ck["rg"] + 1.96 * ck["rg_se"]) < (f["rg"] - 1.96 * f["rg_se"]),
                              diff=diff, se_diff_rho0=se0, z_rho0=z0, p_rho0=p0,
                              se_diff_rho0p5=se05, p_rho0p5=p05, q_bh_rho0="", q_bh_rho0p5=""))
    q0 = bh([r["p_rho0"] for r in crows])
    q05 = bh([r["p_rho0p5"] for r in crows])
    for i, r in enumerate(crows):
        r["q_bh_rho0"] = q0[i]
        r["q_bh_rho0p5"] = q05[i]
    st138 = os.path.join(OUTPKG, "ST138_canonical_rg_contrast.csv")
    with open(st138, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(crows[0].keys()))
        w.writeheader()
        w.writerows(crows)
    print("wrote", st138)

    # comparison vs in-house ST60
    st60 = os.path.join(OUTPKG, "ST60_r6_ldsc_rg.csv")
    diffs = []
    with open(st60, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = (r["trait1"], r["trait2"])
            if key in rg and r["trait1"] in EXPOSURES and r["trait2"] in OUTCOMES:
                d = abs(rg[key]["rg"] - float(r["rg"]))
                diffs.append(d)
    print("canonical vs in-house |diff| rg: median=%.4f max=%.4f" % (median(diffs), max(diffs)))
    # MHC sensitivity deltas
    for out in OUTCOMES:
        free, incl = rg[("Frailty", out)], rg[("Frailty", out, "mhc")]
        print("Frailty-%s: mhc_free=%.4f mhc_included=%.4f delta=%+.4f" % (out, free["rg"], incl["rg"], incl["rg"] - free["rg"]))


if __name__ == "__main__":
    main()
