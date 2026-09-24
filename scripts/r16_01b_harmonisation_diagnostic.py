#!/usr/bin/env python3
"""Diagnose r3-vs-r16 discrepancy for the eosinophil-free arm (asthma)."""
import gzip, os, re
from io import StringIO
import numpy as np, pandas as pd
from scipy import stats

BASE = r"D:\衰老研究\v3_pipeline"
FIN = os.path.join(BASE, "data", "real", "finngen_full")
IDX = os.path.join(BASE, "results", "r3", "tables", "r3_index_instruments.csv")
FG_COLS = ["chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp",
           "beta","se","af_alt","af_alt_cases","af_alt_controls"]

d = pd.read_csv(IDX).drop_duplicates("snp")
d = d[d.component != "eosinophil_pct"].copy()
for c in ("ea","oa"): d[c] = d[c].str.upper()
rs = set(d.snp)
pat = re.compile("|".join(re.escape(x) for x in sorted(rs)))
hits = []
with gzip.open(os.path.join(FIN, "finngen_R12_ALLERG_ASTHMA.gz"), "rt",
               errors="replace") as fh:
    fh.readline(); buf = ""
    while True:
        block = fh.read(64*1024*1024)
        if not block: break
        buf += block
        cut = buf.rfind("\n")
        if cut == -1: continue
        seg, buf = buf[:cut], buf[cut+1:]
        if pat.search(seg):
            hits.extend(ln for ln in seg.split("\n") if pat.search(ln))
    if buf and pat.search(buf): hits.append(buf)
raw = pd.read_csv(StringIO("\n".join(hits)), sep="\t", header=None,
                  names=FG_COLS, low_memory=False)
tok = raw["rsids"].astype(str).str.split(",")
exact = tok.apply(lambda ts: any(t in rs for t in ts))
raw = raw[exact].copy()
raw["snp"] = tok[exact].apply(lambda ts: next(t for t in ts if t in rs))
out = raw.drop_duplicates("snp")
m = d.merge(out, on="snp", how="inner", suffixes=("_inst","_fg"))
print("matched:", len(m), "of", len(d))
m["alt_fg"] = m["alt_fg"].str.upper(); m["ref_fg"] = m["ref_fg"].str.upper()
pal = (m.ea.isin(["A","T"]) & m.oa.isin(["A","T"])) | \
      (m.ea.isin(["C","G"]) & m.oa.isin(["C","G"]))
same = m.ea == m.alt_fg
flip = m.ea == m.ref_fg
m["allele_ok"] = same | flip
m["pal"] = pal
m["amb_fg"] = pal & m["af_alt_fg"].between(0.42,0.58)
print("dropped allele-mismatch:", int((~m.allele_ok).sum()))
print("palindromic:", int(pal.sum()), " ambiguous(FG af 0.42-0.58):", int(m["amb_fg"].sum()))
m["by"] = np.where(same, m["beta_fg"], -m["beta_fg"])
m["wr"] = m["by"]/m["beta_I"]; m["wr_se"] = m["se_fg"]/m["beta_I"].abs()

def ivw(sub):
    w = 1/sub["wr_se"]**2
    mu = (w*sub["wr"]).sum()/w.sum()
    se = np.sqrt(1/w.sum())
    return len(sub), round(float(np.exp(mu)),3), float(2*stats.norm.sf(abs(mu/se)))

print("rule A (r16: drop amb by FG af):", ivw(m[m.allele_ok & ~m.amb_fg]))
b = m[m.allele_ok].copy()
# instrument af_alt_inst = freq of instrument 'alt' allele; ea freq on instrument side:
ea_af_i = np.where(b.ea == b.alt_inst.str.upper(), b.af_alt_inst, 1-b.af_alt_inst)
ea_af_fg = np.where(same[b.index], b.af_alt_fg, 1-b.af_alt_fg)
flip_freq = np.abs(ea_af_i-(1-ea_af_fg)) < np.abs(ea_af_i-ea_af_fg)
b.loc[flip_freq & b.pal, "by"] = -b.loc[flip_freq & b.pal, "by"]
b["wr"] = b["by"]/b["beta_I"]; b["wr_se"] = b["se_fg"]/b["beta_I"].abs()
print("rule B (keep pal, freq-orient):", ivw(b))
print("rule C (keep pal as-is):", ivw(m[m.allele_ok]))
amb = m[m.allele_ok & m.amb_fg]
print("\nambiguous-dropped SNPs (n=%d):" % len(amb))
print(amb[["snp","component","beta_I","wr","af_alt_fg","af_alt_inst"]].to_string(index=False))
mm = m[~m.allele_ok]
print("\nallele-mismatch SNPs (n=%d):" % len(mm))
print(mm[["snp","component","ea","oa","ref_fg","alt_fg"]].to_string(index=False))
