"""
R5-01: 正交衰老工具变量 MR —— 裁决 R4 工具变量循环性危机
================================================================
暴露：与 R3 血液指数构建完全无关的衰老生物学 GWAS（正交工具）
  - GrimAge 加速度   (GCST90014288, McCartney 2021, n=34,467)
  - Hannum 年龄加速度 (GCST90014289, n=34,449)
  - IEAA             (GCST90014290, n=34,461)
  - PhenoAge 加速度  (GCST90014292, n=34,463)
  - 衰弱指数          (GCST90020053, Atkins 2021, n=175,226)
  （端粒长度 ieu-b-4879 因 OpenGWAS JWT 限制无法获取，记为局限）
结局：FinnGen R12 三过敏表型 (ALLERG_ASTHMA / L12_ATOPIC / ALLERG_RHINITIS)

方法：p<5e-8 选点 -> 1KG EUR clump (r2<0.05, 10Mb) -> harmonise -> Steiger 过滤
      -> IVW(fixed/random) + Egger + weighted median + PRESSO + Cochran Q
多重校正：5 工具 × 3 结局 = 15 个 IVW 检验族内 BH-FDR
裁决（预注册 docs/R5/R5_PREREGISTRATION_UPDATE.md）：
  路径 A: >=2 个独立衰老工具对同一结局方向一致且 FDR<0.05
  路径 B: 全部衰老工具 FDR>=0.05（血液指数显著性交由 R5-02/04）
  路径 C: 其余 -> 收窄命题
"""
import os
import sys
import gzip
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def load_finngen_snps(path, rsids):
    """内存安全的 FinnGen 行扫描：仅保留 rsid 集合内的行（避免 pandas 分块 OOM)。"""
    rsids = set(rsids)
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        need = ["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]
        miss = [c for c in need if c not in idx]
        if miss:
            raise ValueError(f"FinnGen 文件缺列 {miss}: {path}")
        for line in f:
            p = line.rstrip("\n").split("\t")
            rs = p[idx["rsids"]]
            hit = rs in rsids or any(part in rsids for part in rs.split(","))
            if not hit:
                continue
            try:
                rows.append({
                    "chrom": p[idx["#chrom"]], "pos": int(p[idx["pos"]]),
                    "ref": p[idx["ref"]], "alt": p[idx["alt"]], "snp": rs,
                    "p": float(p[idx["pval"]]), "beta": float(p[idx["beta"]]),
                    "se": float(p[idx["sebeta"]]), "af_alt": float(p[idx["af_alt"]]),
                })
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates("snp", keep="first").reset_index(drop=True)
    return df

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(BASE, "src", "r1"))
from utils import (  # noqa: E402
    load_finngen, load_finngen_sample_sizes, harmonise_pair, ld_clump,
    steiger_filter, mr_ivw, mr_egger, weighted_median, mr_presso_outliers,
    cochran_q, FINNGEN_FULL_OUTCOMES, DATA_DIR,
)

R5_DATA = os.path.join(BASE, "data", "r5")
OUT_DIR = os.path.join(BASE, "results", "r5", "tables")
os.makedirs(OUT_DIR, exist_ok=True)

P_SELECT = 5e-8
CLUMP_R2 = 0.05
CLUMP_KB = 10000

# 暴露配置：列名映射 -> 统一 (snp, ea, oa, beta, se, p, eaf, n)
AGING = {
    "GrimAge": dict(file="GrimAge_hits_p1e5.tsv", n=34467,
                    ea="A1", oa="A2", beta="Effect", se="SE", p="p_value",
                    eaf="Freq1", n_col="N"),
    "Hannum": dict(file="Hannum_hits_p1e5.tsv", n=34449,
                   ea="A1", oa="A2", beta="Effect", se="SE", p="p_value",
                   eaf="Freq1", n_col="N"),
    "IEAA": dict(file="IEAA_hits_p1e5.tsv", n=34461,
                 ea="A1", oa="A2", beta="Effect", se="SE", p="p_value",
                 eaf="Freq1", n_col="N"),
    "PhenoAge": dict(file="PhenoAge_hits_p1e5.tsv", n=34463,
                     ea="A1", oa="A2", beta="Effect", se="SE", p="p_value",
                     eaf="Freq1", n_col="N"),
    "Frailty": dict(file="Frailty_hits_p1e5.tsv", n=175226,
                    ea="effect_allele", oa="other_allele", beta="beta",
                    se="standard_error", p="p_value",
                    eaf="effect_allele_frequency", n_col=None),
}

OUTCOMES = ["ALLERG_ASTHMA", "L12_ATOPIC", "ALLERG_RHINITIS"]


def load_exposure(name, cfg):
    path = os.path.join(R5_DATA, cfg["file"])
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        print(f"  !! {name} 子集文件缺失或为空 ({path})，跳过")
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:
        print(f"  !! {name} 读取失败: {e}，跳过")
        return pd.DataFrame()
    if df.empty:
        return df
    out = pd.DataFrame({
        "snp": df["variant_id"].astype(str),
        "ea": df[cfg["ea"]].astype(str).str.upper(),
        "oa": df[cfg["oa"]].astype(str).str.upper(),
        "beta": pd.to_numeric(df[cfg["beta"]], errors="coerce"),
        "se": pd.to_numeric(df[cfg["se"]], errors="coerce"),
        "p": pd.to_numeric(df[cfg["p"]], errors="coerce"),
        "eaf": pd.to_numeric(df[cfg["eaf"]], errors="coerce"),
        "n_snp": (pd.to_numeric(df[cfg["n_col"]], errors="coerce")
                  if cfg["n_col"] else cfg["n"]),
    })
    out = out.dropna(subset=["snp", "ea", "oa", "beta", "se", "p"])
    out = out[out["snp"].str.startswith("rs")]
    out = out.drop_duplicates("snp", keep="first")
    sel = out[out["p"] < P_SELECT].copy()
    print(f"  {name}: {len(out)} 行 p<1e-5 子集 -> {len(sel)} 个 p<5e-8 SNP")
    return sel


def run_mr_for_pair(har, n_exp, n_out):
    """对调和后的 SNP 表跑全套 MR，返回 dict。"""
    res = {}
    har = steiger_filter(har, n_exp, n_out)
    res["n_steiger_removed"] = int((~har["steiger_pass"]).sum())
    har = har[har["steiger_pass"]].reset_index(drop=True)
    res["n_snp"] = len(har)
    if len(har) < 3:
        res["note"] = "insufficient SNPs (<3)"
        return res, har
    bx = har["beta"].values
    by = har["beta_outcome"].values
    sy = har["se_outcome"].values
    res["mean_F"] = float(np.mean((bx / har["se"].values) ** 2))
    b, se, p = mr_ivw(bx, by, sy, random=False)
    res.update(ivw_beta=b, ivw_se=se, ivw_p=p, ivw_or=float(np.exp(b)))
    br, ser, pr = mr_ivw(bx, by, sy, random=True)
    res.update(ivw_re_beta=br, ivw_re_p=pr)
    Q, qp = cochran_q(bx, by, sy, b)
    res.update(cochran_q=Q, cochran_q_p=qp)
    try:
        es, ese, esp, ei, eip = mr_egger(bx, by, sy)
        res.update(egger_beta=es, egger_p=esp, egger_intercept=ei, egger_int_p=eip)
    except Exception as e:
        res["egger_err"] = str(e)
    try:
        wm, wms, wmp = weighted_median(bx, by, sy)
        res.update(wmed_beta=wm, wmed_p=wmp)
    except Exception as e:
        res["wmed_err"] = str(e)
    try:
        keep, iters = mr_presso_outliers(bx, by, sy)
        res["n_presso_outliers"] = int((~keep).sum())
        if keep.sum() >= 3:
            bp, sep, pp = mr_ivw(bx[keep], by[keep], sy[keep])
            res.update(presso_beta=bp, presso_p=pp)
    except Exception as e:
        res["presso_err"] = str(e)
    return res, har


def main():
    print("=" * 70)
    print("R5-01 正交衰老工具 MR")
    print("=" * 70)

    # 1) 暴露：选点 + clump
    instruments = {}
    for name, cfg in AGING.items():
        sel = load_exposure(name, cfg)
        if sel.empty:
            print(f"  !! {name} 无 p<5e-8 SNP，跳过")
            continue
        cl = ld_clump(sel, r2_thresh=CLUMP_R2, kb=CLUMP_KB)
        print(f"  {name}: clump 后 {len(cl)} 个独立工具 SNP")
        instruments[name] = cl

    all_snps = set()
    for cl in instruments.values():
        all_snps.update(cl["snp"].tolist())
    print(f"\n  5 工具合计独立 SNP（并集）: {len(all_snps)}")

    cl_rows = []
    for name, cl in instruments.items():
        for _, r in cl.iterrows():
            cl_rows.append(dict(trait=name, snp=r["snp"], beta=r["beta"],
                                se=r["se"], p=r["p"], eaf=r["eaf"], n_snp=r["n_snp"]))
    pd.DataFrame(cl_rows).to_csv(os.path.join(OUT_DIR, "r5_aging_instruments_clumped.csv"), index=False)

    # 2) 结局：FinnGen 三表型，仅取并集 SNP
    n_map = load_finngen_sample_sizes()
    outcome_data = {}
    for oc in OUTCOMES:
        path = FINNGEN_FULL_OUTCOMES[oc]
        df = load_finngen_snps(path, all_snps)
        print(f"  结局 {oc}: 命中 {len(df)}/{len(all_snps)} SNP, N={n_map.get(oc)}", flush=True)
        outcome_data[oc] = df

    # 3) 逐 (工具, 结局) MR
    rows = []
    for name, cl in instruments.items():
        n_exp = int(np.nanmedian(cl["n_snp"])) if len(cl) else AGING[name]["n"]
        for oc in OUTCOMES:
            har = harmonise_pair(cl, outcome_data[oc])
            n_out = n_map.get(oc, 0)
            if har.empty:
                rows.append(dict(trait=name, outcome=oc, n_snp=0, note="no harmonised SNP"))
                continue
            res, har2 = run_mr_for_pair(har, n_exp, n_out)
            res.update(trait=name, outcome=oc, n_out=n_out)
            rows.append(res)
            tag = f"{name:9s} x {oc:15s}"
            if "ivw_p" in res:
                print(f"  {tag} n={res['n_snp']:3d} F={res.get('mean_F', np.nan):7.1f} "
                      f"IVW beta={res['ivw_beta']:+.4f} p={res['ivw_p']:.2e} "
                      f"OR={res['ivw_or']:.3f} Q_p={res.get('cochran_q_p', np.nan):.2e} "
                      f"EggerInt_p={res.get('egger_int_p', np.nan):.3f}")
            else:
                print(f"  {tag} {res.get('note', 'NA')}")

    res_df = pd.DataFrame(rows)

    # 4) 族内 BH-FDR（15 检验）
    mask = res_df["ivw_p"].notna()
    if mask.sum() > 0:
        rej, q, _, _ = multipletests(res_df.loc[mask, "ivw_p"], method="fdr_bh")
        res_df.loc[mask, "ivw_q"] = q
        res_df.loc[mask, "ivw_fdr05"] = rej
    res_df.to_csv(os.path.join(OUT_DIR, "r5_orthogonal_aging_mr.csv"), index=False)

    # 5) 决策树裁决（方向一致性严格执行：显著工具中不得有反向显著者才算"一致")
    verdict = {}
    for oc in OUTCOMES:
        sub = res_df[(res_df["outcome"] == oc) & (res_df["ivw_fdr05"] == True)]
        n_sig = len(sub)
        signs = list(np.sign(sub["ivw_beta"].values)) if n_sig else []
        n_up = int(sum(1 for s in signs if s > 0))
        n_dn = int(sum(1 for s in signs if s < 0))
        consistent = n_sig >= 2 and (n_up == n_sig or n_dn == n_sig)
        verdict[oc] = dict(n_sig=n_sig, n_up=n_up, n_down=n_dn,
                           sig_traits=sub["trait"].tolist(),
                           sig_betas=[float(b) for b in sub["ivw_beta"].values],
                           direction_consistent=consistent)
    n_sig_total = int((res_df["ivw_fdr05"] == True).sum()) if mask.sum() else 0
    path_a_outcomes = [oc for oc, v in verdict.items()
                       if isinstance(v, dict) and v["n_sig"] >= 2 and v["direction_consistent"]]
    if path_a_outcomes:
        path = "A"
    elif n_sig_total == 0:
        path = "B_candidate"  # 需 R5-02/04 血液指数显著性确认
    else:
        path = "C"
    verdict["summary"] = dict(n_sig_total=n_sig_total, path_a_outcomes=path_a_outcomes,
                              decision_path=path)
    with open(os.path.join(BASE, "results", "r5", "r5_01_verdict.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 70)
    print(f"决策树裁决：路径 {path}")
    for oc in OUTCOMES:
        v = verdict[oc]
        print(f"  {oc}: 显著工具 {v['n_sig']} 个(正{v['n_up']}/负{v['n_down']}) "
              f"{v['sig_traits']} 方向一致={v['direction_consistent']}")
    print(f"  路径 A 结局: {path_a_outcomes}")
    print("=" * 70)


if __name__ == "__main__":
    main()
