# Nature 子刊补强包：全量 Immune/Lung 单细胞图谱供体级表达-年龄分析（P0-1）

**日期**：2026-09-25 | **数据**：Tabula Sapiens 2.0（`G:\衰老\Immune.h5ad`，592,317 细胞 / 24 供体 / 19 年龄阶段 / 74 组织；`Lung.h5ad`，65,847 细胞 / 3 供体 59–61 岁）| **计数层**：`layers/decontXcounts`（去污染计数）

## 预声明分析计划（分析前锁定，见 scripts/ 时间戳）

- **推断单位 = 供体**（不将全部细胞并入 GLM；细胞级汇总仅在敏感性层，因其把同一供体的细胞当独立重复、P 值反保守）。
- **主推断模型**：等供体加权 WLS，响应 = log(供体表达率 + 0.5/总计数)，协变量 = 供体年龄（中心化）。
- **FDR 族**（预声明）：Treg 族 5 检验（IL4 主检验 + IL13/GATA3/FOXP3/IL5 背景）；转化族 4 检验（IL18×2、TNFSF14×2 细胞区室）。肺图谱仅描述（3 供体同龄段，无法做年龄检验）。
- **环境 RNA 证伪检验**（预声明）：供体级 Treg 目标基因率 vs 浆细胞同基因率（浆细胞不生物表达 IL4/IL5；相关即提示环境污染）。

## 结果（ST139–ST145，figures/）

### Treg 主族（12 供体有 Treg，共 2,649 细胞；ST140 主表 / ST144 敏感性）

| 基因 | 主推断 RR/十年 (95% CI) | P | 族内 q | Spearman rho (P) | LOO 方向 | 环境代理 |
|---|---|---|---|---|---|---|
| IL5 | **1.72 (1.15–2.58)** | **0.025** | 0.125 | 0.57 (0.054) | 稳定 (>1) | rho=0.48 (0.11) |
| IL4 | 1.97 (0.78–5.02) | 0.184 | 0.46 | 0.44 (0.151) | 稳定 (1.57–2.62) | **rho=0.59 (0.045)** |
| IL13 | 1.14 (0.33–3.98) | 0.840 | 0.84 | 0.57 (0.052) | — | — |
| GATA3 | 0.68 (0.30–1.55) | 0.382 | 0.62 | −0.16 (0.61) | — | — |
| FOXP3 | 0.82 (0.46–1.44) | 0.499 | 0.62 | −0.09 (0.78) | — | — |

解读：**IL5 名义显著**（P=0.025，族内 q=0.13，未过 FDR），方向与 LOO 稳定；**IL4 方向一致为正但供体等权后不显著**，且与环境代理相关（部分环境 RNA）。细胞计数加权的合并 Poisson（P~1e-253）为伪重复放大，仅作敏感性层报告（ST140 已标注层级）。**结论分级：提名级正向证据**——支持"2 型表达随年龄漂移上调"，不构成 DICE IL4 共定位的独立复制。

### 转化层（ST141/ST143/ST145）

- **IL18**：血液中免疫细胞几乎不表达（最高 classical monocyte 0.34% 细胞、率 0.22/1e4）；**肺图谱中高度定位于肺泡巨噬细胞（75% 细胞）、肺泡 I 型（71%）、基底细胞（59%）、肺泡 II 型（55%）**——屏障上皮+巨噬细胞来源，与 IL18 作为上皮/髓系靶点的机制一致（ST145）。
- **TSLP 肺泡 II 型 42%**（已知 AT2/棒细胞来源的验证——证明图谱可检测上皮alarmin，间接支持 IL18 定位的可信度）。
- **TNFSF14**：单核吞噬细胞/组织驻留巨噬细胞 96% 细胞表达——LIGHT 的髓系来源提示。
- IL18/TNFSF14 年龄斜率不可估（供体数 3–17 且计数稀疏，ST143）。

### 环境 RNA 与图谱限度（必须随文声明）

IL4 在浆细胞（13%）、红系祖细胞（41%）中"可检测"——非生物来源，提示残余环境污染；decontXcounts 已校正但仍需声明。Lung 图谱 3 供体均 ~60 岁，无年龄变异。

## 文件清单

| 文件 | 内容 |
|---|---|
| `tables/ST139_treg_donor_pseudobulk.csv` | 12 供体 × 5 基因 Treg 伪 bulk 计数/总计数/年龄 |
| `tables/ST139b_treg_detection_by_donor.csv` | 供体级检测细胞数/检测率 |
| `tables/ST140_treg_age_analysis.csv` | **主分析表**：等权 WLS（主）+ Spearman + 合并 Poisson（敏感性）+ 族内 BH-FDR |
| `tables/ST141_cell_source_expression.csv` / `ST141_cell_source_expression_rate.csv` | Immune 图谱 12 基因 × 45 细胞类型的来源定位（率与检测率两口径） |
| `tables/ST142_*_pseudobulk.csv` | IL18/TNFSF14 主导区室供体伪 bulk |
| `tables/ST143_translational_age_glm_results.csv` | 转化族年龄斜率（多为不可估，如实报告） |
| `tables/ST144_treg_age_sensitivity.csv` | LOO / 大供体子集 / 环境代理证伪检验 |
| `tables/ST145_lung_cell_source_expression.csv` | 肺图谱 IL18/TSLP/IL33/TNFSF14 上皮定位 |
| `figures/FigA_treg_il5_il4_vs_age.{png,pdf}` | 供体级 Treg IL5/IL4 表达率-年龄 |
| `figures/FigB_immune_cell_source.{png,pdf}` | Immune 图谱细胞来源（率口径） |
| `figures/FigC_lung_localization.{png,pdf}` | 肺上皮定位（IL18/TSLP/IL33） |
| `scripts/00–08*.py` | 全部提取与分析脚本（可复跑） |

## 复跑

```
E:\ldsc_venv\Scripts\python.exe scripts\01_extract_obs.py
... 02_extract_counts.py   # ~3 min, 流式 CSR 提取
... 03_donor_analysis.py   # ST139/139b/141/142/143
... 04_treg_sensitivity.py # ST144
... 05_lung_extract.py     # 肺图谱
... 06/08 图
... 07_rebuild_st140.py    # 主表
```
