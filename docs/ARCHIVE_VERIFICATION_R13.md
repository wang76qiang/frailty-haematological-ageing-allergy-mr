# R13 修订执行与归档终检记录（2026-09-23）

**依据**：`编辑部审稿意见_第七轮_R12修订稿复审_2026-09-23.md`（QC1、QC2、文本1、程序1、建议项×2）
**性质**：全部文本/表注级改动，未改动任何分析结果或数值。

---

## 一、改动清单

| 项 | 位置 | 改动 |
|---|---|---|
| QC1（必做） | `01_手稿_manuscript/manuscript_main.md` Results（反向 MR 节，嗜酸分解句） | 括号内加入 F 统计定义：per-SNP 一阶 F = (β/SE)² 计算于各 SNP 被选择的 Pan-UKBB 成分 GWAS；mean F = 1,229 由全基因组显著选择所致（含 winner's-curse 上偏说明），仅证明选择强度，不与 frailty 工具 F ≈ 46（较小源 GWAS）直接可比 |
| QC1（必做） | `03_补充信息_SI/supplementary_information.md` ST133 表注 + Note 10(v) | 同步加入上述 F 约定，明确两口径不可比 |
| QC2（必做） | `01_手稿_manuscript/manuscript_main.md` Discussion（机制段） | 新增解释遗传层（ST133）与表型层（ST124）嗜酸剔除分歧的三个非排他候选机制：①表型重建在固定 100% 白细胞比例内重新分配权重（中性/淋巴百分比机械吸收嗜酸份额），而遗传层直接移除嗜酸载体 SNP；②嗜酸 SNP 标记邻近 2 型调控位点，移除即切除位点级信号；③遗传层预测登记制 GWAS 终点 vs 表型层预测自报过敏/IgE。结论：两层并读为界定嗜酸贡献的边界，而非相互矛盾 |
| 文本1（必做） | `01_手稿_manuscript/manuscript_main.md` 局限第二条 | 删除 "(the independent-exposure analysis is a registered pending analysis)" 承诺句 |
| 文本2（建议） | `01_手稿_manuscript/manuscript_main.md` Discussion | "illustrates a general principle for post-GWAS biology" → "is consistent with a broader caution in post-GWAS biology" |
| 建议项1 | 手稿 Methods（OneK1K 段）+ Results（OneK1K 段）+ SI Note 5(iii) + ST88 表注 | 逐工具 F 标注：39 工具 F = 5.1–23.2；14 个 FDR 显著检验的 F 跨度 9.5–23.2；两个低于 F > 10 惯例的检验点名标注（TNF rs2229094，CD8T，F = 9.5；IL6 rs61480128，CD8T，F = 10.0），数值取自 `ST88_r7_onek1k_celltype_mr.csv` 的 fstat 列 |

## 二、一致性核查

- SI 中已无 "pending analysis" 或 "general principle" 残留（全文检索确认；仅审稿意见档案保留原文记录）。
- ST88 表注新增声明与 ST88 CSV fstat 列逐值一致。
- ST133 F 定义与管线源码 `src/r16/01_eosinophil_free_mr.py`（`mean_F = mean((bx/se_x)²)`；bx/se_x 来自 `r3_index_instruments.csv` 的 beta_I/se_I，即成分 GWAS 中的 SNP 效应）核对一致。

---

## 三、归档终检（程序1）：✅ 本地包完整

| 核验项 | 结果 |
|---|---|
| `04_补充数据_supplementary_data/` 盘上文件数 | 131 = ST01–ST117（117）+ ST116b + ST124–ST135（12）+ ST135b |
| 期望清单比对 | Missing: 无；Unexpected: 无 |
| `supplementary_data_index.csv` 索引行数 | 131，与盘上文件一一对应 |
| `supplementary_tables_csv.zip` 内容 | 131 个条目，与期望清单一致，无缺失 |
| 手稿/投稿信声明文本 | "Additional files 1–117 and 124–135 (ST116b, ST135b companions)" 与本地结构一致 ✅ |
| Zenodo concept DOI | 手稿与投稿信均为 `10.5281/zenodo.22152039`，两处一致 ✅ |
| STROBE-MR 清单 | 已含 ST133/ST134/ST135/ST135b 行 ✅ |

## 四、仍需作者网页操作（投稿前必做，代理无法代办）

- [ ] GitHub 推送 R10–R17 全部新增文件（线上冻结点仍为 ST1–106 时代内容）
- [ ] Zenodo 发布新版本（versioned DOI 不变）；完成前 `Data availability` 声明不为真
- [ ] OSF osf.io/f4v79 / 注册 qm2hy 页面确认公开可见
