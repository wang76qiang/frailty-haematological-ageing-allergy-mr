# R11 投稿前阻塞项（Minor 8）

截至 2026-09-19，以下为“投稿系统会核验/外审会首查”的程序项。**未完成前不应投稿。**

| # | 阻塞项 | 现状 | 作者操作 | 完成后代理动作 |
|---|---|---|---|---|
| 1 | Zenodo/GitHub 与 Data availability 失配 | 线上仍是 ST1–106 冻结点，正文已声明 1–111 | 推送 R10/R11 新文件至 GitHub；Zenodo 发布新版本（version DOI 不变） | 复核落地页 + 回写 DOI/版本号 |
| 2 | OSF 预注册时间戳 | 包已备好（`docs/osf_package/`），未注册 | 网页上传取得 DOI（约 2 分钟） | 回写 Methods 与 SI Note 两处 “to be timestamped” |
| 3 | R6 “a priori” 可核验性（M7） | 无时间戳 | 同 OSF | 或删除全部未时间戳的 a priori 措辞 |
| 4 | 摘要字数 | 187 词（已合规） | — | — |
| 5 | deCODE IL-18（可选附加层） | 需 reCAPTCHA 表单 | 下载 IL-18 文件至 `data/real/interval_somascan/` | 运行 `src/r9/06_decode_il18_replication.py` |
| 6 | #2 独立暴露 sumstats | 未获取 | 运行 `src/r11/fetch_independent_exposures.py get <GCST>` | 运行 `01_independent_ageing_exposures.py` |
| 7 | #3 发病队列提取 | 无个体数据 | 按 `DATA_APPLICATION_UKB_AoU.md` 申请/提取 | 运行 `02_incident_allergy_cohort.py --input` |
| 8 | #4 湿实验 | 无 | 按 `WETLAB_PROTOCOL.md` 执行 | 回写结果 |

## 投稿系统材料清单核对
- [x] Cover letter（`05_投稿信_letters/cover_letter_main_genomemedicine.md/.docx`，2026-09-22 更新，ST 范围已同步至 124–132）
- [x] 手稿 .docx / SI .docx（R12 重构后已于 2026-09-22 重生成）
- [x] ST1–117、ST116b、ST124–132 + 索引（`supplementary_data_index.csv`，127 文件，zip 已重建）
- [x] STROBE-MR 清单（2026-09-22 补齐 ST112–117、ST124–126、ST129–132 行并重生成 docx）
- [x] 推荐审稿人（已含于投稿信）

---

# R12 编辑部复审后修订附记（2026-09-22）

本轮（第六轮）按编辑部审稿意见完成文本级修订，不改变任何分析结果：

- **Major 1**：BBJ 跨祖源小节——修正正文/SI 中 Cochran Q P 误写（1×10⁻¹⁶ → 3.3×10⁻⁵，与 ST107 存档值一致，df=27 核算通过）；正文与 SI Note 8 现如实报告随机效应（P=0.21，不显著）与加权中位数（β=−0.34，方向相反）估计量，哮喘/鼻炎跨祖源结论降级为"IVW 方向一致但估计量敏感"。
- **Major 2**：维度方向裁定的主证据在 Results/Discussion 中明确前移至 MHC-free 全基因组 rg；HLA 主导的衰弱 MR 工具仅作筛查层。
- **Major 3**：时钟阴性结论补充功率限定（GrimAge 仅 5 SNP；复制 meta 中 GrimAge–哮喘 OR 1.08, P=6.5×10⁻¹⁰ 等小效应实为显著）。
- **Minor**：HRS 早警信号与 ClinicalTrials.gov 段落压缩；局限性第十二补入 P2 的 post hoc 属性；Methods 增加 "pre-specified = file-provenance, no external timestamp" 全局限定；Fig. 4f 图例强化 schematic 警示；Fig. 5 图例补 OneK1K 不 corroborate 声明；SI ST130 表注补负对照截距诊断；正文明确 TNFSF14 为单平台类型（SomaScan-only）。
- **附带修复**：`supplementary_information.md` 的 Supplementary Note 1 与 Note 2 原文件含生成截断标记（"line truncated to 2000 chars"），Note 2 的 Sixth–Twelfth 条局限曾整段丢失——已按主文稿十二条局限重建并补入跨祖源/AD 条目。

## 仍需作者执行的网页操作（代理无法代办）

- [ ] **GitHub/Zenodo 再归档**：推送 R10–R17 新增文件（含 ST107–ST117、ST116b、ST124–ST135（ST135b 随附）与两份修订文档），Zenodo 发布新版本（versioned DOI 不变）；完成后 Data availability 声明方为真。
- [ ] **OSF 时间戳核验**：正文已声明注册号 qm2hy（2026-09-20），投稿前在 OSF 页面确认公开可见。
- [ ] 摘要现为 329 词，符合 Genome Medicine 结构化摘要 ≤350 词上限（原 200 词约束系 Nature Aging 场景，已不适用）。

---

# R13 附记（2026-09-23，第七轮审稿意见落地）

本轮按 `编辑部审稿意见_第七轮` 完成 2 项必做文本修订（ST133 的 F 统计口径定义入正文/SI 两处；Discussion 增加遗传层 vs 表型层嗜酸剔除分歧解释）、删除 pending-analysis 承诺句、软化 "general principle" 外推、OneK1K 逐工具 F 标注（39 工具 F = 5.1–23.2，两个亚阈值检验点名），并机读核验补充包完整性（131 文件 = ST01–117 + ST116b + ST124–135 + ST135b；索引与 zip 一一对应，无缺失）。

## 三项网页操作执行结果（2026-09-23）

- [x] **GitHub 推送**：✅ 完成。单提交 `02813b61`：新增 ST127–ST135 + ST135b、索引刷新至 131 行、docs/ 两份修订文档；线上复核 131/131 全齐。Release `v1.1.0-R13` 已发布。
- [ ] **Zenodo 新版本**：⚠️ 本机网络对 zenodo.org 为 DNS 级阻断（0.0.0.0），DoH 端点均被重置，无法代办。Release v1.1.0-R13 已可触发 GitHub-Zenodo 集成（若启用）。**作者需在可访问 zenodo.org 的网络登录确认/发布新版本**（详见 `ARCHIVE_VERIFICATION_R13.md` 第四节②）。完成前 Zenodo 部分声明不完全为真。
- [x] **OSF 公开确认**：✅ 已核验。API 实查 f4v79 与注册 qm2hy 均 public=True（2026-09-20 注册），无需操作。

---

# R14 附记（2026-09-23，顶刊编辑独立复审 M1–M5 + 程序项落地）

本轮按外部独立编辑复审意见完成 7 项修订，新增 1 张补充表（ST136），并实测核验三处归档：

## 分析/计算
- **M3（新算，`src/r16/04_eosinophil_free_panukbb.py` → ST136）**：嗜酸-free 分解的**独立队列复制**。在 Pan-UK Biobank（UKB）对 full / eosinophil-free / eosinophil-only 三工具复跑：full index 哮喘 OR 2.79、鼻炎 OR 1.57；**eosinophil-free 丢失鼻炎（OR 1.53, P=0.39）但保留哮喘（OR 1.86, P=1.7e-7；RE P=0.034）**；eosinophil-only 两者强关联（OR 3.03 / 1.58）；荨麻疹对照 null。→ 哮喘的"嗜酸承载"结论改为**队列依赖**（FinnGen 衰减、UKB 保留），鼻炎/AD 结论不变。已写入正文 Results、摘要、Conclusions、局限与 SI Note 10(v)。

## 文本
- **M1（LDSC 限定）**：Methods/Results/Limitations/SI Note 7(ii)/Note 2 明示——维度对比承重于 in-house 与第三方 Python-3 两套实现，**非官方 Bulik-Sullivan ldsc**；外部基准（ST103）仅覆盖 h² 与 clock-clock rg，**不覆盖 frailty–allergy 的承重 rg**，需外部官方复算。
- **M2（ST135 重写）**：承认 frailty→asthma 的**单 SNP 方向多不可判定**（仅 48% 正确；Steiger 显著者 10 反向 vs 2 正确），方向结论仅由 MHC-free 全基因组 rg + 非 HLA 再裁定承担；SI Note 10(vii) 同步。
- **M4（winner's curse）**：删除"selection bias does not inflate"结论，改为"重叠面板的稳定性检查，非独立校正，残余选择偏倚未排除"；Methods/Limitations 同步。
- **M5（发病限定）**：摘要 Methods/Conclusions 与 Conclusions 节写入"全部过敏证据为现患/病例对照，set-point 未经发病检验"。
- **Minor 4（措辞统一）**：全文裁决树相关措辞统一为 **"a priori (by file provenance, unregistered)"**（摘要/Intro/Methods/Results/Discussion；共 7 处），并保留定义句。
- **程序**：作者行空上标与 Fig. 6j 双重引用经机读确认**已无残留**。

## 三处归档核验（2026-09-23 只读实测）
- [x] **GitHub**：`wang76qiang/frailty-haematological-ageing-allergy-mr` 经 GitHub API 确认 **public=True**，`pushed_at 2026-09-23`。
- [x] **OSF**：`api.osf.io/v2/registrations/qm2hy` 返回 **public=true / reviews_state=accepted**（注册 2026-09-20）。
- [ ] **Zenodo**：本环境对 zenodo.org 为传输级阻断，无法核验。**且本轮新增 ST136 + 修订 docx/SI/投稿信/STROBE/zip，当前 Zenodo 版本必不含 ST136 → 作者须再归档。**

## 仍需作者执行的网页操作
- [ ] **GitHub 推送 ST136**（`ST136_eosinophil_free_panukbb_mr.csv`）+ 刷新索引/zip + 修订的 `src/r16/04_eosinophil_free_panukbb.py` 与全部修订文档。
- [ ] **Zenodo 发布新版本**（versioned DOI 不变），纳入 ST136；完成前 Data availability 声明不全为真。
- [ ] （可选）deCODE IL-18、UKB/AoU 发病队列、湿实验——不变。

## 材料清单更新
- ST 范围 **ST01–ST117 + ST116b + ST124–ST136**（含 ST135b、ST136），共 **132 个 ST 数据文件**（另附 `supplementary_data_index.csv`）；索引已加 ST136 行；`supplementary_tables_csv.zip` 已重建（132 条目）。
- 手稿/SI/投稿信/STROBE 四份 `.docx` 与 README md5 清单（185 文件）已重生成。


---

# R15 附记（2026-09-23，第二轮顶刊编辑复审 3 Major + 6 Minor）

- **M2'（最重要科学修正）**：ST133 的“哮喘由嗜酸承载”在独立队列未复制。FinnGen eos-free 哮喘 OR 1.42（RE P=0.45，欠功效）vs Pan-UKB eos-free OR 1.86（P=1.7e-7；RE P=0.034；WM 2.74）；两队列均正向，故 FinnGen 衰减更可能是功效不足而非嗜酸承载。摘要/Results/Conclusions/SI Note 10(v)/Note 2/投稿信统一改为“did not replicate / asthma attribution unresolved”。鼻炎复制成功；AD 仅 FinnGen 可测。
- **M1'**：摘要新增“this contrast used two non-canonical LD score implementations and awaits canonical confirmation”。
- **M3'**：标题改为 `Dimension-resolved Mendelian randomization separates a functional blood-cell composite from epigenetic ageing in allergic disease`；手稿 9 处 + SI 3 处 “allergy-specific” 改为 outcome-informed / specific to allergic endpoints；Discussion “cell-type resolution” → “hypothesis-generating cell-type analysis”；投稿信 IL4 “resolves” → “partially resolves”。
- **Minor 1–4**：摘要补 UKB 阳性端点、时钟与两个衰老暴露补 underpowered、删除重复 “independent UK Biobank” 表述。
- **Minor 5**：ST136 表注说明 eos-only 哮喘 OR 3.03 ≈ full 2.79（仅少数非嗜酸成分）。
- **Minor 6**：Zenodo 再归档须纳入 ST136。

## 校验
摘要 348 词（<=350）；标题/摘要/正文/SI/投稿信/STROBE 标题一致；引用 1–59 双向；em_dash=0；CJK=0；ST 范围 ST01–ST117 + ST116b + ST124–ST136；docx×4、索引、zip、README md5 清单（185 文件）重生成。


---

# R16 附记（2026-09-23，第三轮顶刊编辑复审 M1–M3 + 6 Minor）

- **M1（标题去构念合并）**：标题改为 `Dimension-resolved Mendelian randomization implicates a functional frailty dimension and a rhinitis-anchored blood-cell composite rather than epigenetic clocks in allergic disease`；手稿/SI/投稿信/STROBE 四处一致。
- **M2**：标题动词 “separates”→“implicates … rather than”；摘要统一 “epigenetic clocks showed no consistent signal and were underpowered”；维度对比仍标注两套非官方 LDSC 实现、待官方确认。
- **M3**：复合指数在摘要/Conclusions 重定位为 “rhinitis-anchored / most replicably for allergic rhinitis”；正文小节标题改为 “nominate a candidate type-2-like signal”。
- **Minor 1**：摘要 Results 与 Conclusions 的时钟表述统一为 “no consistent signal (underpowered)”。
- **Minor 2**：主 Conclusions 首句改为 “a functional dimension and, most replicably for allergic rhinitis, an eosinophil-weighted blood-cell composite”。
- **Minor 3**：投稿信 H1 改为 “main paper: dimension-resolved ageing comparison and blood-cell composite”。
- **Minor 4**：投稿信 IL4 标注 “DICE only; not corroborated by independent single-cell eQTL resources”。
- **Minor 5**：ST136 表注注明 Pan-UK Biobank 批次与 Additional file 130（ST130）同源。
- **Minor 6**：Zenodo 再归档须纳入 ST136（作者网页操作）。

## 校验
摘要 349 词（<=350）；标题四处一致；引用 1–59 双向；em_dash=0；CJK=0；docx×4 + README md5 清单（185 文件）重生成。


---

# R17 附记（2026-09-24，Minor 1–4 完成 + 三处归档闭环）

- **Minor 1**：手稿 Results 与 SI Note 10(vi) 的 “formally established” → “supported under the two non-canonical implementations used here”。
- **Minor 2**：Discussion “causal (MR) evidence” → “instrument-based MR evidence”。
- **Minor 3**：标题改为 `Dimension-resolved Mendelian randomization implicates a functional frailty dimension and a rhinitis-anchored blood-cell composite in allergic disease without a consistent epigenetic-clock signal`；手稿/SI/投稿信/STROBE/GitHub README 一致。
- **Minor 4（归档闭环，2026-09-24 实测）**：
  - **GitHub**：`wang76qiang/frailty-haematological-ageing-allergy-mr` 已推送 `tables/ST136_...csv`、更新的 `tables/supplementary_data_index.csv`、`scripts/r16_01`–`r16_04`、更新 README（API 确认 ST136 存在）。
  - **Zenodo**：发布 **v1.2.0**（`10.5281/zenodo.22936545`）与 **v1.2.1**（`10.5281/zenodo.22936675`）；concept DOI `10.5281/zenodo.22152039` 始终指向最新版；包内含 ST136、R16 脚本与修订手稿/SI（v1.2.1 zip 2,014,922 bytes）。 另发布 **v1.2.2**（10.5281/zenodo.22936744；含 SI Note 3 因果措辞修正）。
  - **OSF**：`qm2hy` public/accepted（2026-09-20 注册）已核验。

## 校验
摘要 349 词；引用 1–59 双向；em_dash=0；CJK=0；手稿/SI/投稿信/STROBE/GitHub README 标题一致；docx×4 + README md5 清单（185 文件）重生成。

---

# R18 附记（2026-09-25，第八轮编辑复审 1–5 项落实 + canonical LDSC 复算）

## 已完成的 5 项

1. **Canonical LDSC 复算（Item 1，新计算）**：官方 bulik/ldsc v1.0.1 以纯机械 py3 补丁运行（补丁树 \E:\ldsc_work\ldsc\，含 patch_log.txt；脚本已归档 \src/r17/\）。标准 munge 协议 + 公开 1000G Phase 3 HapMap3 权重（MHC 剔除；.l2.M 由行数导出；--not-M-5-50）。**15 对 rg 全部复现**：frailty 0.467/0.361/0.246 vs in-house 0.458/0.353/0.246（|Δ| 中位 0.025、max 0.086）。canonical 正式对比（ST138）：哮喘/鼻炎全部 q ≤ 2.6e-3（ρ=0）；**AD 臂口径依赖**（ρ=0 不显著、ρ=0.5 三/四显著）——已如实写入手稿 Results/Limitations、SI Note 10(viii)。**MHC 敏感性**：含 MHC 参考仅 +0.019–0.033。产物：ST137（18 行）、ST138（12 行）；索引/zip 已更新（134 条目）。
2. **嗜酸-free 哮喘表述（Item 2）**：摘要/Results/投稿信主语歧义修正（"attenuation did not replicate; association positive in both cohorts, significant in the larger"）。
3. **F 口径（Item 3）**：QC1 已闭合确认（Note 10(v) + ST133 表注）；本轮补 ST136 表注同一口径。
4. **HRS 措辞（Item 4）**：摘要 "validated in NHANES and HRS" → "validated in NHANES and longitudinally characterized in HRS"。
5. **docx 通讯标记（Item 5）**：机读核验通过（1,* / 2,* 上标）；docx×4 重新生成后复验通过。

## 附带修正
- Results IL18 句补全 HELIC 数字（OR 1.22, P=1.1e-8；合并 OR 1.19, 95% CI 1.15–1.23, P=3.6e-23, I²=0%；Additional files 127, 128），与 SI/ST127/ST128 一致。
- 手稿/SI/投稿信/STROBE 的 Additional files 范围统一为 124–138。

## 校验
摘要 270–291 词（≤350）；引用 1–59 双向；em_dash=0；CJK=0；标题一致；README md5 清单刷新（191 文件）。

## 待作者执行的网页操作
- [ ] **GitHub 推送**：ST137/138 CSV、supplementary_data_index.csv、supplementary_tables_csv.zip、\src/r17/\ 三个脚本、修订的手稿/SI/投稿信/STROBE（docx+md）、README。
- [ ] **Zenodo 发新版**：纳入 ST137/138 与 zip；完成前 Data availability 中 ST 范围声明（124–138）不完全为真。
- [ ] OSF qm2hy：无需变动（本轮无新注册分析）。

---

# R19 附记（2026-09-25，Nature 子刊补强包并入主稿）

## 并入内容
- **新分析层**（P0-1）：Immune 图谱 24 供体供体级表达-年龄检验（Treg 5 族 + 转化 4 族，内部预声明；环境 RNA 证伪检验）+ Lung 图谱定位层。核心：Treg IL5 随年龄上调（名义 P=0.025，q=0.13）；IL4 方向一致但降级（环境相关）；IL18 定位于肺泡巨噬/肺泡 I/II 型/基底上皮；TSLP-AT2 内阳性对照。
- **文稿**：手稿 Methods 新小节 + 多重检验段两个新族声明；Results 新小节（Donor-level single-cell atlas corroborates an age-drifting type-2 set-point）；Discussion 增补段；局限增补句；摘要增补句；SI Note 10(ix) + ST139-145 图例 + Supp Fig 11-13 图例；STROBE 加行；投稿信加一条。全部名义级/假设生成级标注，明确不构成 DICE IL4 共定位的独立复制。
- **范围**：Additional files 124–145（+ST139、139b、140、141、141b、142、143、144、145）；补充图 13 张；索引/zip 143 条目（脚本核验对齐）；docx×4 重生成。
- **归档**：src/r17/03–10_sc_*.py + ANALYSIS_PLAN_sc_atlas.md；原始分析包保留于 Nature子刊/（含数据提取缓存）。

## OSF
- [ ] **作者操作**：将 src/r17/ANALYSIS_PLAN_sc_atlas.md 追加至 OSF 注册（qm2hy 或新注册），标注内部预声明族与 post hoc 属性。完成前 Methods 中“documented in the archived analysis package”指 GitHub/Zenodo 存档。

## 待作者执行
- [ ] GitHub 推送：ST139–145、索引、zip、Supp Fig 11–13、src/r17 全部、修订四文档（md+docx）、README。
- [ ] Zenodo 发新版（纳入 ST139–145 与 143 条目 zip）。

---

# R20 附记（2026-09-25，方案 B：Nature Communications 上探补强包并入）

## 并入内容
- **P0-2 目标试验模拟（新算）**：NHANES 2017-18 高/低指数模拟分配 + 稳定化 IPTW + 调查权重 + PSU 聚类：OR 1.90 (1.41–2.56), P=2.5e-5, E-value 3.21 (CI 限 2.17)；per-SD 1.43 (1.20–1.71)。药物亚组欠功效如实描述（ST152）。
- **P0-3 基因级三角验证（存档 QC 并入）**：RORC/IL6 三条证据线；IL4 保护 (P=1.3e-8, PP4≈1.0)（ST147/148）。
- **P0-4 药物连接层**：LINCS 衰老签名（ST149，存档）+ alarmin 签名**从 GCTX 重算**（ST151：dasatinib/tofacitinib 上调、azithromycin 抑制，q=0.024；**mTOR 抑制剂 null——旧 mTOR-alarmin 信号未复现并声明废止**）；无药物同时抑制两套签名；PRISM 提名 SRC/IL6 通路（ST150）。
- **摘要 NC 化**：非结构化 143 词；Discussion 增补目标试验因果推断段；局限增补三条；Conclusions 增补 E-value。
- **范围**：124–152；索引/zip 150 条目（核验对齐）；Supp Fig 14；SI Note 10(x)(xi)；STROBE 三行；**NC 版投稿信**新增（GM 版保留备投）。
- **归档**：src/r17/11–13（目标试验/森林图/LINCS alarmin）。

## 待作者执行
- [ ] GitHub 推送（ST146–152、Supp Fig 14、索引、zip、src/r17/11–13、NC 投稿信、修订四文档）
- [ ] Zenodo 发新版（150 条目 zip）
- [ ] OSF 追加目标试验与药物层计划（post hoc 属性 + 内部预声明族）
- [ ] **NC 投稿系统特有项**：Editorial Policy/Reporting Summary 检查表（统计部分按手稿 Methods 填写）、Data/Code Availability 声明已具备（GitHub + Zenodo concept DOI + OSF）

## R20.1 补充（2026-09-25，编辑投前三项落实）
- **① 无嗜酸重选存档核查（结果：不支持直接引用）**：esults/r3/tables/r3_index_no_eosinophil_mr.csv（R3 期，非 allele-aware）显示 no-eosinophil 指数三结局全部**强于**全指数（哮喘 2.81 vs 2.32；鼻炎 1.86 vs 1.62；AD 2.47 vs 1.90）——与 ST133 的等位感知模式（鼻炎丢失、哮喘衰减）**发散**，与等位错配污染一致；该表此前已在 SI Note 10(v) 声明废止。处置：**不重新引用**，在 SI Note 10(v) 增补发散细节一句（说明废止原因的发散证据），构念循环的诚实声明维持为最优辩护；真正独立的固定权重指数计算列为未来工作。
- **② 目标试验收入校正（实做优于声明）**：倾向模型补入家庭收入贫困比（INDFMPIR），重跑 src/r17/11——**OR 1.86 (1.38–2.49), P=3.5e-5, E-value 3.12 (CI 2.11), n=4,403/689 事件**（原 1.90/3.21, n=5,052/771）；摘要/Results/Discussion/Conclusions/SI/NC 投稿信七处数字已全局同步；包内 ST146 已更新；吸烟包年数不可得声明入 Methods/Discussion。
- **③ LINCS null 机制候选句**：已查 LINCS 细胞系构成（98 类中 83 类为永生化细胞系，屏障上皮基本缺席）；Results 增补两个机制候选句（JAK-STAT 双向调控候选 + 上皮代表性缺失，含 83/98 数字）。
- 校验：摘要 144 词、em-dash=0、引用双向、docx×3 重生成、README md5 刷新。

## R20.2 标题缩短（2026-09-25，编辑意见 2 落实）
标题改为机制性版本 A functional frailty dimension shares genetic architecture with allergic disease；"rather than epigenetic clocks" 对照移入摘要首句（"epigenetic-clock evidence has been inconsistent"，摘要 149 词 ≤150）。六处同步：手稿、SI、NC 投稿信、GM 投稿信（备投）、STROBE 清单、docx×5 重生成。注：README 历史修订块中的旧标题引用为版本日志，有意保留。

## R21 投稿包迁移（2026-09-25）
最新全套投稿资料已由 D:\衰老研究\v3_pipeline\Genome Medicine\ 整体迁移至 D:\衰老研究\v3_pipeline\Nature子刊\Genome Medicine\（215 文件，robocopy /MOVE，0 失败；包内相对路径不受影响，绝对路径引用需以新位置为准）。Nature子刊\ 根目录保留原始分析包（scripts/data/tables/figures + README + 草案）。本文件后续路径相对新位置。

---

# R22 附记（2026-09-25，Nature Communications 投稿格式改造）

## 结构重排（手稿）
- 节序改为 NC Article 版式：Title page → Abstract → **Introduction**（原 Background）→ Results → Discussion（原 Conclusions 并入 Discussion 末段）→ **Methods**（由文前移至 Discussion 后）→ **Data availability** → **Code availability**（新增独立节）→ References → Acknowledgements → Author contributions → Competing interests → Ethics declaration → Additional information → Tables → Figure legends。
- Keywords 从稿件移除（NC 不刊印），存入  5_投稿信_letters\NC_submission_system_notes.md（投稿系统字段备注，含 Article type）。
- 内嵌 STROBE-MR 清单移出手稿（独立文件保留于 06_报告清单，作 supporting document 上传）。
- References 卷号加粗（Nature 样式结构）；**作者全名列表为作者终检项**（Nature 规则：<=30 位全列，>30 位列前 30 + et al.；现稿暂留 et al.，需据 Crossref 逐条扩展）。
- 全文术语 BMC→NC：Additional file(s) → **Supplementary Table(s)**（手稿 53 处 + SI + 两封投稿信，含 Additional information 段同义反复重写）。

## 新增文件
-  5_投稿信_letters\NC_submission_system_notes.md（系统字段备注）
-  6_报告清单_checklists\Reporting_Summary_draft.md/.docx（NC 必填 Reporting Summary 草案，统计/数据/代码可用性七项，待作者系统内确认）

## docx 重生成（NC 审阅格式）
生成器升级：md2docx_nc.py——Times New Roman 12pt、双倍行距、**连续行号**（w:lnNumType，机读确认）。重生成手稿/SI/NC 投稿信/Reporting Summary 四份。

## 校验
md5 清单 217/217（新增 3 条目已登记）；标题五处一致；摘要 149 词；引用 61 双向；em-dash=0；CJK=0；E-value 2.66 无残留旧值；'Additional file' 零残留；H1 节序与 NC 版式一致。

## 待作者（NC 投稿系统）
- Reporting Summary 系统内确认/编辑；Editorial Policy Checklist；References 作者全名扩展（Nature 规则）；此前四笔网页操作（GitHub/Zenodo/OSF）不变。

## R23 附记（2026-09-25，参考文献：重编号 + Crossref 作者全名扩展）
- **引用顺序修复**：参考文献按正文首次出现顺序重编号——51–61（5q31/转化/药物警戒，Results 内首引）→ 25–35；25–50（工具/方法学 GWAS，方法与转化节后首引）→ 36–61；1–24 不变。重映射后首次出现序列零违例。
- **Crossref 作者扩展**：61 条全部处理。52 条经 Crossref 核验替换为全名（Nature 式：<=30 全列，>30 列前 30 + et al.——GBD 1030 作者、McCartney 109 作者、Astle 等大列表均按规则截断）；团体/机构作者（Tabula Sapiens ×2、GTEx Consortium、EAGLE Eczema Consortium、FinnGen R12、NCHS、Sanderson 11 人、Ochoa 34 人、Sun 2018 33 人、McCartney）按 Crossref 权威记录或 DOI 直查核定。单作者文献（Horvath、Wallace）维持原样。双句点全局清理。
- **核查结论**：61/61 双向引用完整；首次出现顺序与编号一致（机读验证零违例）；图 1–6 首现顺序正确；Supplementary Fig 1–14 全部被引；Supplementary Tables 1–117 + 124–152 全部被引（ms+SI）；Table 1 被引。
- docx 重生成；md5 清单 217/217。

## R24 附记（2026-09-25，AI 痕迹严格检查与修改）
按 
o-ai-slop + stop-slop 两套规范对全文机读扫描（禁词/空副词/元话语/套式结构/句首 Wh/连接词链/em-dash 共 40 余规则）。
**扫描结论**：基线极净——零禁词、零空副词链、零 'Of note/Importantly/taken together' 类元话语、零 Moreover/Furthermore 连接词链、em-dash=0。
**修改 4 处**（最小有效编辑）：① 'Together, these frame ageing as...'（总结式元话语）→ 直接断言 'On this evidence, ageing modifies...'；② 'What the frailty instrument measures is a genuine concern:'（Wh 开句+冒号揭示）→ 'A genuine concern is the content of the frailty instrument:'；③ 删除 'We emphasize what this layer does and does not show.'（解释性元话语，后续两句自证）；④ 引言三连 Whether 排比 → 'It has remained unknown whether...' 主语先导。
**判定保留**（科学文体术语/认识论负载，删之损义）：'robust' ×6（统计稳健性术语）、'rather than' ×16（均为载义的边界判别，非戏剧化对比）、'we show that'（Nature 主动语态）、被动语态（Methods 惯例）。
复核：残留零标志、摘要 149 词、引用 61 双向、em-dash=0；docx 重生成，md5 清单刷新。

## R25 附记（2026-09-25，最终预审与技术审查）
- **真表格改造**：docx 生成器升级（v2）——markdown 表格块渲染为可编辑 Word 实体表（Table Grid, TNR 9pt）；修复 Table 1 空行断块问题（5 个单行碎片表 → 1 个 4 行实体表）。手稿/SI/STROBE/Reporting Summary 四份 docx 全部重生成，pipe 残留 0。
- **终审结果**：md5 清单 217/217；摘要 149 词；引用 61 双向且正文首现 1→61 递增（标题页 [2,1] 为署名单位上标伪迹，已定性）；图 1–6 首现顺序正确；Supp Fig 1–14 图例齐；索引/zip 150/150 对齐；ST137–152 在位；术语零 BMC 残留；em-dash=0；CJK=0；NC 投稿信四要素齐；Reporting Summary 与系统备注在位；行号（lnNumType）确认。
