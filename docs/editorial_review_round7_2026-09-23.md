# 编辑部第七轮审稿评估报告（R12 修订稿 · 独立复审 + 数据级实证核验）

**稿件**：Dimension-resolved Mendelian randomization identifies an eosinophil-dominated blood-cell (type-2-like) axis in allergic disease
**评估日期**：2026-09-23 | **前轮**：第六轮（Major Revision：M1–M5）→ R12 修订稿
**本轮性质**：第六轮意见的**实证核验复审**（verification re-review）——对 M1–M5 的关闭判定不采信文本声明，逐一回溯到底层数据文件核验；同时进行新一轮独立挑错
**决定建议**：**Minor Revision（接受轨道）**——第六轮全部 5 项 Major 经验证为真闭合，无新 Major。剩余 2 项 QC 澄清 + 4 项文本/程序尾项，预计 2 个工作日内完成。本稿现阶段的证据-主张匹配度已达到 Genome Medicine 外审安全线，建议按期投稿。

---

## 一、第六轮 Major 意见的实证核验（逐条回溯数据文件，非纸面核验）

| 意见 | 核验方式 | 裁定 |
|---|---|---|
| **M1**：rg 对比需正式差异检验 | 调取 `ST134_frailty_vs_clock_rg_contrast.csv` 逐行复核：哮喘 4 个对比、鼻炎 4 个对比，in-house（ST60）与独立实现（ST126）两套估计全部 ci_nonoverlap = True，两套 q 值均 ≤ 1.4×10⁻³，与摘要 "q ≤ 1.4 × 10⁻³" 精确一致；AD 臂 8 个对比全部不显著（含 PhenoAge rg 0.274 > frailty 0.250），与正文 "not atopic dermatitis" 一致。正文 Results 已引用（第 74 段 "A formal pairwise contrast confirmed…"）。ρ=0 与 ρ=0.5 两种 SE 口径并列，保守口径下哮喘/鼻炎仍全部 q < 3×10⁻⁴。 | ✅ **真闭合，且执行质量高于要求** |
| **M2**：遗传层嗜酸剔除分析 | 调取 `ST133_eosinophil_free_instrument_mr.csv`：嗜酸-free 工具保留 AD（OR 2.02，P = 2.0e-9，RE P = 5.9e-3）、哮喘衰减至名义（OR 1.42，P = 0.046，RE P = 0.45）、鼻炎 null（OR 0.80，P = 0.17）；嗜酸-only 工具三结局全强（OR 1.85–2.59）。正文（Results 反向 MR 节末段 + 局限第二条）如实报告 "predominantly eosinophil-carried… component-level"，标题 "eosinophil-dominated axis" 措辞与数据一致。SI Note 10(v) 对既有非 allele-aware 内部估计的废止声明齐备。 | ✅ **真闭合** |
| **M3**：8.9 倍对比非方向检验 | 正文已将 8.9 倍明确标注 "(a cross-scale magnitude comparison, reported descriptively)"（第 81 段）；Steiger 逐 SNP 审计（`ST135_steiger_directionality.csv` + ST135b）入正文：4 个钟 100% 方向正确、frailty 48–80%、错误方向 SNP 集中于免疫位点（70–100%）。 | ✅ **真闭合** |
| **M4**：IL4 小节标题断言化 | 小节标题已改为 "The IL4 paradox **partially resolves** into a whole-blood artefact and a candidate cell-type-specific signal"，与正文分级同步。Discussion "illustrates a general principle" 一句保留——**轻微残留**：建议软化（见二.2）。 | ✅ 基本闭合（1 处文本建议） |
| **M5**：摘要分级同步 | 摘要 Results 现逐句携带分级："marginally direction-consistent but estimator-sensitive in Biobank Japan"、"replicated (post hoc)"、"q ≤ 1.4 × 10⁻³… but not atopic dermatitis"。Conclusions 同步。 | ✅ **真闭合** |

---

## 二、本轮新发现问题（2 项 QC + 4 项文本/程序，无 Major）

### QC1（必答）：ST133 的 mean F = 1,229 与全文 F 统计体系不自洽

ST133 中 eosinophil_free / eosinophil_only / full_index 三个工具的 mean_F 分别为 1,229 / 955 / 1,182，而正文其余各处报告的工具 F 均在 10–130 区间（frailty F ≈ 46、嗜酸 count F = 85.5、BBJ 工具 F ≈ 9–23）。F > 1,000 对一个由弱效应血细胞 SNP 构成的工具在常规 F = β²/SE² 定义下不可解释。**外审统计审稿人必抓此项**。要求：在 ST133 表注或 Methods 明确该 F 的计算定义（疑似基于成分加权后的复合暴露近似，而非单 SNP 一阶 F），并与全文口径调和；若确为不同定义，改名（如 F_composite）避免与弱工具阈值 F > 10 的判据混读。此项不影响效应量本身，但必须在投稿前消除。

### QC2（建议）：遗传层 vs 表型层嗜酸剔除结果的分歧需要一段解释

ST133（遗传层）：嗜酸-free 工具丢失哮喘/鼻炎关联；ST124（表型层）：嗜酸-free 指数仍与自报过敏强相关（OR 1.63，P = 3.7e-17）。正文两处都如实报告了，但目前只并置、未解释。这恰是外审会问的第一个问题（为什么遗传剔除后哮喘关联消失而表型剔除后不消失——候选解释：NHANES 复合指数的权重在成分间重新分配、遗传工具剔除的是携带嗜酸等位效应的 SNP 而非全部嗜酸变异、或构念重叠集中于遗传层）。不要求补分析，要求 Discussion 增加 2–3 句对这一分歧的机理解释候选，并将其列入局限。

### 文本 1（必答）：删除或降级 "registered pending analysis" 表述

局限第二条末句 "the independent-exposure analysis is a registered pending analysis"——承诺一项未完成的注册分析是外审红灯（要么做、要么删）。建议删除该承诺句，或改为 "an independently weighted index is a natural next analysis"。

### 文本 2（建议）：Discussion "illustrates a general principle for post-GWAS biology"（第 109 段）

第六轮 M4 的残留。IL4 悖论解析的证据基础是 DICE 单资源 + significance-only eQTL（PP.H4 上偏已声明）+ OneK1K 不复制，不足以支撑 "general principle" 级别的外推。改为 "is consistent with a broader caution…" 或限定在 bulk-tissue colocalization 语境。

### 程序 1（投稿前必做）：归档与声明一致性终检

主文稿 Data availability 现声明 Zenodo doi:10.5281/zenodo.22152039（第五轮记录中曾出现 .22152133）。投稿前必须核验：① 该 Zenodo 记录与 GitHub 仓库实际包含 ST107–ST135、ST116b、ST135b 及 SI Note 10 全部产物；② OSF 注册 qm2hy 页面公开可访问；③ SI Note 10（ST124–ST135）本身已列入 Note 1 索引表（本轮核查：Note 1 目录含 Note 10，但 ST118–ST123 归属伴随论文的划分说明需与伴随论文同步）。三者缺一，Data availability 声明不为真。

### 程序 2（建议）：OneK1K 工具的 F ≈ 9–23 逐工具标注

第六轮 Minor 3 的残留。正文 Methods（第 52 段）只给出范围 "F ≈ 9–23"，而文稿自设 F > 10 门槛——部分工具（F = 9）低于门槛。在 ST88 表注或正文括号中逐工具列出 F 值，避免 "自定门槛自破" 的观感。

---

## 三、对四问的正式回答（本轮独立判断，不沿用前轮结论）

### 1. 现有研究是否完整？—— **约 95%**

计算证据体系**完整且超额完成**：五工具裁决树、六源复制 + meta、双实现 MHC-free LDSC（in-house + 独立 Python-3，ST126/134 互证）、免疫位点剔除 + 非 HLA 再裁定（ST131）、healthspan 外部锚（ST132）、NHANES 人群验证（含权重敏感性与表型层嗜酸剔除，ST124）、HRS 16 波纵向构念验证（含均值匹配早警信号）、5q31 多信号 coloc + 先验网格、pQTL 三平台复制（IL18）、安全性扫描、跨祖源边界（BBJ）、Steiger 逐 SNP 方向审计。每一层都配有一个本可推翻它的检验，且每一次检验失败（EAGLE、BBJ-AD、OneK1K-IL4）都被转化为边界而非隐藏——这一点在全部七轮评审中保持一致记录，是本稿最突出的方法论资产。

剩余缺口（均已诚实声明）：① **发病队列**——全部过敏证据为终生病例-对照/现患数据，set-point 模型本质是发病命题，用现患数据检验，这是本稿与"决定论级"结论之间最后的概念距离；② **功能验证**（湿实验）；③ 程序归档项。无未声明的计算缺口。

### 2. 现有计算结果是否足以支撑中心结论？—— **足以支撑文稿实际主张的三个层级，不足以支撑任何更强主张（而文稿未作更强主张）**

| 子主张 | 支撑度 | 本轮核验要点 |
|---|---|---|
| (a) 维度不对称：功能性维度 > 表观钟（哮喘、鼻炎） | **高**（M1 闭合后） | ST134 两套实现 16 个正式对比全部 q ≤ 1.4e-3，且 CI 不重叠；AD 臂不显著已被正确排除出主张。功效不对称的替代解释被 scale-free rg 层面的正式检验排除 |
| (b) 嗜酸主导血细胞轴与过敏关联 | **高（成分级限定后）** | 多层多队列一致；M2 闭合后构念循环在遗传层与表型层均被直接检验——哮喘/鼻炎为嗜酸承载（component-level）、AD 含非嗜酸成分，措辞已同步降级 |
| (c) 因果/可干预性 | **未主张——正确** | 正文 "not shown to be caused by"、反向链明确声明未直接检验衰老端点；Steiger 审计与反向 MR 仅界定方向 |

### 3. 结论是否可靠？—— **主结论可靠性高，分级可靠性格局清晰**

- 方向性（指数 ↔ 过敏）：**高**。五估计量、多队列、权重敏感性（+8.5% < 30% 阈值）、负对照、表型与遗传两层嗜酸剔除交叉验证。
- 维度特异性：**哮喘/鼻炎高，AD 中（EUR 内中-高、跨祖源不成立且已收缩）**。
- 量级：**已正确放弃**（I² = 90–99%、工具尺度 OR 477/6.6×10⁵ 自证不可解读）。
- BBJ 哮喘/鼻炎：正文分级 "marginal and estimator-sensitive" 与 ST107 数据一致（IVW q ≈ 0.04 但 RE P = 0.21 / 加权中位数方向相反或 null）——这是**诚实分级**，但请注意：外审可能认为该层不足以计入"跨祖源一致性"论据；建议 Conclusions 保持现措辞（"directionally consistent… marginal and estimator-sensitive"），不要在外审回应中升级为 replication。
- 机制（Treg-IL4）：**提名级**，DICE-only + significance-only eQTL 上偏 + 先验敏感性下 AD 端 PP.H4 在 p12 = 1e-6 时降至 0.54——分级正确。
- 转化（IL18 > TNFSF14）：IL18 为全文最硬靶点（三平台 + coloc + 先验网格 + 安全扫描方向一致）；TNFSF14 单工具主导（rs413141 占 71% 权重）已如实降级。

### 4. 证据链是否完整？—— **闭合。第六轮的两处结构性接缝已用计算补救，本轮无新接缝**

```
衰老(5工具) ──裁决树──> 维度特异性(frailty>钟, 哮喘/鼻炎)   ✅ M1 闭合：ST134 正式 z 检验, 双实现
     │                       │
     ├── MHC-free rg(双实现 ST126/134) ✅ ── 非HLA再裁定(ST131) ✅ ── healthspan锚(ST132, 名义级) ✅
     │
血液2型指数 ──MR/复制/跨祖源──> 过敏(3结局, 分级)              ✅ M2 闭合：ST133 遗传层剔除 + ST124 表型层剔除
     │                       │                                  ⚠ 残余：两层剔除结果分歧待解释(QC2, 文本级)
     ├── NHANES人群验证(年龄不变, 权重敏感) ✅ ── 反向界定(Steiger审计 ST135) ✅
     ├── HRS纵向构念(不测过敏, 已声明) ✅ ── 早警信号(均值匹配后存活) ✅
     ├── 机制: IL4悖论部分解析(DICE, 标题已同步) ✅(M4) ── OneK1K 邻接架构(非复制) ✅
     └── 转化: IL18(三平台) ✅ / TNFSF14(单工具, 已降级) ✅
```

证据链的每一个终止节点都是"经检验后收缩主张"而非"断点"，这是本稿与普通 MR 投稿的本质区别。

---

## 四、编辑总结

七轮评审弧线：过度包装（R10 前）→ 逐轮自我收缩 → 第六轮识别出两根新承重柱的形式化缺口 → **R12 修订以超额质量闭合**（M1 的 rg 对比在双实现下执行、M2 的遗传层剔除连同非 allele-aware 旧估计的废止声明一起交付）。本轮对全部五项 Major 做了数据级回溯核验，判定为真闭合。

当前文本的局限声明（12 条主局限 + SI Note 2）与证据分级标签体系，已达到本领域少见的成熟度——外审的意见将集中在 QC1（F 统计口径）与 QC2（两层嗜酸剔除分歧）这类可回答的问题，而非生存级质疑。

**决定：Minor Revision（接受轨道）。** 完成 QC1、QC2、文本 1、程序 1 四项必做项（预计 2 个工作日）+ 两项建议项后，本稿可投 Genome Medicine。若外审顺利，该稿的证据-主张匹配度足以在首轮外审中占据主动。
