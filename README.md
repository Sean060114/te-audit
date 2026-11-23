# TE 注释完整性检查工具链

一套回答**「这份转座子注释，哪些地方不能信」**的工具，配套的判据检验记录与真实案例。

不是又一个注释流程，而是**注释之后**的一层事实核查：拿 EDTA 自己产出的
`intact.gff3` / `TEanno.gff3` / 分类表互相对照，找出文件内部自相矛盾、
预期该有却没有、以及"数字看起来正常但其实测的是软件参数"的地方。

全部代码只依赖 Python 标准库。

---

## 一、模块地图

| 模块 | 回答什么问题 | 入口 | 测试 | 文档 |
|---|---|---|---:|---|
| **0 名称规范化** | EDTA 的 `DNA/hAT` 和 TEsorter 的 `hAT_TIR_transposon` 到底是不是一回事？ | `te-normalizer/normalize.py` | 90 | `te-normalizer/README.md` |
| **1 完整性检查** | 这份注释的哪些结论站不住？ | `te_check.py`（C1–C6）<br>`bench_multi.py`（多数据集横向） | 77 | `te_check_README.md`<br>`TE_case/报表_*.md` |
| **2 交叉验证** | EDTA 与独立注释差在哪，谁更可能对？ | `compare_araport11.py`（拟南芥 vs Araport11）<br>`compare_rice_te.py`（水稻 vs MSU RGAP7）<br>`candidates_novel.py` | — | `EDTA_vs_Araport11_评估报告.md`（含 k-mer 数据与发表级缺口附录）<br>`报表_水稻交叉验证与跨染色体基准.md` |
| **交付** | 把一份注释整理成可交付、可自校验的包 | `build_delivery_package.py`<br>`validate_delivery.py` | — | `delivery/`（不入库，可由 tar 包解出） |
| **上游贡献** | 把工具链里发现的命名缺陷回馈给 EDTA | `te-normalizer/make_alias_pr.py` | — | PR #655 |

模块 1 的六项检查：

| 编号 | 查什么 | 当前结论 |
|---|---|---|
| C1 | 哪一类别整类没有 identity 字段 | helitron，五份数据全一致（格式固有，非缺陷） |
| C2 | TIR 内部区过短（更像 MITE 或假阳性） | 有；短者不应与完整元件同列 |
| C3 | 反向重复匹配率能否当保守性指标 | **不能**，测的是 TIR-Learner 的 80% 阈值 |
| C4 | 注释里的元件是否都进了分类器 | 需按超家族分层，未分层会算出假缺口 |
| C5 | 库层面是否混入其他类别 | 有，且多数在过滤后被正确排除 |
| C6 | 高 LTR 一致性但缺结构域证据 | 有；高一致性 ≠ 近期插入 |

---

## 二、当前状态

**模块 0 已收尾**（代码与文档层面）：809 个别名、0 冲突，90 条测试；
能出规范化名 + 置信度 + 差异报告，并能把不一致率放回 EDTA 的分母上
（17 / 1593 = 1.07%，与 EDTA 自报值一致，见 `te-normalizer/examples/rate_demo.out`）。
**唯一未闭合项是上游 PR #655**（21 条别名）待维护者回复——不由我们决定，
代码侧能做的已经做完（补丁、理由、逐条排除项都写在 PR 描述里）。

判据在 **六份独立运行、两个科** 的数据上验证过（详见 `TE_case/`）：

| 数据集 | 物种 | intact 元件 | C4 已定超家族 LTR |
|---|---|---:|---|
| Chr1 | A. thaliana Col-0 | 218 | 40/40（100%） |
| col1m | A. thaliana Col-0 | 69 | 4/4（100%） |
| alyrata | A. lyrata | 20 | 3/3（100%） |
| rice_pre | O. sativa（单子叶） | 114 | 14/14（100%） |
| rice_full | O. sativa，Chr1 单条 | 1,669 | 187/198（94.4%） |
| rice_multi | O. sativa，Chr1+Chr12 | 2,733 | 350/373（93.8%）<br>逐条：Chr1 94.4% / Chr12 93.1% |

**C1 跨科成立**（六份一致）；**C3 指标无效**（阈值投影，不可跨数据集比较）；
**C4 分层后成立**（此前"阈值被证伪"的结论是分母算错，已推翻）；
**C6 跨科复现**（Chr1 / alyrata / 水稻若干例，成因需逐个区分）。

**run-to-run 已实测**：Chr1 在单染色体库与双染色体库下各跑一次，
六项检查逐项相同（LTR 218→218、C4 187/198、C6 23、C2 24、C3 中位 0.8182），
但同源注释覆盖从 37.83% 升到 39.42%（召回 +1.7 pp、精确率 −1.5 pp）。
**结构判据的噪声底≈0，召回/精确率则随库规模变动。**

**模块 2 已扩到水稻两条染色体**：以 MSU RGAP7（同一 IRGSP-1.0 组装）为独立真值集，
Chr1+Chr12 合并成一次运行（共用一套库），召回 87.0% / 87.4%——**密度升高 9.25 pp，召回没有下降**；
LTR 组双向一致度 78–85%，细桶 Gypsy 97.9% / 98.4%。EDTA 独有部分 85–90% 是 TIR + Helitron。
**同一段 Chr1 只是换了个库，召回就 +1.7 pp、精确率 −1.5 pp**——"库多大"比"哪条染色体"更关键。
详见 `TE_case/报表_水稻交叉验证与跨染色体基准.md`。

报告里保留了完整的自我修正记录，未删改——C3 的表述修正、C4 的分母修正、
`.mod.mod` 归因的两次撤回、EDTA 独有成分里 LTR 份额的分桶口径错（4.1% → 7.6%），
都原样留在文档里。错误的分母和它导致的错误结论，
本身就是这套方法最值得记下的部分。

---

## 三、复现入口

```bash
# 单份数据：六项检查（零依赖）
python te_check.py --edta-dir rice/edta_out --json rice/full_report.json

# 多数据集横向对比 + 判据检验（产出与 TE_case/bench_six.tsv 逐字节一致）
python bench_multi.py \
  "chr1=chr1_edta_out" "col1m=bench_archive/col1m.fa.mod" \
  "alyrata=bench_archive/alyrata.fa.mod" "rice_pre=rice/edta_out_preflight" \
  "rice_full=rice/edta_out" "rice_multi=rice/edta_out_multi" \
  --tsv TE_case/bench_six.tsv

# 模块 0：名称规范化（单列 / 混淆矩阵配对模式）
cd te-normalizer && python normalize.py examples/test_names.txt --verbose
python tests.py

# 跑一次新的 EDTA（需 WSL + conda 环境 edta，不要在 Windows 侧跑）
THREADS=16 bash run_bench_rice.sh              # 水稻 Chr1，约 73 分钟
PREFLIGHT=1 bash run_bench_rice.sh             # 先截 3 Mb 试跑，约 5 分钟
# 脚本已参数化，可跑任意基因组：
GENOME=/path/to/genome.fa TAG=Name SPECIES=others \
OUTDIR=$HOME/te_bench_x DEST=/path/to/out THREADS=16 bash run_bench_rice.sh
```

EDTA 运行必须在 WSL 内发起，`run_bench_rice.sh` 只在该环境的 bash 下可执行。

---

## 四、目录说明

```
te-normalizer/          模块 0：名称规范化（独立可跑，含自己的测试与文档）
TE_case/                模块 1 的基线、案例与本项目全部判据检验报告
  （已移至 archive/）报表_多基因组判据检验.md   拟南芥属内三份数据，C4 结论已推翻
  报表_水稻跨科判据检验.md     水稻跨科 + Chr1 重跑验证（判据结论以此为准）
  报表_水稻交叉验证与跨染色体基准.md  模块 2：单条库 → 双条同库，含库大小效应与 run-to-run 实测
  报表_结构注释质量基线.md     最早的单点基线
  基准_字段语义校准.md         建基准第 1 步：LTR 七条字段语义 + 三项独立复算
                              （配套 verify_ltr_structure.py，--tsv 输出逐元件判定；
                               TSD 99.0%、末端基序 98.8%、两翼相似度 95.0%）
  ltr_structure.tsv            420 条逐元件结构复算（四项门槛 + 失败项）
  基准_LTR逐元件台账.md        建基准第 2 步：420 条结构 LTR 与 MSU 逐元件对账
                              （配套 link_msu.py + ltr_msu_ledger.tsv；
                               A/B/C = 183/178/59，结构门槛落在逐条而非批次数，
                               配套 grade_ltr_abc.py + ltr_abc.tsv；
                               §11 用 MSU 基因流程做第三票，配套 gene_arbitration.py）
  方法学_已有QC工具与自建基准的分工.md  HiTE/TE_Bench/TEtrimmer/TE-Aid 各自能做什么
                              （含"LTR 侧与 EDTA 共用 LTR_retriever"的源码证据；
                               §2.0 更正：位点级前人做过，PolyA 2021）
  方法学_与PolyA与Carey的判据逐项对照.md  与三篇位点级文献的逐项对照：哪一条在重复、
                              哪一条没人做过（结论：只剩"跨来源 + 解析层 + 归因"三条）
  文献核对_PolyA置信度公式源码级.md   PolyA 的置信度公式用源码核对（PDF 被 429）：
                              `2^s/Σ2^s` 不除 λ；skip 下限硬编码 1%；两套置信度
                              不在同一尺度，百分数不可互比
  基准_TIR字段语义校准.md      建基准第 3 步：TIR/MITE 字段约定 + 结构证据强度
                              （配套 verify_tir_structure.py；1,824 条里仅 62 条硬证据）
  基准_TIR与MSU逐元件对账.md   建基准第 4 步：62 条 TIR 锚与 MSU 对账
                              （配套 link_msu.py --kind tir + tir_msu_ledger.tsv；
                               34 条可判、27 条一致 = 79.4%，28 条 MSU 无整段对应）
  基准_TIR身份核查_28条无对应元件.md  建基准第 5 步：那 28 条的定性
                              （Repbase/Dfam 双路排除 + 家族聚类 + MSU 基因模型；
                               3 条可命名，A 集 27 → 30）
archive/                已被取代的运行脚本（保留是为了留下失败记录）
bench_archive/          col1m / alyrata 的历史运行回传（仅此一份）
chr1_edta_out/          TAIR10 Chr1 基准运行回传
rice/                   水稻运行回传 + 报告 JSON（序列文件不入库）
ref/                    参考序列（不入库，可从 Ensembl Plants release-57 重下）
delivery/               交付包（不入库，用 make_deliverables.sh 重建）
papers/                 立项时读的文献
选题与发展规划_v3定稿.md   执行版：选题定稿（叫什么、不叫什么）、四阶段路线图、
                              产出清单、保研材料三件套（简历/陈述/面试稿）
选题说明书_详细论证版.md   立项论证，含最初的判据设想（v2，保留作存档）
```

`.gitignore` 里写明了哪些内容故意不入库、以及判断标准（重跑一次要多久）。

---

## 五、产出清单与处置

根目录 41 个条目，总体积约 230 MB（其中 `rice/` 98 MB、`ref/` 37 MB、`delivery/` 35 MB、
`TE_case/` 30 MB 都是序列、交付物或运行环境，代码与报告合计不到 1 MB）。

### 保留：核心代码（改了就坏）

| 条目 | 作用 |
|---|---|
| `te_check.py` `test_te_check.py` `bench_multi.py` | 模块 1 全部逻辑与 **77 条**测试 |
| `run_bench_rice.sh` | 现行 EDTA 运行脚本（参数化，含 `.mod.mod` 修复与 tee 日志） |
| `te-normalizer/` | 模块 0，独立可跑（90 条测试） |
| `compare_araport11.py` `compare_rice_te.py` `candidates_novel.py` `make_publication_tables.py` | 模块 2 |
| `build_delivery_package.py` `validate_delivery.py` `make_deliverables.sh` `make_union_softmask.py` `get_dfam_lib.sh` `run_validation.sh` | 交付包与验证 |
| `archive/` | 三版历史运行脚本 + 一份被推翻的判据报告，失败记录 |

### 保留：分析证据（报告里每个数字的依据）

| 条目 | 体积 | 为什么不能删 |
|---|---:|---|
| `chr1_edta_out/` | 1.7 M | TAIR10 Chr1 基准，重跑需 47 分钟 |
| `rice/edta_out*` | 入库 20 M | 水稻三条运行回传（Chr1 单条 / 预飞 / Chr1+Chr12），重跑需 73–124 分钟；序列文件不入库 |
| `bench_archive/` | 251 K | col1m / alyrata 历史运行，**仅此一份** |
| `TE_case/` | 30 M<br>（入库 <1 M） | 报告、基线、病例；`.venv`（29 M）只装了 `pyhmmer==0.12.3`，`classify_element.py` 的运行环境，重建命令见下 |
| `araport11_vs_edta.txt` `.per_class.tsv` | 12 K | 模块 2 报告的逐类数字来源 |
| `Araport11_TE_family_summary.tsv` `EDTA_B_level_candidates.tsv` `EDTA_candidate_novel.tsv` | 144 K | 模块 2 的三张成果表（家族汇总 / B 级质控 43 条 / 918 条分级，其中 B=43、C=4、D=871）。前两张是同一批 43 个元件的两种视图：B 级表坐标 1-based、带 TIR 错配细节，分级表坐标 0-based、带 reason 列，**起始坐标系统性差 1，不是数据冲突** |
| `chr1_report.json` `rice/full_report.json` `rice/preflight_report.json` `rice/multi_report.json` | 30 K | 四次运行的六项检查原始结果 |
| `rice/msu_vs_edta.json` `rice/msu_vs_edta.per_class.tsv` `TE_case/bench_six.tsv` | 30 K | 跨染色体交叉验证与六数据集对比的逐类数据（`bench_six.tsv` 已用现行数据重放，逐字节可复现） |

`TE_case/.venv` 若丢失或换机，重建只需两行（pyhmmer 有 py3.13 轮子，无需编译）：

```bash
python -m venv TE_case/.venv
TE_case/.venv/Scripts/pip install pyhmmer==0.12.3
```

注意 `classify_element.py` 里写死的 HMM 库路径指向 WSL 内的 TEsorter 安装
（`REXdb_protein_database_viridiplantae_v4.0.hmm`），换机后需同步改第 17–19 行。

### 保留：交付物与立项文档

| 条目 | 说明 |
|---|---|
| `delivery/`（35 M） | 交付包解包版，12 文件 + MD5SUMS，与下面那行 tar 包**逐字节相同**；`validate_delivery.py` 可直接校验 |
| `TAIR10.Chr1.TE_annotation_v1.tar.gz`（9.8 M） | 同一交付包的分发版，删掉 `delivery/` 后 `tar xzf` 即可还原 |
| `TE_case/报表_水稻交叉验证与跨染色体基准.md` | 模块 2 合并版（原两份合并）：单条库 → 双条同库、库大小效应 1.7 pp、run-to-run 实测 |
| `EDTA_vs_Araport11_评估报告.md` | 模块 2 报告（拟南芥）。含四个附录：k-mer 原始数据、发表级缺口、RepeatMasker 库处置、交付统计口径 |
| `te_check_README.md` | 模块 1 工具文档 |
| `选题说明书_详细论证版.md` | 立项论证与两段过程记录（附录 A 实战踩坑、附录 B 注释阶段四个发现，均为正文的一部分） |
| `papers/` `papers_digest.txt` | 19 篇立项文献全文与提要 |
| `EDTA_issue_台账.csv` | 上游 issue 与 PR 记录 |

### 不入库的外部数据（体积大，可重新下载）

| 条目 | 体积 | 怎么回来 |
|---|---:|---|
| `ref/TAIR10.Chr1.fa` | 30 M | `curl` Ensembl Plants release-57（命令见第 3 节） |
| `rice/Oryza_sativa.Chr1.fa` `Chr12.fa` | 70 M | 同上，逐条下载后 `tr -d '\r'` 清行尾、seqid 改为 `OsChr1` / `OsChr12` |
| `araport11_TE.gff3` | 9.2 M | Araport11 官方发布 |
| `ref/msu/osa1_r7.asm.repeat_masked.gff3.gz` | 7.5 M | MSU/UGA RGAP7（服务器慢，用 `curl -C -`） |

Chr1+Chr12 合并基因组不单独保存，需要时两秒重建（已验证与原文件 md5 相同）：
`cat rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa > rice/Oryza_sativa.Chr1_Chr12.fa`
（合计 70,802,779 bp，md5 `05b9072127587f084344991013ff3758`）。
原始下载包 `Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa.gz`（13.3 MB，
md5 `c9fd0909a4c6106a9b81d00399797eba`）也已删除——它只用于生成 `Chr1.fa`，
而 `Chr1.fa` 就在库里（43,270,923 bp，与包内序列长度一致）。

### 2026-11-17 清理记录

| 动作 | 对象 | 依据 |
|---|---|---|
| 合并 | `EDTA_LTR丢失溯源报告.md` + `TE_case/病例_7472bp_LTR_unknown.md` → `TE_case/病例_001_LTR_静默丢弃.md` | 同一元件（17,528,667–17,536,138）；一个记管线链路、一个记序列证据，分开写只看到半个因果 |
| 合并 | `EDTA_交叉验证报告.md` → `EDTA_vs_Araport11_评估报告.md` 附录 A–D | 该报告第一节结论已被推翻；其 k-mer 数据、发表级缺口、库处置、交付口径仍有效，不该随错误结论一起丢 |
| 归档 | `TE_case/报表_多基因组判据检验.md` → `archive/` | C4 结论已推翻（分母算错）；保留为 col1m / alyrata 的数据出处 |
| 删除 | `附录_实战记录.md` | 27/30 行已包含在 `选题说明书_详细论证版.md` 的附录中 |
| 删除 | `TE_case/bench_rice_pre.tsv` `bench_rice_full.tsv` | 分别是 3 数据集、4 数据集版本，已被五数据集版覆盖，零引用 |
| 删除 | `TE_case/te_check_report.json` | 旧 schema（无 `C4_ltr_tiers`、`severe_total=5`）且 `dir` 指向已不存在的 WSL 路径；同一数据的现版本是 `chr1_report.json` |
| 删除 | `papers/*.log`（19 个） | 抓取时的提示（"6 picture(s) were left out"），每个约 100 字节，无信息量 |

**第二轮（同日，逐条核验后执行）**——判断标准只有三条：内容重复、可由现存文件秒级重建、或只是下载缓存。

| 动作 | 对象 | 依据 | 省 |
|---|---|---|---:|
| 合并 | `报表_水稻独立真值集交叉验证.md` + `报表_水稻跨染色体基准.md` → `TE_case/报表_水稻交叉验证与跨染色体基准.md` | 同一套分析的上下两半，第二份每处对比都要写"上一版报告是……"，读者必须两份并看才知道 85.3% 与 87.0% 差在哪；而那条差（库大小效应 1.7 pp）恰是本轮最关键的发现。合并时按统一口径重算了 EDTA 独有成分，**发现原报告"LTR 只占 4.1%"只数了 Gypsy，实为 7.6%**，已在文中留档 | — |
| 删除 | `rice/Oryza_sativa.Chr1_Chr12.fa` | 与 `cat Chr1.fa Chr12.fa` **md5 完全相同**，重建两秒 | 68 M |
| 删除 | `rice/Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa.gz` | 只是下载缓存，唯一用途是生成 `Chr1.fa`，而后者在库；md5 已记在上文 | 13 M |
| 删除 | `rice/msu_vs_edta.multi.json` `.multi.per_class.tsv` | 与无后缀版 **md5 完全相同**，同一命令跑了两遍留下的副本 | 30 K |
| 删除 | `TE_case/bench_five.tsv` | 逐字节等于 `bench_six.tsv` 的前 5 行；`bench_six` 已用现行数据重放校验 | 0.6 K |
| 删除 | `附录B_anno阶段实测.md` | 与 `选题说明书_详细论证版.md` 附录 B 同文（134 行，差异仅表格对齐与下划线），同一逻辑已在第一轮处理过附录 A | 16 K |
| 保留 | `delivery/`（35 M）+ tar 包（9.8 M） | 两者逐字节相同，本可只留其一；`delivery/` 方便直接翻阅与跑 `validate_delivery.py`，用户选择都留 | — |
| 保留 | `TE_case/.venv`（29 M） | 看着像垃圾，实际是 `classify_element.py` 的运行环境；重建命令已写进第 5 节 | — |

## 六、已知空白

1. ~~**`.mod.mod` 双前缀的成因未定。**~~ **已查明。** EDTA.pl 已把基因组归一化成
   `X.fa.mod`；v2.2.2 的 `EDTA_raw.pl:440` 调用 LTR_retriever 时未传
   `-convert_seq_name`（默认 1），后端又归一化一次，`encode_fasta` 的输出路径是
   `"$genome.mod"`，于是得到 `X.fa.mod.mod`。**与序列名长短无关**，
   v2.3.0 起补上该参数后不再出现。详见 `TE_case/报表_水稻跨科判据检验.md` 5.3 节。
   当前仍用改名修复而不升级——五份数据都在 v2.2.2 下跑，升级会破坏可比性。
2. **水稻有已定超家族的 LTR 未进分类表**（Chr1 11 个 / Chr12 12 个，全为 Gypsy）。
   结构完整、超家族已定。原先怀疑与 1 同源，但归因查实后该推论已下调
   （短 ID 的归一化不改变序列内容），目前无解释。两份库下数量一致，不是随机波动。
3. **Chr1 的输入来源与旧运行无法逐字节核对。** 旧值取自 EDTA 自带 test 集，
   本次用 Ensembl release-57；两个计数（218 / 41）吻合，但原始文件已失。
4. **模块 2 只覆盖水稻两条染色体（19% 基因组）**，无中间密度点，
   且密度与组成同时变化（Chr12 的 Gypsy 份额约为 Chr1 的两倍），二因素无法分离。
   模拟基因组真值集与跨版本对比（2.2.2 / 2.3.0）均未开始。

## 七、运行环境

| 组件 | 版本 | 说明 |
|---|---|---|
| Python | 3.13 | 全部代码零第三方依赖 |
| EDTA | 2.2.2（conda env `edta`，WSL） | 全部数据均此版本，升级会让新旧数据不可比 |
| 参考组装 | IRGSP-1.0 / TAIR10 | 均取自 Ensembl Plants release-57 |
