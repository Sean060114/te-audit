# 方法学：已有的注释质检工具，与自建基准的分工

> 写这份文档的原因：项目做到第 5 步时，"HiTE / Earl Grey 不是已经有注释质量检查了吗"
> 是个必然会再被问到的问题（审稿人也会问）。这里把调研结果和边界一次写清，
> 免得每轮重新论证。

## 一、结论

**有，而且比这个问题里提到的多——七个。但它们解决的都不是本项目的问题。** 更关键的一条：
**HiTE 的 LTR 结构判定与 EDTA 调用的是同一份实现**（都是 LTR_retriever），
拿它给 EDTA 的 LTR 做质检，等于自己检查自己。

## 二、七个现成工具，四种粒度

| 工具 | 输出粒度 | 需要什么 | 本项目的可用性 |
|---|---|---|---|
| **TE-Aid**（Mobile DNA 2022 配套） | **逐条可视化**：blastn hits+divergence / 覆盖度 / self dotplot / TIR·LTR 建议 / ORF / domain | blastn + R + EMBOSS，一条 consensus + 基因组 | **可用**（本机缺 blastn/R，需 WSL）。但它只出图，判定仍靠人 |
| **TEtrimmer**（Nat Commun 2025） | **family 级自动分级**（5 项：末端重复 / 分类状态 / MSA 序列数 / 全长 BLASTN hits / PFAM domain）+ GUI | blastn + mafft + cd-hit + bedtools + hmmer | **可用**，本机缺依赖 |
| **TEtrimmer** `TE_anno_benchmarker.py` | **注释级混淆矩阵** + Sankey（1,372 行，改编自 EDTA `lib-test.pl`） | 两个 **RepeatMasker `.out`**（query + 人工 curated 参考） | **可用且值得做**：MSU 的 `repeat_masked.gff3` 可当参考（需转 `.out`） |
| **MCHelper** | family 级自动 curation + 人工检查模块（内嵌 TE-Aid 图） | 同上 | 不适用（评的是家族） |
| **HiTE** `module/lib_evaluation.py` | 碱基覆盖（TP/FP/FN → sensitivity / precision / FDR / F1） | **standard_lib**（gold standard 共识库）+ 待测库 + 基因组 | **不可用**：无 gold standard |
| **TE_Bench**（Mobile DNA 2026） | 同上，但带模拟真值、按 TE 类与嵌套分层出图 | Snakemake + SLURM + Dfam（或 Repbase） | 框架可用，**本机跑不了**（需 HPC） |
| **Earl Grey** BEAT / `--curate` | family 级 consensus 边界与冗余 | 人工指定要 curate 的家族 | 不适用 |
| **LAI**（LTR_retriever） | 组装质量（用 LTR 注释反推） | LTR 注释 | 方向相反：它是注释的**消费者** |

四种粒度要分清：

- **库级**（HiTE 评估、TE_Bench）：回答"这套库整体好不好"，输出一个 F1。
- **注释级·大类**（TEtrimmer benchmarker）：回答"这套注释在 LTR / LINE / SINE /
  DNA-TIR / Helitron / Not-TE 六类上的碱基混淆矩阵长什么样"，出 Sankey 图。
- **family 级**（TEtrimmer 分级、Earl Grey curate、MCHelper）：回答"这个家族的
  consensus 两端对不对、拷贝够不够"。
- **逐条人工**（TE-Aid）：不判定，只把证据画成 4 张图给人看。

**位点级自动判定仍然没有。** 原因不是没人做，而是位点级的真值只能来自人工 curation
或模拟基因组——而这两条路在本项目里分别被"工作量"和"没有水稻模拟真值"堵住。

### 2.1 最接近的一步：TEtrimmer 的分级，但它评的是 consensus

TEtrimmer 给每个输出 consensus 分评估等级，依据五项：末端重复是否存在、分类状态、
MSA 序列数、全长 BLASTN hits 数、PFAM domain。这是目前**自动化程度最高的"质控"**。

差异在对象：它给**一条家族代表序列**打分，我们给**基因组上每一个插入位点**定性。
一条 consensus 可以等级很高，而它在基因组上的 32 个拷贝里有 3 个被标成了别的类别
（本项目实测：`TE_00001527` 就是这种情况）。**consensus 的质量不传递到位点。**

### 2.2 TEtrimmer benchmarker 能做什么、不能做什么

`--query` 与 `--reference` 都要 RepeatMasker `.out`，`reference` 必须是"人工 curated
库"跑出来的。它把两侧的 bed 取交集与补集（`bed_merge_and_cross_intersect`），按类别
汇总碱基，出混淆矩阵与 Sankey。

- **能做**：跨大类的错误能看见。23 条 AnacB 冲突（DNA-TIR ↔ LTR）正好是跨大类，
  这张图能反映。
- **不能做**：类别被压成六类，**超家族级错误全部消失**——Gypsy 标成 Copia、DTM 标成
  DTC，在矩阵里都是同一个格子。而本项目的 7 条 TIR 冲突、185 条 LTR 一致率里，
  超家族才是主战场。
- **不能做**：边界。±5 bp 的 TSD 归属差异在这套指标里是零。

**它的正确用法是当"对外背书"，不是当"我们的方法"。** 拿 MSU 的 `repeat_masked.gff3`
转成 `.out` 当 reference，跑一次，得到 Sankey + 混淆矩阵——这张图有现成工具与
Ou et al. 2019 的方法学背书，审稿人认。而元件级的结论仍由我们的台账提供。

## 三、HiTE 的评估到底在算什么

`lib_evaluation.py`（282 行）读一遍就清楚：

1. 把 `--standard_lib` 与 `--test_lib` 各自 blastn 到基因组，按
   `--coverage_threshold`（0.80 / 0.95 / 0.99）保留全长拷贝；
2. 把每条染色体**切成 100 kb 段**，把两侧的片段登记到段里；
3. TP = 真值标了、待测也标的碱基；FN = 真值标了、待测没标；FP = 反之；
4. 输出 `sensitivity / precision / FDR / F1`。

两个性质决定了它不能用来做本项目的质检：

- **它只比位置覆盖，不比标签。** 一条元件被 EDTA 标成 DTM、真值是 DTC，在这套
  指标里是完整的 TP。而本项目的核心争议恰恰是**标签**（23 条 AnacB、28 条 TIR）。
- **粒度是 100 kb 段 + 碱基覆盖。** 用它评 4–9 kb 的元件，分辨率差两三个数量级；
  它连"这一条元件边界偏了 5 bp"都表达不出来。

## 四、独立性：LTR 同源，TIR 异构

这条比"粒度不匹配"更要紧，因为它决定"能不能当第三方"。

| | EDTA | HiTE | 是否独立 |
|---|---|---|---|
| LTR 结构判定 | LTR_FINDER + LTRharvest → **LTR_retriever**（`EDTA.pl` L290–295 直接 `command -v LTR_retriever`，L475 传路径给 `EDTA_raw`） | LTR_FINDER + LTRharvest → **LTR_retriever**（`judge_LTR_transposons.py` L98 `run_judge_LTR_detection`，内部调 `run_LTR_detection` / `run_LTR_retriever` / `run_FiLTR`） | **否，同一实现** |
| TIR 结构判定 | **TIR-Learner** | 自己的 TSD 搜索（`multi_process_tsd_v1`）+ **itrsearch** + 侧翼比对 | **是，不同实现** |

所以：

- **LTR 侧**，HiTE / EDTA / 任何调 LTR_retriever 的流程，共用同一套 TSD + 两翼
  相似度 + motif 判据。它们会**一起**在同一个地方出错——比如刚查实的 5 条"边界落在
  2 nt 交替微卫星里、TSD 有多个等效位置"，三个工具会犯同样的错。这类**判据级盲区**
  只能靠换**证据类型**（同源比对、基因模型），换工具没用。
- **TIR 侧**，HiTE 是真正独立的第二实现。这是目前唯一能给那 28 条无对应元件
  提供"第二个结构判据"的现成资源。**值得在 WSL 里跑一次**（本机缺 blastn 与
  itrsearch，沙箱禁 `wsl.exe`）。

## 五、那本项目自建的到底是哪一块

自建的部分只有一件是现成工具替代不了的：**跨证据类型的位点级仲裁**。

| 证据 | 来源 | 是否与 EDTA 同源 |
|---|---|---|
| 结构自证 | 我们自己的复算（TSD / motif / 两翼相似度） | 判据同源（LTR_retriever 的规则），但**实现独立**，能查实现 bug |
| 同源比对的位点标签 | MSU `repeat_masked.gff3`（不同流程、不同证据） | **不同源** |
| 转座酶基因模型 | MSU `all_models.gff3` | **不同源** |
| 家族性 / 拷贝数 | 我们的聚类 | 不同源（但证据弱） |

"实现独立"和"证据独立"是两件事。前者能查出 bug——motif 那次链方向错误就是靠
逐元件复算查出来的，任何接受"97.6% 通过"的人都不会去查。后者才能查出**判据**错。

## 六、该借用的，不该重写的

按性价比排序，四条：

1. **TEtrimmer 的 benchmarker 跑一次（最值得做）**。把 MSU `repeat_masked.gff3` 转成
   RepeatMasker `.out` 当 reference，EDTA 的 `TEanno.gff3` 当 query，得到 Sankey +
   六类混淆矩阵。**理由不是它更准，而是它有现成工具与 Ou et al. 2019 的方法学背书**，
   审稿人不会质疑。我们的元件级台账负责给出它给不出的那一层。格式转换正好是
   `te-normalizer` 的活儿。
2. **TIR 侧跑一次 HiTE**，当第二结构判据裁决 28 条——唯一有增量的现成结构证据。
3. **评估框架不要自己写。** 要跑分就用 TE_Bench（Snakemake + SLURM，需 HPC）；
   要做单条人工核查就用 TE-Aid（4 张图，标准做法）。两者都不必重实现。
4. **库级指标（F1、Perfect 比例）不要自建。** 我们没有 gold standard，自建的库级
   指标一定是自证的，写进论文会被一眼看穿。LAI 也已在 LTR_retriever 里，直接用。

## 七、一句话回答审稿人

> 现成工具分四档：TE-Aid 出图给人看、TEtrimmer 给 consensus 分级并出六类混淆矩阵、
> HiTE 与 TE_Bench 做库级 F1、Earl Grey / MCHelper 做 family curation。**没有一档
> 判定"基因组上某一条插入的标签是否成立"**，且前三档都需要 gold standard 或人工
> curated 参考。本研究要的正是这一档，而在本项目条件下不存在可用的 gold standard
> （Repbase 商业授权、Dfam 无水稻家族）。因此自建的是"跨证据类型的位点级仲裁"；
> 注释级混淆矩阵借 TEtrimmer benchmarker 出图，库级指标不重复实现。

## 八、出处

- HiTE：Nat Commun 2024, `s41467-024-49912-8`；仓库 `CSU-KangHu/HiTE`（master）。
  `module/lib_evaluation.py` 282 行、`module/judge_LTR_transposons.py` 392 行、
  `module/judge_TIR_transposons.py` 212 行。
- TEtrimmer：Nat Commun 2025, `s41467-025-63889-y`；bioRxiv 2024.06.27.600963；
  仓库 `qjiangzhao/TEtrimmer`，`tetrimmer/TE_anno_benchmarker.py` 1,372 行
  （类别映射取自 EDTA `lib-test.pl`，v1.0.1 2024-03-21）。
- TE-Aid：Mobile DNA 2022（"A beginner's guide to manual curation of TEs" 配套）；
  仓库 `clemgoub/TE-Aid`。输出 4 图：hits+divergence / coverage / self dotplot /
  TIR·LTR + ORF + domain。
- Earl Grey：Mol Biol Evol 2024, `msae068`；bioRxiv 2022.06.30.498289。
- TE_Bench：Mobile DNA 2026, `s13100-026-00405-z`；仓库 `hkania/TE_Bench`。
- 上述结论均据源码与论文原文，未凭印象。
