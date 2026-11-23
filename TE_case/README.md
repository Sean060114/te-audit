# TE_case：判据检验与基准建设

> 这里放的是"每一个数字的依据"。12 份文档**不是并列的**，有阅读顺序。
> 对外说明看 [../README.md](../README.md)。

---

## 一、只想看结论：读三份

| 顺序 | 文档 | 一句话 |
|---:|---|---|
| 1 | [基准_LTR逐元件台账.md](基准_LTR逐元件台账.md) | 420 条结构 LTR 与 MSU 逐条对账，A/B/C = 183/178/59 |
| 2 | [基准_TIR与MSU逐元件对账.md](基准_TIR与MSU逐元件对账.md) | 62 条 TIR 锚只有 34 条可判、27 条一致 |
| 3 | [方法学_与PolyA与Carey的判据逐项对照.md](方法学_与PolyA与Carey的判据逐项对照.md) | 哪一条在重复前人、哪一条没人做过 |

读完这三份能回答："这套判据凭什么成立，与文献比新在哪。"

---

## 二、建基准五步（按顺序）

每一步的产出都是下一步的输入，跳读会看到未定义的数字。

| 步 | 文档 | 脚本 | 产物 | 结论 |
|---:|---|---|---|---|
| 1 | [基准_字段语义校准.md](基准_字段语义校准.md) | `verify_ltr_structure.py` | `ltr_structure.tsv` | LTR 七条字段语义；三项独立复算 TSD 99.0% / motif 98.8% / 相似度 95.0% |
| 2 | [基准_LTR逐元件台账.md](基准_LTR逐元件台账.md) | `link_msu.py` `grade_ltr_abc.py` | `ltr_msu_ledger.tsv` `ltr_abc.tsv` | 420 条逐元件分级；§11 用 MSU 基因流程做第三票 |
| 3 | [基准_TIR字段语义校准.md](基准_TIR字段语义校准.md) | `verify_tir_structure.py` | `tir_structure.tsv` | 三条与 LTR 相反的约定；1,824 条里仅 62 条硬证据 |
| 4 | [基准_TIR与MSU逐元件对账.md](基准_TIR与MSU逐元件对账.md) | `link_msu.py --kind tir` | `tir_msu_ledger.tsv` | 34 条可判、一致 27；28 条 MSU 无整段对应 |
| 5 | [基准_TIR身份核查_28条无对应元件.md](基准_TIR身份核查_28条无对应元件.md) | `tir_family_cluster.py` `gene_occupancy.py` `pack_mule_test.py` | `tir_family.tsv` | 那 28 条不同质：3 条可命名、5 条裹宿主基因、20 条待定 |

**五步里最有价值的是第 1 步与第 3 步**——它们确定了字段语义，而后面所有比较都
建立在这上面。第 1 步之前，跨工具比较的错误率是"读错就全错"级别的。

---

## 三、方法学三份

| 文档 | 回答什么 |
|---|---|
| [方法学_已有QC工具与自建基准的分工.md](方法学_已有QC工具与自建基准的分工.md) | 七个现成工具分四档粒度；**HiTE 的 LTR 判定与 EDTA 共用 LTR_retriever**（源码证据） |
| [方法学_与PolyA与Carey的判据逐项对照.md](方法学_与PolyA与Carey的判据逐项对照.md) | 位点级前人做过（PolyA 2021）；本项目只剩三条差异；九维度对照表 |
| [文献核对_PolyA置信度公式源码级.md](文献核对_PolyA置信度公式源码级.md) | PolyA 的公式用源码核对：`2^s/Σ2^s` 不除 λ；与 Carey 的置信度不在同一尺度 |

第 1 份的结论比"功能重叠"更硬：**拿 HiTE 质检 EDTA 的 LTR 是同源自证**，因为两条
流程的 LTR 结构判据是同一份实现（LTR_FINDER + LTRharvest + LTR_retriever）。

---

## 四、报表：判据在六份数据上的检验

| 文档 | 数据 | 用途 |
|---|---|---|
| [报表_结构注释质量基线.md](报表_结构注释质量基线.md) | Chr1 | 最早的单点基线 |
| [报表_水稻跨科判据检验.md](报表_水稻跨科判据检验.md) | 水稻 + Chr1 重跑 | 跨科检验；`.mod.mod` 成因查明 |
| [报表_水稻交叉验证与跨染色体基准.md](报表_水稻交叉验证与跨染色体基准.md) | 水稻 Chr1 / Chr1+Chr12 | 模块 2；库大小效应 1.7 pp |

三份报表都保留了完整的自我修正记录，未删改——C3 的表述修正、C4 的分母修正、
`.mod.mod` 归因的两次撤回、EDTA 独有成分里 LTR 份额的分桶口径错（4.1% → 7.6%），
都原样留着。**错误的分母和它导致的错误结论，本身就是这套方法最值得记下的部分。**

---

## 五、台账（tsv）说明

所有台账都是 UTF-8、LF 行尾、制表符分隔，可被 pandas 直接读。

| 文件 | 行数 | 一行是什么 | 生成方式 |
|---|---:|---|---|
| `ltr_structure.tsv` | 420 | 一个 LTR 元件的四项结构复算与失败项 | `verify_ltr_structure.py --tsv` |
| `ltr_msu_ledger.tsv` | 420 | 一个 LTR 元件与 MSU 重复注释的判定（`verdict_code`） | `link_msu.py --kind ltr` |
| `ltr_abc.tsv` | 420 | 一个 LTR 元件的 A/B/C 分级与 `structure_fail` | `grade_ltr_abc.py` |
| `tir_structure.tsv` | 1,824 | 一条 TIR 的臂长、实测相似度、分级 | `verify_tir_structure.py` |
| `tir_msu_ledger.tsv` | 62 | 一条强/中档 TIR 锚与 MSU 的判定 | `link_msu.py --kind tir --tier 强,中` |
| `tir_family.tsv` | 62 | 一条锚的家族大小（家族性聚类结果） | `tir_family_cluster.py` |
| `bench_six.tsv` | 6 | 一份数据集的六项检查结果 | `bench_multi.py --tsv` |

`verdict_code` 的取值与含义见 `link_msu.py` 头部注释，共九种：
`agree_superfamily` / `agree_class` / `conflict_class` / `msu_ltr_only` /
`msu_nonltr_label` / `ed_unknown_msu_resolved` / `both_unknown` / `insufficient_cov` / `empty`。

---

## 六、脚本

| 脚本 | 做什么 | 依赖 |
|---|---|---|
| `verify_ltr_structure.py` | LTR 的 TSD / motif / 两翼相似度三项独立复算 | 标准库 + FASTA |
| `verify_tir_structure.py` | TIR 的字段校准、臂长分级、E-value | 标准库 + FASTA |
| `link_msu.py` | EDTA 元件与 MSU 重复注释逐条对账（`--kind ltr\|tir`） | 标准库 + MSU gff |
| `grade_ltr_abc.py` | 结构判定 × 跨源判定 → A/B/C | 标准库 |
| `gene_occupancy.py` | 元件内是否有完整基因模型（含零模型） | 标准库 + MSU `all_models.gff3` |
| `gene_arbitration.py` | 用基因流程给重复标签做第三票（含标签审计） | 标准库 + 同上 |
| `tir_family_cluster.py` | TIR 家族性聚类（k-mer + 无缺口比对） | 标准库 |
| `pack_mule_test.py` | Pack-MULE 判别（旁系位点检索）—— **判据不成立，见文档** | 标准库 + FASTA |
| `classify_element.py` | 基于蛋白域的元件分类 | pyhmmer（见 `docs/维护记录.md` 重建命令） |
| `intact_baseline.py` | 结构注释质量基线统计 | 标准库 |

`pack_mule_test.py` 保留的原因：它的**基线**（随机基因 50% 是单拷贝）推翻了原判据。
这个负结果比一个"看起来能跑"的判据更有价值。

---

## 七、病例

[病例_001_LTR_静默丢弃.md](病例_001_LTR_静默丢弃.md)：一个 7,472 bp 的 LTR 元件
在 EDTA 的 LTR 流程里被判出、却在分类环节消失。从管线链路与序列证据两侧落到同一
个坐标（17,528,667–17,536,138），是"静默丢弃"这条判据的来源。

---

## 八、复现顺序

```bash
# 第 1 步（需要基因组 FASTA）
python TE_case/verify_ltr_structure.py \
  --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
  --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa \
  --tsv TE_case/ltr_structure.tsv

# 第 2 步
python TE_case/link_msu.py --kind ltr \
  --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
  --msu ref/msu/osa1_r7.asm.repeat_masked.gff3.gz \
  --tsv TE_case/ltr_msu_ledger.tsv
python TE_case/grade_ltr_abc.py --structure TE_case/ltr_structure.tsv \
  --ledger TE_case/ltr_msu_ledger.tsv --tsv TE_case/ltr_abc.tsv

# 第 3-5 步见各文档末尾的复现命令
```

MSU 的 `repeat_masked.gff3.gz` 与 `all_models.gff3.gz` 均来自 rice.uga.edu，
下载方式见 `docs/维护记录.md` 第五节。
