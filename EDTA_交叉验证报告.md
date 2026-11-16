> ## ⚠️ 本报告第一节的结论已被推翻，请以 `EDTA_vs_Araport11_评估报告.md` 为准
>
> 本报告第一节写的「9.24% 不偏低」**是错的**。当时只有无库 k-mer 一种证据，把
> 「k=20 可重复碱基 12.96%」误当成了 TE 总量的上界。实际上 k=20 只能捕捉**近期**重复，
> 高度退化的古老 TE 残迹（20-mer 已不再精确重复）完全不计入，所以 12.96% 只是**下限**。
>
> 用 Araport11 官方 TE 注释（TAIR10 同坐标）做参照后，真实值是：
> **Chr1 的 TE 含量 16.70%（5.08 Mb），EDTA 注释 9.20%（2.80 Mb），碱基级召回仅 36.6%。**
>
> 第二节及之后的**数据本身仍然有效**（k-mer 分层、拷贝数区间、分类一致性），
> 只是对数据的解释需要按新报告修正。下面保留原文以便对照。

# EDTA 注释交叉验证结论（9.24% 是否偏低）

- 基因组：`chr1.fa.mod`，30,427,671 bp，单条染色体，GC 35.9%，AT 64.1%，N 164,359 bp
- EDTA v2.2.2，参数 `--step anno --anno 1 --species others -t 8`（默认**非** sensitive）
- EDTA 最终注释：**2,799,643 bp / 9.20%**，3,502 条

---

## 一、结论

**9.24% 不偏低，处于合理范围。** 真正的问题不是总量，而是**退化片段的召回率**和**个别类别的分类一致性**。

三条独立证据：

1. **无库 k-mer 法**（k=20，纯序列自比对）：全基因组可重复碱基（20-mer 出现 ≥2 次）占 **12.96%**，比 EDTA 的 9.20% 只高 3.8 个百分点。差额基本是 EDTA 按设计**不注释**的串联重复（端粒 `CCCTAAA` 阵列、卫星、rDNA）。
2. **按拷贝数分层**：中低拷贝（3–49 拷贝，TE 家族的典型区间）EDTA 覆盖 40–55%；高拷贝（≥50，串联重复区间）只覆盖 1.6–3.9% —— 这正是 EDTA 的正常行为（EDTA 不注释卫星，只少量输出 `repeat_fragment`）。
3. **单拷贝区被注释 1,594,615 bp**（占 EDTA 掩码 57%）：说明 EDTA 把大量**已无重复拷贝证据的退化 TE 片段**也标了出来，这是“宁多勿漏”的方向，不是漏注释。

所以：**总量没有系统性低估**；但中低拷贝 TE 有约一半的退化片段没被识别，这部分只在“发表级完整性”要求下才需要补。

---

## 二、关键数据

### 2.1 全基因组重复空间 vs EDTA

| 证据 | 碱基 | 占基因组 |
|---|---|---|
| 20-mer 出现 ≥2 次（可重复） | 3,943,463 | **12.96%** |
| 20-mer 出现 ≥3 次 | 2,044,470 | 6.72% |
| EDTA 注释 | 2,799,643 | **9.20%** |

### 2.2 按拷贝数分层（每碱基取覆盖它的 20-mer 最大拷贝数）

| 拷贝数 | 碱基 | 占比 | 被 EDTA 覆盖 | 覆盖率 |
|---|---|---|---|---|
| 2 | 1,898,993 | 6.24% | 401,290 | **21.1%** |
| 3–4 | 816,582 | 2.68% | 333,287 | **40.8%** |
| 5–9 | 455,259 | 1.50% | 251,162 | **55.2%** |
| 10–49 | 413,266 | 1.36% | 210,199 | **50.9%** |
| 50–99 | 138,310 | 0.45% | 5,437 | 3.9% |
| ≥100 | 221,053 | 0.73% | 3,642 | 1.6% |
| **合计** | **3,943,463** | **12.96%** | **1,205,017** | **30.6%** |

读法：
- **行 2–5（3–49 拷贝）= TE 家族主战场**，覆盖 41–55%，尚可但有提升空间
- **行 6–7（≥50 拷贝）几乎不覆盖**，属串联重复/卫星，符合 EDTA 设计
- 行 1（恰好 2 拷贝，1.90 Mb）覆盖率仅 21%，其中混杂了**片段重复（segmental duplication）**和基因家族，不能全算 TE 漏注释

### 2.3 分类一致性（EDTA 自带混淆矩阵，`chr1.fa.mod.EDTA.anno/`）

| 类别 | 条目 | 不一致数 | 不一致率 | 备注 |
|---|---|---|---|---|
| **DNA/Tc1_Mariner** | 48 | 7 | **14.89%** | 7 条被判成 PIF_Harbinger ← **需处理** |
| DNA/CACTA（嵌套） | 13 | 2 | 15.38% | 样本小，影响有限 |
| LTR/Copia | 68 | 2 | 2.94% | 2 条被判成 Mutator |
| LTR/Gypsy | 623 | 2 | 0.32% | 良好 |
| 其余各类 | — | — | <1% | 良好 |

---

## 三、发表级注释还差什么

按必要性排序：

### 必做

**① 修 `rename_LTR_skim.pl` 的坐标顺序 bug 并重跑 LTR**
影响 1/41 条 intact LTR（`1:17536138..17528667`，7,472 bp）。对发表级注释，`intact LTR` 集会用于 **LTR 插入年龄分布**（5'/3' LTR 差异），少一条就是数据缺失。补丁见 `EDTA_LTR丢失溯源报告.md` 第五节。

**② 处理 Tc1_Mariner 14.89% 误分类**
用 `--sensitive 1` 重跑 TIR 模块（TIR-Learner 会做更严格的末端/靶位点验证）。绝对数量只有 72 条，但发表时审稿人会看混淆矩阵。

**③ 产出下游交付文件**
- 软屏蔽基因组（当前只有硬屏蔽 `MAKER.masked`，新屏蔽 1.70 Mb / 5.6%）
- BED12
- 家族级 bp/占比/拷贝数表

### 建议做

**④ 用独立库（Dfam/RepBase）跑一次 RepeatMasker 交叉验证**
上面的 k-mer 法只能证明“总量不离谱”，不能证明“类别对”。`run_validation.sh` 已备好，跑完把 `.tbl` 结果拿来对比，可以给出“EDTA vs Dfam”的类别级差异表——这是发表时最常被要求补充的对照。

**⑤ 决定是否补注释退化片段**
若要求高完整性，可在 EDTA 基础上用 `RepeatMasker -lib EDTA.TElib.fa -e ncbi` 放宽阈值重跑注释层（不动 TElib），把中低拷贝区覆盖率从 ~50% 提上去。代价是假阳性上升，需要人工抽查。

---

## 四、方法学说明（写 methods 时可用）

k-mer 法细节：
- k = 20，base-4 编码（A/C/G/T → 0/1/2/3），滚动编码为 int64 整数
- 含 N 的 k-mer（窗口内任一非 ACGT）剔除，占 0.54%
- 唯一 k-mer 数 28,684,567（30,427,652 个窗口中）
- 每碱基重复度 = 覆盖该碱基的所有 k-mer 的最大出现次数
- EDTA 掩码由 `chr1.fa.mod.EDTA.TEanno.gff3` 构建，排除 `target_site_duplication`、`long_terminal_repeat`、`primer_binding_site` 三类附属特征（避免重复计数）

局限：
- k=20 只能捕获较近期重复；20% 分歧的两拷贝共享精确 20-mer 的概率仅 0.8²⁰ ≈ 1.2%，因此**高分歧的古 TE 不被 k-mer 法感知**，本方法对“漏注释”的判断是**保守下限**
- 2 拷贝层混有片段重复与多拷贝基因，不能全归为 TE

---

## 五、环境问题：conda 版 RepeatMasker 的库是空的

首次执行 `run_validation.sh` 时 RepeatMasker 直接退出，报错：

```
Using Master RepeatMasker Database: .../Libraries/RepeatMaskerLib.h5
Title    :
Version  :
Date     :
Families :
Species "homo sapiens" is not known to RepeatMasker.
```

排查结果（`~/miniconda3/envs/edta/share/RepeatMasker/Libraries/`）：

| 文件 | 实际状态 | 影响 |
|---|---|---|
| `Dfam.h5` | **11 字节**，内容就是字符串 `placeholder` | conda 包不附带 Dfam，需自行下载 |
| `RepeatMaskerLib.h5` | 权限 `?????????`，读取报 `Input/output error` | 主库不可用 → `Families:` 为空 |
| `RepeatMasker.lib` | 32 MB / 19,768 条，但**植物类群条目 = 0** | 全是昆虫、鱼、哺乳动物（Sitophilus、Halyomorpha、Danio…），对植物基因组无效 |

因此对植物基因组，这套库**根本没法用来做交叉验证**（EDTA 自身不受影响，因为它用的是自己建的 `EDTA.TElib.fa`）。

### 解决方案

Dfam 3.9 提供 `Dfam-RepeatMasker.lib.gz`（**15.8 MB**），是 FASTA 格式的共识库，可直接 `-lib` 传给 RepeatMasker 4.1.5，**无需升级版本、无需下载 GB 级 h5**：

```bash
mkdir -p ~/edta_run/dfam_lib && cd ~/edta_run/dfam_lib
wget https://www.dfam.org/releases/Dfam_3.9/families/Dfam-RepeatMasker.lib.gz
gunzip -k Dfam-RepeatMasker.lib.gz
grep -c '^>' Dfam-RepeatMasker.lib        # 确认家族数
```

> 注意：Dfam 4.0 的 `dfam40.*.h5.gz` 分区格式要求 **RepeatMasker ≥ 4.2.4**，当前环境是 4.1.5，不要用。若要升级，需 `conda install -c bioconda "repeatmasker=4.2.4"`，但会动到 EDTA 依赖的环境，发表阶段不建议。

---

## 六、交付文件的统计口径修正

`make_deliverables.sh` v1 的家族表按“逐条长度相加”，嵌套元件会重复计数，合计 **10.10%**，与软屏蔽实测 **9.20%** 对不上——投稿时这是硬伤。v2 同时给出两种口径：

| 口径 | 含义 | 合计 | 用途 |
|---|---|---|---|
| `bp_redun` | 逐条相加（含嵌套/重叠重复计数） | 3,073,010 bp = **10.10%** | “元件总数”描述 |
| `bp_union` | 同类内合并去重 | — | “该类占基因组比例” |
| 全局去重 | 所有类别合并去重 | 2,799,643 bp = **9.20%** | **与软屏蔽逐字节一致，引用此值** |

修正后的家族表（`bp_union` 口径，与 `TEanno.sum` 的 `bpMasked` 基本吻合）：

| 类别 | 条目 | 逐条 bp | 去重 bp |
|---|---|---|---|
| nonTIR/helitron | 1136 | 998,041 | 996,497 |
| LTR/Gypsy | 885 | 942,977 | 847,763 |
| LTR/Copia | 266 | 412,091 | 282,379 |
| TIR/Mutator (DTM) | 296 + 25 MITE | 218,081 | 211,542 |
| TIR/CACTA (DTC) | 235 | 170,729 | 166,576 |
| LINE/L1 | 186 | 120,970 | 120,970 |
| TIR/PIF_Harbinger (DTH) | 303 + 17 MITE | 121,263 | 110,342 |
| SINE/tRNA | 25 | 31,149 | 31,149 |
| TIR/Tc1_Mariner (DTT) | 74 | 24,376 | 24,376 |
| LTR/unknown | 6 | 15,957 | 8,485 |
| TIR/hAT (DTA) | 19 | 4,115 | 4,115 |
| unknown | 29 | 4,121 | 4,121 |

TEsorter 三字母代码对应关系（按条目数与 `TEanno.sum` 逐类核对）：DTA=hAT，DTC=CACTA，DTH=PIF_Harbinger，DTM=Mutator，DTT=Tc1_Mariner。

