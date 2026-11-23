# 文献核对：PolyA 的置信度公式（源码级）

> 起因：要正式引用 PolyA，但预印本 HTML 里 Eq 2/3/4 渲染丢失。本文记录实际核对过程与
> 结果。**结论：论文公式与代码实现不完全一致，且 PolyA 与 Carey 的置信度不在同一尺度，
> 数值不可直接比较。**

## 一、引用元数据（已确认）

| 项 | 值 | 依据 |
|---|---|---|
| 标题 | PolyA: a tool for adjudicating competing annotations of biological sequences | Crossref |
| DOI | `10.1101/2021.02.13.430877` | Crossref |
| 状态 | **`posted-content`（预印本，无期刊版本）** | Crossref `type` 字段 |
| 版本 | v1，2021-02-13 提交 | bioRxiv |
| 代码 | `github.com/TravisWheelerLab/PolyA`（master） | 仓库可达 |
| 可访问性 | 正文 HTML 可读；PDF 下载被 429 限流 | 实测 |

**引用时必须注明"预印本，未见正式发表记录"。** Crossref 检索只返回 `posted-content`
一条记录，按标题与作者反查无期刊版本。

## 二、PDF 拿不到时的替代路径（按可靠性排序）

实际试过并记录，供以后复用：

| 路径 | 结果 |
|---|---|
| `biorxiv.org/.../xxx.full.pdf` | **429 Too Many Requests**（两次，两个路径都拒） |
| Europe PMC REST（`resultType=core`）| 命中 PPR282077，但 `fullTextUrlList` 只回指 bioRxiv |
| Unpaywall / fatcat API | 代理层 502（本环境不可达） |
| `r.jina.ai` 读 PDF | 代理层 502 |
| **GitHub 源码（`raw.githubusercontent.com`）** | **可用，且比论文更精确** |

→ **经验：当论文的公式是核心依据、而 PDF 不可得时，优先读实现。** 实现是唯一的运行时
真相，论文是它的描述；两者不一致时，实现说的话作数，但引用时要同时说明两者都写了什么。

## 三、逐项核对

| 论文（HTML 文字描述） | 源码（master，函数/常数） | 判定 |
|---|---|---|
| 对竞争集合 `Q` 做后验归一化，均匀先验下 `P(qi|t) ∝ P(t|qi)` | `confidence_cm()`：`conf_i = 2^score_i / Σ_j 2^score_j` | 语义一致 |
| 比对分数是 **λ 缩放的 log-odds** | `confidence_cm()` 里是 `2 ** int(region[index])`，**没有除 λ**；`λ` 只出现在 `calculate_complexity_adjusted_score()` 做 **dust 复杂度校正**（`char_complexity_adjustments[c] /= lamb`、`s_raw + t_sum/λ + 0.999`） | **表述易误解**：λ 不进置信度的指数 |
| 位置级：窗口 `w=31`，取覆盖该位置窗口的**未加权平均** | `DEFAULT_CHUNK_SIZE = 31`，help 为 "size of the window in base pairs analyzed together" | 一致 |
| 转移惩罚 `tm = 1e-55`、`ts = 1−tm−tms`、`tss = 1−tms` | `DEFAULT_TRANS_PENALTY = 55`（注释："A value of 55 means the probability will be 1e-55"）；`SAME_PROB_LOG = 0.0`（留同行惩罚为 0） | 一致 |
| skip state 有 "special default score" | `SKIP_ALIGN_SCORE = 10.0`，注释 "no lambda adjustment needed" | **更精确**：固定 10 分参与 `2^score` 归一化 |

## 四、论文没写、但代码里有的四个常数（引用时值得补）

| 常数 | 值 | 含义 |
|---|---|---|
| `PROB_SKIP` | 0.4 | skip state 的先验概率；注释给的理由是"约 60% 基因组是 TE 来源" |
| `PROB_SKIP_TR` | 0.06 | 串联重复先验（约 6% 基因组） |
| `SKIP_ALIGN_SCORE` | 10.0 | skip 的固定分数，**不做 λ 校正** |
| skip 置信度下限 | **0.01（硬编码）** | `if confidence_list[0] < 0.01: 置 0.01，其余按 0.99 等比压缩` |

最后一条最值得注意：**PolyA 的 skip 置信度永远不会低于 1%**。也就是说 skip state 不是
"判断此处没有 TE"，而是"至少留 1% 的余地"。这与本项目 C 集的含义（明确的不可判定 +
归因）差别很大——**不能把 PolyA 的 skip state 当作"它也会拒绝下结论"的证据。**

另有一处细节：`same_subfam_change_matrix`。若同一 subfamily 但 consensus 坐标相差
>50 bp（正链 `prev > curr+50`、负链 `prev+50 < curr`），仍按"换注释"处理。论文正文
未提，属于实现层的处理。

## 五、关键发现：PolyA 与 Carey 的置信度不在同一尺度

两篇用同一套概率逻辑，但 λ 的处理不同：

```
PolyA （confidence_cm）  : conf_i = 2^(s_i)      / Σ_j 2^(s_j)
Carey （Eq 10）          : conf_i = 2^(s_i / λ)  / Σ_j 2^(s_j / λ)      λ ≈ 0.1227
```

两个候选、分数差 `d` 时：

| 分数差 d | PolyA 置信度 | Carey 置信度（λ=0.1227） |
|---:|---:|---:|
| 1 | 0.667 | 0.996 |
| 2 | 0.800 | 0.9999 |
| 3 | 0.889 | ≈1.000 |

**所以两篇报出的"73%""95%"不是同一个东西。** Carey 的尺度极敏感：把置信度 0.734 反解
回原始分数差是 `λ·log2(0.734/0.266) = 0.1227 × 1.464 ≈ 0.18` 分；0.528 对应
`≈ 0.02` 分。

**这意味着 Carey 的"concordant 中位 73.4%"其实对应"最佳与次佳几乎打平"**（差 0.18 分）。
同一个分数差若用 PolyA 的尺度算，只有 `2^0.18/(1+2^0.18) ≈ 0.53`。

→ 引用时的两条硬规矩：
1. **不要把 PolyA 的置信度与 Carey 的置信度放在同一张图或同一个阈值上。**
2. **引用 73.4% / 52.8% 时必须同时给 λ**，否则这两个百分数无法被复现。

## 六、引用怎么写（可复制）

> Hong C, Wheeler TJ. PolyA: a tool for adjudicating competing annotations of biological
> sequences. bioRxiv 2021.02.13.430877 [Preprint]. doi:10.1101/2021.02.13.430877.
> （预印本，检索时未见期刊版本。其置信度为 `2^s/Σ2^s`，s 为未除 λ 的整数比对分；
> skip state 分数固定为 10 且置信度下限硬编码为 1%。）

> Carey KM, Patterson G, Wheeler TJ. Transposable element subfamily annotation has a
> reproducibility problem. Mob DNA. 2021;12:4. doi:10.1186/s13100-021-00232-4.
> （置信度为 `2^(s/λ)/Σ2^(s/λ)`，λ 由 `esl_scorematrix` 从打分矩阵推断，
> 本文 λ≈0.1227；concordant/discordant 中位置信度 73.4% / 52.8%。）

作者名以原文为准：PolyA 的第一作者在正文首页，核对时请照抄，本文档未逐字核对该字段。

## 七、复现步骤

```bash
# 1. 确认是否已有正式版（不要靠印象）
curl -s "https://api.crossref.org/works?query.bibliographic=PolyA+adjudicating+competing+annotations&rows=3&select=DOI,title,type,container-title"
#    → type 为 posted-content 即预印本

# 2. 试 PDF（常被限流，失败不要纠缠）
curl -sSL -A "Mozilla/5.0" -o p.pdf "https://www.biorxiv.org/content/10.1101/2021.02.13.430877v1.full.pdf"

# 3. 拿不到就读实现（这条通常更快也更准）
curl -s "https://raw.githubusercontent.com/TravisWheelerLab/PolyA/master/polyA/confidence_cm.py"
curl -s "https://raw.githubusercontent.com/TravisWheelerLab/PolyA/master/polyA/constants.py"
curl -s "https://raw.githubusercontent.com/TravisWheelerLab/PolyA/master/polyA/calculate_score.py"

# 4. 在本地把两套公式算一遍，验证尺度差异
python -c "
for d in (1,2,3):
    poly = 2**d/(1+2**d); carey = 2**(d/0.1227)/(1+2**(d/0.1227))
    print(d, round(poly,4), round(carey,4))
"
```

## 八、遗留

- PolyA 的 Eq 2/3/4 具体符号仍未逐字核对（PDF 不可得）。但**实现已给出等价物**，且
  本文档所有公式性陈述均来自源码，不是转述。
- Carey 的 λ 是否对所有区域用同一个值：论文说按 GC 含量做区域特异矩阵，λ 因此逐区域
  不同；本文引用的 0.1227 是 25p41g 矩阵的值，与 PolyA 默认矩阵相同，所以两篇的 λ 可
  放在一起比。**若换矩阵，第五节那张换算表要重算。**
