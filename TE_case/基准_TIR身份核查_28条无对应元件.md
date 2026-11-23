# 28 条"MSU 无对应"TIR 元件的身份核查

**要查的问题**：上一节留下的 28 条元件，用 Repbase/DFam 的水稻 Mutator、CACTA consensus
比对，看它们是真转座子还是拼接出来的假结构。

**结论先行**：

1. **Repbase 拿不到**（GIRI 商业授权，不能下载、不能随仓库分发）。
2. **Dfam 里没有水稻**——不是"没搜到"，是库里一条都没有（§1 有可复现的取证）。
3. 换了一条不依赖外部共识库的通路：**MSU 自己的基因注释流程**。它给出了 3 条的名字，
   其中 1 条明确是 **Mutator sub-class**。
4. 剩下 25 条里，5 条完整裹住宿主基因（两种读法分不开），20 条无判据可用。

**所以这 28 条不同质，不能整体判真或判假。** 能写进报告的只有 §4 的那 3 条 +
§5 的 MSU 内部矛盾。

---

## 一、Repbase / Dfam：两条路都不通，而且可以证明

### 1.1 Repbase

GIRI 的商业数据库，需订阅。不能下载，也不能随仓库分发——即使拿到了，把 consensus
序列提交进仓库同样是许可问题。**这条路在本项目的条件下不存在。**

### 1.2 Dfam：不是没找到，是库里没有

Dfam 是 Repbase 的公开继任者，本来是最合理的替代。实测 Dfam 3.9（30,646 个家族）
**没有任何一个水稻家族**。取证用的是公开 API，逐条可复现：

```
$ curl -s "https://dfam.org/api/families?name=Oryza&limit=200"      → total_count 0
$ curl -s "https://dfam.org/api/families?name=Dasheng&limit=200"    → 0   （水稻 CACTA 代表家族）
$ curl -s "https://dfam.org/api/families?name=RIRE&limit=200"       → 0   （水稻 LTR）
$ curl -s "https://dfam.org/api/families?name=Tos17&limit=200"      → 0   （水稻 MITE）
$ curl -s "https://dfam.org/api/families?name=Mutator&limit=200"    → 0
$ curl -s "https://dfam.org/api/families?name=CACTA&limit=200"      → 0
$ curl -s "https://dfam.org/api/families?name=Tourist&limit=200"    → 0
$ curl -s "https://dfam.org/api/families?name=Stowaway&limit=200"   → 0
```

看似命中的两条，查下去都不是植物：

| 检索词 | 命中 | 实际是什么 |
|---|---|---|
| `name=Ory` | 85 | 全是 `OryAfe-1.*` —— **Orycteropus afer（土豚）**，Zoonomia 哺乳动物项目的产物 |
| `name=MULE` | 43 | `MULE_B.platifrons_25`（双壳类）、`MULE_S.constricta_33`（海胆）、`MULE_T.granosa_11`（牡蛎）—— 全是海生无脊椎动物 |
| `name=EnSpm` | 48 | `EnSpm-N16_DR`（斑马鱼）、`CMC-EnSpm-13_AstCal` |

`clade=` 参数不可靠，不能用来佐证：`clade=human` / `clade=Oryza_sativa` / `clade=plant`
都原样返回全部 30,646 条（即参数被忽略），而 `clade=Viridiplantae` 返回 0。**唯一可信的
是 name 检索，而 name 检索给出的答案是零。**

> 这一条同时给仓库里 `get_dfam_lib.sh` 的那句"下载的库经实测不含任何植物序列"补上了
> 凭据。那个脚本当时被标为废弃，理由是经验性的；现在可以说清楚原因：**Dfam 的植物
> 转座子收藏本来就只有零星几个类群，水稻整体缺席。**

## 二、替代判据一：家族性聚类（结果弱，但要说）

如果这 28 条是真的 TIR 家族成员，同族的其它拷贝应该也在同一份基因组里、且大概率也
被 EDTA 检出。做法：把 1,824 条结构 TIR 的左臂两两比对（8-mer 预筛 + −8..+8 位移的
无缺口比对，相似度 ≥0.70、臂长比 ≥0.75 记同族，并查集），脚本
`TE_case/tir_family_cluster.py`，结果 `TE_case/tir_family.tsv`。

聚类本身是可信的：1,824 条聚成 1,367 簇，最大簇 314 条、臂长 20–22 bp、全部 DTT
——就是 Tc1/Mariner 那个短 TIR 大家族，没有过度合并。

| 组 | n | 家族大小 中位 | 孤立（=1）占比 |
|---|---|---|---|
| 一致（agree_superfamily + agree_class） | 27 | 2 | 36% |
| 冲突 | 7 | 1 | 57% |
| **无对应（28 条）** | 28 | **1** | **61%** |

只看臂长 ≥100 的子集（与"长臂"质疑直接相关）：无对应 61% 孤立，一致组 47% 孤立。

**方向对，区分度不够。** 61% vs 47% 在 n=28/19 上说明不了什么，而且"孤立"本身有两种
解释（真·单拷贝古老元件 / 拼出来的结构），这个判据分不开。**不作定性依据。**

## 三、替代判据二：重复含量零模型（结果也是弱的，且有一个循环）

"这 4 kb 在 MSU 里几乎没有重复注释"听起来很可疑。但要看零模型：随机取 8,000 个同长
窗口，MSU TE 覆盖率的**中位是 0.168，P25 是 0.000**——四分之一的随机窗口本来就是零重复。
水稻的重复分布是双峰的，4 kb 窗口很容易落在基因间区。

| | MSU TE 覆盖率 中位 | 低于该值的随机窗口占比 |
|---|---|---|
| 随机同长窗口（n=8,000） | 0.168 | — |
| 无对应 28 条 | 0.012 | 27% |
| 一致 27 条 | 0.649 | 79% |

28 条落在最低 27% 里，**不算异常**。

更要紧的是：一致组那个 0.649 是**循环的**——"MSU 有对应"这个分组本身就是用 MSU 的
重复注释定义的，再拿重复覆盖率去比较它，等于把定义当结论。**这条也不作证据。**

## 四、替代判据三：MSU 基因模型（这条给出了名字）

关键转折：MSU 的 `repeat_masked.gff3`（重复流程）和 `all_models.gff3`（基因流程）是
**两条互不通气的流水线**。重复流程说"这 4 kb 没有 TE"，不等于 MSU 对这块地方一无所知。

下载 `osa1_r7.all_models.gff3.gz`（rice.uga.edu，8.9 MB）后，逐个元件统计"两端都落在
元件内"的完整基因模型。脚本 `TE_case/gene_occupancy.py`。

零模型：8,000 个同长随机窗口里，**11%** 含至少一个完整基因模型，外显子占比中位 0.160。

| 组 | 去重后 n | 含 ≥1 个完整基因 | vs 随机 | 含完整**宿主**基因 | 含完整**转座子蛋白**基因 |
|---|---|---|---|---|---|
| 一致组 | 25 | 16% | +5 pp | 12% | 4% |
| **无对应** | **25** | **32%** | **+21 pp** | **20%** | **12%** |
| 冲突 | 7 | 14% | +3 pp | 0% | 14% |

无对应组含完整基因的比例是随机窗口的 3 倍、一致组的 2 倍。以随机窗口的 11% 为零假设：
无对应组 8/25，二项检验 **p = 0.004**（期望 2.75 个）；作为对照，一致组 4/25 的
p = 0.29，**与随机无差别**。也就是说"含完整基因"这件事只在无对应组里显著。

**三条拿到了名字：**

| 位点 | 长度 | 两臂 | TSD | MSU 基因模型 | 注释原文 |
|---|---|---|---|---|---|
| OsChr12:2628957 | 5,021 | 251/252 bp | 9 bp | LOC_Os12g05700，6 外显子，CDS 2,841 bp | transposon protein, putative, **Mutator sub-class**, expressed |
| OsChr1:36609926 | 5,016 | 156/156 bp | 7 bp | LOC_Os01g63170，6 外显子，CDS 2,892 bp | transposon protein, putative, unclassified, expressed |
| OsChr12:1202488 | 4,066 | 222/222 bp | 7 bp | LOC_Os12g03160，5 外显子，CDS 2,061 bp | transposon protein, putative, unclassified, expressed |

这三条**同时满足 MULE 的架构特征**，而这一点完全是独立于基因注释的：

- 臂长 156–252 bp —— MULE 的 TIR 典型 150–250 bp，而本项目里 Tc1/Mariner 那个
  314 条的大簇臂长只有 20–22 bp；臂长本身就足以把这三条排到 Mutator 一侧；
- TSD 7–9 bp —— Mutator 的 TSD 通常 8–10 bp，其中命名到 Mutator 的那条正好 9 bp；
- 内部 2.0–2.9 kb 的单一转座酶 CDS —— 自主元件的标准配置。

**所以至少这 3 条是真转座子，而且能被命名到 Mutator（MULE）一级。** 用户要的
"Mutator consensus 比对"没有做成，但得到的结论比"像不像某个 consensus"更硬：
MSU 的基因流程在同一个坐标上放了一个完整的、表达了的、6 外显子的转座酶基因。

## 五、由此暴露的 MSU 内部矛盾

上面这 3 条，在 MSU 自己的两份文件里说法相反：

```
repeat_masked.gff3  → 这 4–5 kb 上只有零散小 MITE/SINE（覆盖率 0.9%–3.4%）
all_models.gff3     → 同一坐标上有一个完整、表达、剪接的 Mutator 转座酶基因
```

两条流水线不互相校验。这与上一节在 EDTA 侧发现的矛盾是同一类问题
（EDTA 的 TIR 流程与自己 LTR 流程打架，MSU 站 LTR），但方向相反：
**这次是 MSU 的基因流程替 EDTA 说了话。**

对基准的含义：**"MSU 没有注释"不能当作"这里没有 TE"的证据。** 上一节把 28 条记成
"MSU 不支持整段"是对的，但不能进一步推成"EDTA 报错了"。

## 六、剩下 22 条

**5 条完整裹住宿主基因**（`expressed protein`，750–3,763 bp，含一条 8 外显子 3.7 kb 的
LOC_Os01g50800）：Chr1:4586706、Chr1:29178746、Chr1:37694047、Chr1:39320381、
Chr1:39321160。

这 5 条有两种读法，**用现有数据分不开**：

- **A. Pack-MULE**：水稻 Mutator 家族以捕获宿主基因著称（Pack-MULE），一个 4–5 kb 的
  MULE 内部装一个完整基因是**该家族的正常形态**。若如此，这 5 条是真的。
- **B. 拼接假阳性**：TIR-Learner 把两个相距 4 kb 的反向重复接起来，中间正好裹进一个基因。
  若如此，这 5 条是假的。

要分开，需要看那个基因在元件外有没有旁系拷贝。**这一步已经做了，结论是判据不成立**
（见 §六之二）。

### 六之二、Pack-MULE 判别：判据本身不成立

判别式本来很干净：真 Pack-MULE 意味着元件里装的是**捕获副本**，那么 donor 应该还在
基因组里；如果那个基因全基因组只有一份，"捕获"无从谈起，只能是元件边界把唯一的原
基因裹了进去。

做法（`TE_case/pack_mule_test.py`）：取每个外显子的 20 bp 种子在全基因组精确枚举，
±6 位移内无缺口比对、相似度 ≥0.85 记同源，5 kb 内聚成位点。**零模型是必须的**——
水稻本身就多片段重复，随便挑个基因也常有旁系拷贝，所以随机抽 10 个同规模基因跑同一
流程做基线。

| | 旁系位点数 | 含"近同源副本"（过半外显子、平均 ≥0.95） |
|---|---|---|
| 随机 10 个基因（基线） | 中位 1，**50% 为 0** | **20%** |

| 元件 | 基因 | 外显子 | 旁系位点 | 最佳位点 |
|---|---|---|---|---|
| OsChr1:4586706 | LOC_Os01g09130 | 1 | **0** | — |
| OsChr1:29178746 | LOC_Os01g50800 | 8 | **14** | 6/8 @ 0.988 |
| OsChr1:37694047 | LOC_Os01g64949 | 2 | 2 | 1/2 @ 1.000 |
| OsChr1:39320381 | LOC_Os01g67660 | 3 | **0** | — |
| OsChr1:39321160 | LOC_Os01g67660 | 3 | **0** | — |

**判据不成立的原因是基线**：随机基因里就有 **50%** 是单拷贝。所以"这个基因在元件外
没有拷贝"丝毫不反常，"唯一 → 边界画错"的推理链断在第一步。这条判据**没有区分力，
5 条一条都没判掉。**

唯一有信息量的是 LOC_Os01g50800：14 个旁系位点，其中一处 6/8 外显子平均相似度
0.988。这与"近期捕获副本"的形态吻合（基线里只有 20% 的基因有这种近同源副本），
但反过来也与"普通的多拷贝基因家族"吻合，单向读不出来。

三条辅助事实都指向"元件裹住的是一个普通基因"，但都只是倾向、不是判决：

1. 这 5 个基因都**含内含子且有完整 UTR**（如 LOC_Os01g50800：3,763 bp 跨度为
   2,037 bp CDS），不是加工过的短副本；
2. MSU 的功能注释把它们写作 `expressed protein`（普通基因），而不是那 3 条转座酶
   基因的 `transposon protein`（TE 相关模型）；
3. 5 条所在的元件长度 1.4–4.6 kb，基因几乎占满整个内部，两侧留给 TIR 之外的
   侧翼序列很少。

**结论：5 条留在 B 集，不因这条判据升降。** 要真正分开，只能做共线性：拿 donor 的
侧翼序列比对元件两侧的侧翼，看是不是同一段基因组上下文。

**其余 17 条**：既没有完整的转座子蛋白基因，也没有完整的宿主基因，MSU 侧只有零散小片段。
**无判据可用**，只能留在"待定"。

## 七、对 A/B/C 的影响

上一节把 28 条挡在 A 集外面。现在要拆开：

| | 条数（去重） | 归入 | 理由 |
|---|---|---|---|
| 有完整转座酶基因 | 3 | **A** | 独立流程 + 架构双重支持 |
| 完整裹住宿主基因 | 5 | **B** | Pack-MULE 与拼接口分不开；旁系检索的**判据本身不成立**（基线 50% 基因就是单拷贝） |
| 无任何判据 | 17 | **C / 待定** | 只有结构证据 |

TIR 的 A 集因此从 27 条变成 **30 条**（27 + 3），B 集 +5。这比"62 条全靠结构证据"
或"28 条一律可疑"都更接近实情。

## 八、还能做什么（按性价比排序）

1. **Repbase**：有机构订阅的话，直接拿 `MULE1_OS`、`CATCA_OS`、`Dasheng` 这些水稻条目
   做 RepeatMasker。这是唯一能给出"家族名 + consensus 相似度"的路径。
2. ~~**Pack-MULE 判别**（paralog 检索）~~ → **已做，判据不成立**（§六之二）。随机基因
   里 50% 就是单拷贝，"唯一基因 → 边界画错"的推理断在第一步。**别再重跑这一条。**
   要分开只剩共线性一条路：拿 donor 的侧翼序列比对元件两侧侧翼，看是否为同一段基因组
   上下文——需要比对工具，本机没有。
3. **转座酶系统发育**：3 条阳性元件的 CDS 与已知 MULE 转座酶（MURA）建树，确认是否
   落在 Mutator 支。需要比对工具，本机没有（无 blastn/HMMER），只能在 WSL 的 edta
   环境里做。
4. **把 MSU 基因流程当独立证据用起来**：`all_models.gff3` 里 1,364/6,529 个 Chr1 基因
   模型是 `transposon protein`，这是一份被本项目忽略的、与重复注释独立的 TE 证据源。
   建议后续对全部 420 条 LTR 也跑一遍同样的比对。

## 九、复现

```bash
# 1. Dfam 取证（无需下载）
curl -s "https://dfam.org/api/families?name=Mutator&limit=200"

# 2. 家族性聚类
python TE_case/tir_family_cluster.py \
  --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
  --ledger TE_case/tir_msu_ledger.tsv --tsv TE_case/tir_family.tsv

# 3. Pack-MULE 判别（paralog 检索，约 2 分钟）
python TE_case/pack_mule_test.py \
  --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa --controls 10

# 4. 基因模型占用（需先下载 MSU 基因模型，8.9 MB）
curl -O https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.gff3.gz
python TE_case/gene_occupancy.py --gff osa1_r7.all_models.gff3.gz \
  --ledger TE_case/tir_msu_ledger.tsv --dedup --draws 4000
```

`osa1_r7.all_models.gff3.gz` 与 `ref/` 下的序列一样**不入库**（8.9 MB，rice.uga.edu
可直接下）。`gene_occupancy.py` 对截断的 gz 会警告而不是崩，便于边下边跑。
