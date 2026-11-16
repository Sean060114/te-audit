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
| **0 名称规范化** | EDTA 的 `DNA/hAT` 和 TEsorter 的 `hAT_TIR_transposon` 到底是不是一回事？ | `te-normalizer/normalize.py` | 57 | `te-normalizer/README.md` |
| **1 完整性检查** | 这份注释的哪些结论站不住？ | `te_check.py`（C1–C6）<br>`bench_multi.py`（多数据集横向） | 74 | `te_check_README.md`<br>`TE_case/报表_*.md` |
| **2 交叉验证** | EDTA 与独立注释（Araport11）差在哪，谁更可能对？ | `compare_araport11.py`<br>`candidates_novel.py` | — | `EDTA_vs_Araport11_评估报告.md`<br>`EDTA_交叉验证报告.md` |
| **交付** | 把一份注释整理成可交付、可自校验的包 | `build_delivery_package.py`<br>`validate_delivery.py` | — | `delivery/`（不入库，可重建） |
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

判据在 **五份独立运行、两个科** 的数据上验证过（详见 `TE_case/`）：

| 数据集 | 物种 | intact 元件 | C4 已定超家族 LTR |
|---|---|---:|---|
| Chr1 | A. thaliana Col-0 | 218 | 40/40（100%） |
| col1m | A. thaliana Col-0 | 69 | 4/4（100%） |
| alyrata | A. lyrata | 20 | 3/3（100%） |
| rice_pre | O. sativa（单子叶） | 114 | 14/14（100%） |
| rice_full | O. sativa（单子叶） | 1,669 | 187/198（94.4%） |

**C1 跨科成立**（五份一致）；**C3 指标无效**（阈值投影，不可跨数据集比较）；
**C4 分层后成立**（此前"阈值被证伪"的结论是分母算错，已推翻）；
**C6 跨科复现**（Chr1 / alyrata / 水稻各若干例，成因需逐个区分）。

报告里保留了完整的自我修正记录，未删改——C3 的表述修正、C4 的分母修正、
`.mod.mod` 归因的两次撤回，都原样留在文档里。错误的分母和它导致的错误结论，
本身就是这套方法最值得记下的部分。

---

## 三、复现入口

```bash
# 单份数据：六项检查（零依赖）
python te_check.py --edta-dir rice/edta_out --json rice/full_report.json

# 多数据集横向对比 + 判据检验
python bench_multi.py \
  "chr1=chr1_edta_out" "col1m=bench_archive/col1m.fa.mod" \
  "alyrata=bench_archive/alyrata.fa.mod" "rice_full=rice/edta_out" \
  --tsv TE_case/bench_five.tsv

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
  报表_多基因组判据检验.md     拟南芥属内三份数据（⚠ 顶部注明 C4 结论已推翻）
  报表_水稻跨科判据检验.md     水稻跨科 + Chr1 重跑验证（现行结论以此为准）
  报表_结构注释质量基线.md     最早的单点基线
archive/                已被取代的运行脚本（保留是为了留下失败记录）
bench_archive/          col1m / alyrata 的历史运行回传（仅此一份）
chr1_edta_out/          TAIR10 Chr1 基准运行回传
rice/                   水稻运行回传 + 报告 JSON（序列文件不入库）
ref/                    参考序列（不入库，可从 Ensembl Plants release-57 重下）
delivery/               交付包（不入库，用 make_deliverables.sh 重建）
papers/                 立项时读的文献
选题说明书_详细论证版.md   立项论证，含最初的判据设想
```

`.gitignore` 里写明了哪些内容故意不入库、以及判断标准（重跑一次要多久）。

---

## 五、产出清单与处置

根目录 40 个条目，按「还需不需要」分四类。当前全部体积 187 MB。

### 保留：核心代码（改了就坏）

| 条目 | 作用 |
|---|---|
| `te_check.py` `test_te_check.py` `bench_multi.py` | 模块 1 全部逻辑与 74 条测试 |
| `run_bench_rice.sh` | 现行 EDTA 运行脚本（参数化，含命名修复） |
| `te-normalizer/` | 模块 0，独立可跑 |
| `compare_araport11.py` `candidates_novel.py` `make_publication_tables.py` | 模块 2 |
| `build_delivery_package.py` `validate_delivery.py` `make_deliverables.sh` `make_union_softmask.py` `get_dfam_lib.sh` `run_validation.sh` | 交付包与验证 |
| `archive/` | 三版历史运行脚本，失败记录 |

### 保留：分析证据（报告里每个数字的依据）

| 条目 | 体积 | 为什么不能删 |
|---|---:|---|
| `chr1_edta_out/` | 1.7 M | TAIR10 Chr1 基准，重跑需 47 分钟 |
| `rice/` | 66 M<br>（入库 11 M） | 水稻全量回传，重跑需 73 分钟；序列文件不入库 |
| `bench_archive/` | 251 K | col1m / alyrata 历史运行，**仅此一份** |
| `TE_case/` | 30 M<br>（入库 <1 M） | 报告、基线、案例；其中 `.venv` 是 classify_element.py 的环境，文档里写了重建命令 |
| `araport11_vs_edta.txt` `.per_class.tsv` | 12 K | 模块 2 报告的逐类数字来源 |
| `Araport11_TE_family_summary.tsv` `EDTA_B_level_candidates.tsv` `EDTA_candidate_novel.tsv` | 144 K | 模块 2 的三张成果表 |
| `chr1_report.json` | 8 K | Chr1 六项检查结果 |

### 保留：交付物与立项文档

| 条目 | 说明 |
|---|---|
| `delivery/`（35 M）| 交付包解包版，12 文件 + MD5SUMS，可用 `make_deliverables.sh` 重建 |
| `TAIR10.Chr1.TE_annotation_v1.tar.gz`（9.8 M）| 同一交付包的分发版 |
| `EDTA_vs_Araport11_评估报告.md` `EDTA_交叉验证报告.md` `EDTA_LTR丢失溯源报告.md` | 模块 2 报告 |
| `选题说明书_详细论证版.md` `附录B_anno阶段实测.md` `附录_实战记录.md` | 立项与过程记录 |
| `papers/` `papers_digest.txt` | 立项文献与提要 |
| `EDTA_issue_台账.csv` | 上游 issue 与 PR 记录 |

### 不入库的外部数据（体积大，可重新下载）

| 条目 | 体积 | 怎么回来 |
|---|---:|---|
| `ref/TAIR10.Chr1.fa` | 30 M | `curl` Ensembl Plants release-57（命令见第 3 节） |
| `araport11_TE.gff3` | 9.2 M | Araport11 官方发布 |
| `rice/*.fa` `*.fa.gz` | 55 M | 同上，Ensembl release-57 |

### 已删除

| 条目 | 原体积 | 为什么可以删 |
|---|---:|---|
| `chr1.fa.mod.union.softmasked.fa` | 30 M | EDTA ∪ Araport11 并集软屏蔽的中间产物，交付包里的 `TAIR10.Chr1.genome.softmasked.fa` 已取代它。重建方式见下 |
| `__pycache__/` | 16 K | Python 缓存，自动重建 |

重建那个软屏蔽文件（输入需先改名对齐）：

```bash
cp chr1_edta_out/Chr1.fa.mod.EDTA.TEanno.gff3 chr1.fa.mod.EDTA.TEanno.gff3
python make_union_softmask.py    # 输出 chr1.fa.mod.union.softmasked.fa
```

## 六、已知空白

1. ~~**`.mod.mod` 双前缀的成因未定。**~~ **已查明。** EDTA.pl 已把基因组归一化成
   `X.fa.mod`；v2.2.2 的 `EDTA_raw.pl:440` 调用 LTR_retriever 时未传
   `-convert_seq_name`（默认 1），后端又归一化一次，`encode_fasta` 的输出路径是
   `"$genome.mod"`，于是得到 `X.fa.mod.mod`。**与序列名长短无关**，
   v2.3.0 起补上该参数后不再出现。详见 `TE_case/报表_水稻跨科判据检验.md` 5.3 节。
   当前仍用改名修复而不升级——五份数据都在 v2.2.2 下跑，升级会破坏可比性。
2. **水稻全量有 11 个已定 Gypsy 未进分类表。** 结构完整、超家族已定。
   原先怀疑与 1 同源，但归因查实后该推论已下调（短 ID 的归一化不改变序列内容），
   目前无解释。
3. **Chr1 的输入来源与旧运行无法逐字节核对。** 旧值取自 EDTA 自带 test 集，
   本次用 Ensembl release-57；两个计数（218 / 41）吻合，但原始文件已失。
4. **模块 2 只做过 Chr1 一份**，未跨到水稻，也尚无独立真值集（RITE / RAP-DB）。

## 七、运行环境

| 组件 | 版本 | 说明 |
|---|---|---|
| Python | 3.13 | 全部代码零第三方依赖 |
| EDTA | 2.2.2（conda env `edta`，WSL） | 全部数据均此版本，升级会让新旧数据不可比 |
| 参考组装 | IRGSP-1.0 / TAIR10 | 均取自 Ensembl Plants release-57 |
