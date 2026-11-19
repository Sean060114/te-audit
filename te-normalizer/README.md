# te-normalizer — TE 注释质量事实检查器（模块 0 + 模块 2）

两个脚本，共用同一套"名字 → 类别/角色/超家族"判断，纯 Python 标准库，零第三方依赖。

| 脚本 | 模块 | 回答什么问题 |
|---|---|---|
| `normalize.py` | 模块 0 | 这个名字该叫什么？两个名字是不是同一件事？ |
| `check_counts.py` | 模块 2 | 这份注释"有几个 TE""占多少"？各家用的是哪种数法？ |

`normalize.py` 输出：SO 标准名 + SO ID（命中时）、格式合法性判断、
类别/角色/超家族三级信息、未知名的 top-3 相似建议。
另有 `--pairs` 模式：输入两列名字，判定"等价 / 同目不同科 / 跨类冲突 / 角色不同"。

## 目录结构

```
te-normalizer/
├── normalize.py                 # 模块 0 核心（七个函数，全中文注释）
├── check_counts.py              # 模块 2 核心（口径核查）
├── tests.py                     # 33 条断言（标准库 unittest，python tests.py 直接跑）
├── make_alias_pr.py             # 从 #542 缺失清单生成 PR 素材（别名补丁 + 分诊表 + PR 正文）
├── data/
│   ├── TE_Sequence_Ontology.txt # EDTA master 的 SO 映射表（修复分隔符解析后 809 个别名 / 105 个 SO 名）
│   ├── aliases_issue542.csv     # 从 #542 全帖程序化提取的 73 个名称比对结果
│   ├── issue542_comments.json   # #542 全部 15 条评论存档
├── pr_issue542/                 # #542 别名 PR 素材（make_alias_pr.py 的产物）
│   ├── aliases_triage.tsv       # 61 条分诊表：目标 SO / 决策 / 分组 / 理由 / 风险
│   ├── aliases_patch_main.diff  # 21 个（主补丁，建议提交）
│   ├── aliases_patch_main_plus.diff  # 25 个（主补丁 + 4 条判据型）
│   ├── aliases_patch_clade.diff # 30 个（clade 级，备选，独立应用）
│   ├── TE_Sequence_Ontology.patched.txt     # 主补丁+判据型的结果文件
│   ├── TE_Sequence_Ontology.with_clade.txt  # 全量 55 条的结果文件
│   ├── verify_before.csv        # 未打补丁时的工具实测输出（基线 46/61）
│   ├── verify_main.csv / verify_main_plus.csv / verify_clade.csv
│   │                            # 三档补丁各自的实测输出（51 / 55 / 46）
│   ├── PR_BODY.md               # 可直接粘贴的英文 PR 描述
│   └── SUBMIT.md                # 提交操作清单（fork→分支→应用→自查→push→开 PR）
└── examples/
    ├── test_names.txt           # 测试集（覆盖四条路径）
    ├── test_names_counts.txt    # 计数列两种排法的混合测试
    ├── edta_chr1_TElib_names.txt# 真实 EDTA 库里的 14 个分类名
    ├── edta_chr1_report.csv     # 上者的映射报告（v1 起含 类别/角色/超家族）
    ├── pairs_matrix.txt         # 真实分歧对：从 EDTA 混淆矩阵提取的 17 条（8 对）
    ├── pairs_verdicts.txt       # 四种判定的回归用例
    ├── class_table.csv          # 模块 2 逐类对照表（真实数据产物）
    ├── counts_demo.out          # 模块 2 完整运行输出（可当报告截图用）
    ├── verbose_demo.out         # test_names.txt 跑 --verbose 的输出快照
    ├── weight_demo.out          # pairs_matrix.txt 跑 --pairs --weight 的输出快照
    ├── rate_demo.out            # 再给 --total 1593 --baseline 1.07：17/1593 = 1.07% 与 EDTA 对齐
    └── report.csv               # 早期运行产物示例

    # 约定：examples/ 下的 .txt 一律是「输入」（可直接喂给 normalize.py），
    #       .out / .csv 一律是「产物」（输出快照或报告，不要当输入）。
    #       此前 verbose_demo / counts_demo / weight_demo 三个产物用了 .txt 后缀，
    #       误当输入时会把报告文字当名字读，全部判成"未命中"而不报错。
```

## 运行

```bash
# 跑测试（改代码后第一件事）
python tests.py
python tests.py -v                              # 逐条用例名

# 模块 0
python normalize.py examples/test_names.txt
python normalize.py examples/test_names.txt --verbose          # 逐条：原始名 -> 命中键 -> SO 名
python normalize.py examples/test_names.txt --out report.csv
python normalize.py examples/pairs_matrix.txt --pairs          # 分歧对判定
python normalize.py examples/pairs_matrix.txt --pairs --verbose # 配对 + 归一明细
python normalize.py examples/pairs_matrix.txt --pairs --out pairs_report.csv
python normalize.py examples/pairs_matrix.txt --pairs --weight   # 按次数加权
python normalize.py examples/pairs_matrix.txt --pairs --weight \
       --total 1593 --baseline 1.07    # 再加：与 EDTA 同分母的加权不一致率

# 模块 2
python check_counts.py --gff3 chr1.fa.mod.EDTA.TEanno.gff3 \
                       --sum  chr1.fa.mod.EDTA.TEanno.sum \
                       --bed  chr1.fa.mod.EDTA.TEanno.bed \
                       --out  class_table.csv
```

## 四条查找路径（normalize.py 的核心逻辑）

1. 精确命中别名表 → SO 标准名
2. 大小写不敏感命中 → SO 标准名（并标注真实写法）
3. 三级命名 → 剥一层再查。**两种层级顺序都试**：`LTR/Copia/SIRE` 砍最后一级，
   `ClassI/LTR/Gypsy` 取最后两级
4. 连字符/下划线互换再查
5. 全部失败 → 编辑距离 top-3 建议（绝不静默丢弃）

## `--verbose`：为什么必须显示"命中键"

默认只给一行 `[matched_ci] sine/alu -> Alu_SINE_retrotransposon`。这句话有个陷阱：
看起来像"表里有 sine/alu"，实际表里写的是 `SINE/Alu`，只差大小写。

`--verbose` 把中间那一层拆出来，输出三段：

```
路径          原始名                    命中键                    SO 标准名 (SO ID)
路径2 大小写  sine/alu                  SINE/Alu                  Alu_SINE_retrotransposon (SO:0000206)
路径3 剥级    LTR/Copia/SIRE            LTR/Copia                 Copia_LTR_retrotransposon (SO:0002264)
路径4 连字符  DNA/CMC_Transib           DNA/CMC-Transib           Transib_TIR_transposon (SO:0002282)
路径5 未命中  nMITE                     —                         建议：MITE(0.89); TIR/MITE(0.62); MITE/hAT(0.62)
```

**"命中键"这一列是判断"要不要修上游库"的唯一依据**：如果大量输入都靠路径 2 命中，
说明库里的写法与用户习惯不一致，应该提 PR 加别名，而不是让每个用户自己改输入。

两处实现细节：`pad()` / `clip()` 按**显示宽度**而不是字符数补齐（中文占 2 列，
用 `len()` 会把表格排歪）；CSV 报告也同步加了"命中键"列，方便事后审计。

## 为什么要"类别 / 角色 / 超家族"（2026-11-15 实测依据）

EDTA 的 `TEanno.gff3` 里，同一个 `classification=LTR/Copia` 同时挂着 4 个 SO ID：

| SO ID | 是什么 | 角色 |
|---|---|---|
| SO:0002264 | Copia 元件本体 | element |
| SO:0000286 | 它的长末端重复 | sub_feature |
| SO:0000434 | 它的靶位点重复 | sub_feature |
| SO:0000657 | 装它的容器 | region |

不看"角色"，就没法把 classification 映射到唯一的 SO ID。同理，注释行 3666 行
vs 元素级 3461 行，差的 205 个正是 `sub_feature` / `region` 行——不会区分角色，
就会把零件当整机数。所以 `describe_so()` 补的这三列不是装饰，是后续模块计数正确的前提。

## `--pairs` 的五种判定

| 判定 | 条件 | 含义 |
|---|---|---|
| 等价 | 归一后同一 SO ID | 同概念不同写法，可直接合并 |
| 同科异名 | 同目、同超家族、不同 SO ID | 同一批家族的两种 SO 写法 |
| 同目不同科 | 同目、不同超家族 | **EDTA 只记录、不处置的那一类** |
| 跨类冲突 | Class I vs Class II | 必有一方判错 |
| 角色不同 | element vs sub_feature/region | 不可比，先按角色过滤 |

对今天实测的 17 条真实分歧：**6 条属"跨类冲突"（EDTA 的 order 级过滤器会丢弃），
11 条属"同目不同科"（记录后无任何处置，占 65%）**。这个 65% 就是模块 1 要接手的部分。

## 模块 2：`check_counts.py` 查的是什么

同一份 EDTA 注释（拟南芥 Chr1，30,263,312 bp），只换数法，答案就不同：

| 问题 | 数法 | 结果 |
|---|---|---|
| 有几个 TE | gff3 全部行 | **3666** |
| | 元素级条目（扣掉 205 行非元件） | **3461** |
| | RepeatMasker 汇总 | **3453** |
| 占多少 | bed 逐行累加（未合并） | 2,837,354 bp / **9.38%** ← 假数 |
| | bed 合并去重叠 | 2,796,034 bp / **9.24%** |
| | RepeatMasker 汇总 | 2,797,686 bp / **9.24%** |

**两个结论，一个是真分歧，一个是假分歧：**

1. **计数是真分歧**（3666 / 3461 / 3453）。205 行差在零件：`repeat_region` 41 +
   `target_site_duplication` 82 + `long_terminal_repeat` 82；剩下 8 个差是 RepeatMasker
   合并了相邻碎片。→ 报数时必须写清数的是整机还是含零件。
2. **覆盖是假分歧**。9.38% 是"没合并就累加"造成的：那 3461 行 bed 里有 236 处区间
   相交或首尾相接（嵌套元件 + 相邻碎片），重叠段被算两遍，虚增 41,320 bp。
   合并后 9.24%，与 RepeatMasker 只差 1,652 bp（0.005 个百分点，源于坐标约定）。

这条自身修正值得单独记住：**"两个官方数字不一致"这句话，必须先排除自己算错**。
第一版脚本没合并区间，得出的 9.38% 看起来像个重大发现，其实是我自己的 bug。

逐类对照还会暴露工具链内部的真实分歧（`class_table.csv`）：
`TIR/PIF_Harbinger` 汇总 311 vs gff3 320（+9）、`nonTIR/helitron` 1143 vs 1136（−7）
——这两项才是需要追的，因为它们不是零件、不是合并，而是分类归属不一致。

## 已有数据结论（2026-11-13 生成）

- #542 全帖共提取 73 个被报告名称，其中 12 个作者已补进 master，61 个仍缺失；
- 缺失清单 = `aliases_issue542.csv` 中状态为"缺失（候选 PR 内容）"的行，
  这就是未来向 EDTA 提"数据 PR"的候选素材（PR 前需逐个人工复核，
  参考 #542 中作者对 ARTEFACT 的处理：语义不明则维持 repeat_fragment 兜底）。

## 变更记录

- 2026-11-19（v1.5）修镜像组合并的两处判定缺陷（`rank_pairs`）。
  1. **身份分命名空间。** 旧实现用 `a_id or a` 当身份：未命中时退回原始名，而 SO ID
     本身也是字符串，两套东西共用一个命名空间。输入里出现 `SO:0002278` 这种名字时，
     它会与真正的 SO ID 撞成同一个键，把"等价"（两条都归到 SO:0002278）与
     "无法判定"（首列未命中）并成一组。现在退回时加 `name:` 前缀，两者永不相等。
  2. **组内判定不再取首条。** 旧实现 `setdefault(key, {"verdict": r["verdict"]})`
     只保留**最先进入那一行**的判定，同一份数据换行序就换标签，且毫无痕迹。
     现在收集全部判定：只有一种时用它，多于一种时标为 `混合`，并在行尾列出
     `判定混合：等价×1 / 无法判定×2`。**不替数据选边。**
  3. 这两条互为保险：命名空间分开后撞键几乎不可能，但"混合"检测留着——
     将来 `relate` 若引入带方向的判定（"谁更可信"），它会立刻把问题暴露出来，
     而不是悄悄取错。
  4. 干净数据上的输出逐字节不变（`examples/weight_demo.out` 未改动）。
  5. `tests.py` 69 → 73 条断言（撞键与混合判定各 2 条）。
- 2026-11-17（v1.4）`--total` / `--baseline`：把不一致率放回 EDTA 的分母上。
  1. 背景：`--weight` 报的是「分歧对占多少对 / 涉及多少条」，分母是**出现分歧的那部分**；
     EDTA 打印的 `Inconsistent_rate` 分母是**库里全部序列**（本次 17 / 1593 = 1.07%）。
     两个数碰巧都在 1% 附近，一旦并排引用就会得出"我们和 EDTA 一致"的错误结论。
  2. 新增 `print_rate()`：给了 `--total N` 才算，输出 非等价条数 / 加权不一致率；
     `--baseline X` 再与 EDTA 自报值并排，差值 < 0.005 pp 时写"一致"而不写 `-0.00`。
  3. **分子分两档**：`全部非等价`（可比 EDTA 口径）与`只算真冲突`（剔除"角色不同"——
     `element` vs `sub_feature` 不是判断分歧）。两个都报，读者自己选，不替他决定哪个算错。
  4. 分母给错（`--total` 小于实际涉及的序列数）时**不输出比率**，只告警——
     宁可少一个数，也不给一个分母错了的百分数。
  5. `--total` / `--baseline` 不带 `--pairs`、`--baseline` 不带 `--total`、`--total <= 0`
     三种误用均报错退出（exit 2），不静默忽略。
  6. `tests.py` 57 → 69 条断言（新增 12 条：`print_rate` 单元 6 条 + CLI 5 条 + 用例数守卫 1 条）。
     **守卫**：`TestSuiteGuard` 断言用例总数不少于 69，防止以后改代码时顺手删断言。
- 2026-11-16 `load_te_so` 增加别名冲突记账；修文件句柄泄漏。
  1. 返回值改为三元组 `(alias_to_so, so_names, conflicts)`。旧实现遇到「同一别名指向不同
     SO 名」时后者静默覆盖前者，查表结果取决于文件里的行走顺序，出了问题毫无痕迹。
     现在冲突被收集起来，CLI 启动时打到 stderr（含最终生效值与前 10 个冲突项）。
     **保留策略不变**（仍是后者覆盖），改它会影响已有输出；报告只是让覆盖可见。
     实测随仓库分发的表 **809 个别名、0 冲突**，但仍加了回归锁——表一旦被改出重复定义，
     测试会立刻失败，而不是等到某次注释被归错类才发现。
  2. 修掉 `for line in open(...)` 的文件句柄泄漏：改为 `with` 读入。原写法会抛
     `ResourceWarning`，在 Windows 上还可能占着文件锁。
  3. `read_names` 的三种静默失败全部改为出声（新增可选参数 `issues`，不传仍保持旧行为）。
     - **纯数字行**：旧实现用 `tokens[-1:]` 兜底，把 `12345` 当成一个名字造出来。
       该兜底声称处理"数量在前的格式"，而那种行里一定有名字、`picked` 不会为空；
       `picked` 为空的唯一原因恰恰是"没有名字"。前提与场景自相矛盾，现已改为跳过并记录。
     - **配对模式遇到单名字行**：旧实现 `if len(picked) >= 2` 不成立就什么都不做，
       整行静默消失——文件看着有内容，统计里却没有它。
     - **单列模式遇到两列名字**：旧实现只取第一个，第二个凭空丢弃。
     三种情况现在都报到 stderr（含行号、原因、原文摘录、前 10 条上限），
     单列模式下还会额外提示"本文件像两列名字，试试 --pairs"。
  4. `tests.py` 42 → 57 条断言（新增 15 条：别名相关 9 条 + 解析可见性 6 条）。
- 2026-11-15（v1.3）`--pairs` 支持 `--weight`。
  1. `read_names` 配对模式返回三元组 `(名字A, 名字B, 次数)`：第三列是权重，缺省 1。
     之前把次数当噪声丢掉，导致「1 条分歧」与「7 条分歧」被当成同等严重。
  2. 新增 `rank_pairs()`：**合并镜像对**再排序。混淆矩阵有方向，`Tc1→PIF 7` 与 `PIF→Tc1 1` 会各占一行，它们是同一件事；
     合并后权重相加（8），方向明细保留（谁把谁吸走，决定该修哪一方）。
  3. 新增 `print_weighted()`：报「行数 / 对数 / 序列数」三个口径，并给出未加权与加权的判定占比对比。
  4. 配对 CSV 增加「次数」列（10 → 11 列）；CSV 保留原始行不合并，合并只出现在控制台。
  5. `--weight` 不带 `--pairs` 时报错退出，不静默忽略。
  6. `tests.py` 42 条断言（新增 9 条：权重解析 3、镜像合并 4、CLI 2）。
- 2026-11-15（v1.2）新增 `tests.py`，33 条断言全通过。**写测试时抓到 3 个真 bug**：
  1. `describe_so("non_LTR_retrotransposon")` 的 superfamily 仍是 `"non"`。根因：它同时以
     `_LTR_retrotransposon` 结尾，先被 LTR 分支截走，之前放在 LINE 分支里的"修复"是死代码。
     已改为在 LTR 分支之前单独处理。
  2. `classify_format("Unknown//ClassI")` 返回 `clade_style` 而非 `malformed`。根因同上：
     空级检查写在 `len(levels)==3` 之后，永远执行不到。已把空级检查提到最前。
  3. 路径 5 在"连 60% 相似都没有"时返回空建议——与"绝不静默丢弃"的原则冲突。
     现退回最接近的一个并标注 `← 低置信，最高相似度不足 0.60`。
  另修一处表述：没有输出 SO 名时不再打印"全部可回查 ✓"（把真空当通过），改为"跳过"。
  测试覆盖：格式判定 / 五条查找路径 / 层级归类 / 关系判定 / 行解析 / 不变量
  （命中键必须真在别名表里、输出必须可回查）/ 7 个真实文件回归 / CLI 输出与 CSV 列数 / 显示工具。
- 2026-11-15（v1.1）新增 `--verbose` / `-v`。
  1. `lookup()` 每条返回值新增 `hit` 字段＝"表里实际命中的那个键"。路径 1 等于输入本身；
     路径 2/3/4 与输入不同——这个差异才是需要复查的对象（`sine/alu` 实际命中 `SINE/Alu`）。
  2. `main` 在 `--verbose` 下输出三段式表格：路径名 / 原始名 / 命中键 / SO 名(SO ID)，
     并在表后给出五条路径各自的计数与说明；`--pairs --verbose` 同样列出配对中非精确命中的成员。
  3. 新增 `disp_width()` / `pad()` / `clip()`：按终端显示宽度对齐（中文占 2 列），
     超长名字截断为 `…`，避免长名把整列挤歪。
  4. CSV 报告新增"命中键"列（位于"状态"之后），列数 9 → 10。
  5. 默认输出格式不变，7 个测试文件回归结果与 v1 完全一致。
- 2026-11-15 `check_counts.py` 新增（模块 2）。第一版把 bed 直接累加得 9.38%，误判为"官方两个数字"；
  复核发现 EDTA 的 `TEanno.bed` 未合并（3461 行有 236 处相交/相接），改为按染色体排序合并后
  得 2,796,034 bp / 9.24%，与 RepeatMasker 只差 1,652 bp。现在脚本同时输出 raw 与 union 两列，
  并明确标注哪个是假数。**教训：宣称"两个权威数字矛盾"之前，必须先排除自己算错。**
- 2026-11-15（v1，依据 EDTA 注释阶段实测）
  1. 新增 `describe_so()`：给 105 个 SO 名补上「类别 / 角色 / 超家族」，未归类 0 个。
     顺带修掉两个归类错误：`non_LTR_retrotransposon` 被砍成超家族 "non"；
     7 个 rRNA 相关条目（18S/25S/5.8S/ETS/ITS）漏收进 RNA_GENE。
  2. 新增 `--pairs` 模式与 `relate()`：五种判定，见上表。
  3. 新增自反性检查：输出的 SO 名必须能被自己再吃回去，每次运行打印一行结论。
  4. 修 `lookup` 路径 3：原实现只砍最后一级，导致 TEsorter 风格
     `ClassI/LTR/Gypsy` 变成 `ClassI/LTR` 永远查不到。现同时尝试"砍最后一级"
     与"取最后两级"，TEsorter 五种写法 5/5 命中。
  5. CSV 报告新增 类别 / 角色 / 超家族 三列。
- 2026-11-15 `read_names` 修复计数列歧义：原实现固定取行首 token，对 `uniq -c`
  输出格式（数量在前，如 `1461 DNA/hAT-hAT19`）会把数字当名字。现改为跳过纯数字
  token（`isdigit()`），取第一个非数字 token，全部为数字时取行尾兜底。
  新增混合格式测试文件 `examples/test_names_counts.txt`；`with open` 规范写法同步应用。
- 2026-11-14 `load_te_so` 分隔符修复：TE_SO 文件中有 7 行（rRNA 相关，行 27-30/96-100）
  用连续空格而非制表符分隔，原实现按 `\t` 只切出 1 列而被静默丢弃。修复后
  `load_te_so` 别名数 776 → 809（去重后；含重复条目为 811，文件内自带
  2 个跨条目重复：`DNA/Kolobok`、`tRNA`），`rDNA/ITS1`、`cytosolic_25S_rRNA`
  等别名可正常命中。

## #542 别名 PR 素材（2026-11-16 生成）

61 个缺失名的分诊结果：**55 条可映射（按提交风险分三档），6 条明确不收**。

- **主补丁 21 条**：20 条结构后缀（`DNA/CACTA/MITE`、`DNA/hAT/nMITE`、`MITE/NA`
  等）+ `TIR/Novosib`（表里已有 `Novosib_TIR_transposon` 行，只缺 `TIR/` 写法）。
  这一档全是"剥掉末段就能落到已有超家族"，机械可定、与维护者立场无冲突。
- **判据型 4 条**：`nMITE`、`nMITE/NA`、`Unknown//ClassI`、`Unknown//ClassII`。
- **clade 备选 30 条**：28 条 LTR 亚家族（`LTR/Copia/Ale`、`LTR/Gypsy/CRM` 等）
  + `LTR/Copia/mixture`、`LTR/Gypsy/mixture`。单独成补丁，不并入主补丁。
- **6 条不收**：`ARTEFACT`（维护者原话 "the current default is probably good
  for now"，即维持兜底、不新增别名）、`LTR/PIF_Harbinger`（报告方把 TIR 超家族
  报成了 LTR）、`mixture` / `mixture/NA`（嵌合标记，不是分类名）、
  `nLTR` / `nLTR/NA`（non-LTR 与 nested LTR 歧义，需报告方澄清）。

**为什么 clade 那 30 条要单独拎出来。** 维护者 2026-03-02 回 yaoxkkk（61 条里 53 条的
来源）时写得很明确："clade-level classifications … is not the convention of RepeatMasker
library naming. You may want to strip off or reformat … before feeding the library to
EDTA2." 先前列为"高置信度"的依据只是表内已有 `LTR/Gypsy-Cigr`、`LINE/L1-Tx1` 等先例，
但先例不能推翻他刚表过的态。所以主补丁只放他不可能反对的部分，clade 另出一份备选，
`PR_BODY.md` 里按"可选、可独立退回"来写。

**改动原因留痕**：`ARTEFACT` 最初被我判为"可纳入"，理由是"沿用作者处理口径"；
回查 #542 第 11 条评论后发现原话是"ambiguous to me … current default is probably
good for now"，作者并未同意加别名，故改判为不收。**二手转述会走形，结论要回到原帖核对。**

**关键发现一：PR 的门槛不在代码，在"是否新增 SO 条目行"。**
`TE_Sequence_Ontology.txt` 头部写明：新增 SO 条目（新行）必须同步改
`bin/gff2bed.pl`。本次 55 个别名全部追加到**已有行的第三列**，不新增行，
所以 `gff2bed.pl` 不用动。这条必须写进 PR 正文，否则维护者会要求补改后打回。

**关键发现二：三档补丁分别实测，命中基线是 46/61，不是 0。**

| 补丁 | 别名数 | 工具命中 | 增量 |
|---|---|---|---|
| 不打补丁 | — | 46/61 | — |
| 主补丁 | 21 | 51/61 | +5 |
| 主补丁 + 判据型 | 25 | 55/61 | +9 |
| clade 备选 | 30 | 46/61 | **+0** |

clade 那 30 条对 `normalize.py` **零增量**——剥结构后缀、按亚家族前缀回退的兜底路径
本来就把它们全猜中了，且猜出的目标与本补丁逐条一致、无一冲突（这反过来验证了映射正确）。
也就是说：

- 对本工具，主补丁的增量是 5 条；clade 部分是纯冗余；
- 对**只按表精确查名**的消费者（EDTA 自身、TEsorter 等），增量才是 21 / 30 / 55 条，
  它们没有这层兜底。

值得记住的差别：**"我把名字认对了"和"上游表里有这个名字"是两件事。**
前者靠工具兜底，后者只能改数据。

三份补丁互不重叠：主补丁动 11 行、clade 动 2 行、交集 0，各自都能单独 `patch` 到原文件；
`main_plus + clade` 叠加后的结果与 `TE_Sequence_Ontology.with_clade.txt` 逐字节一致。

## 你的练习任务（学 Python 的正确姿势：改真实代码）

1. ~~【30 分钟】给 `lookup` 加连字符/下划线互换路径~~（已完成，见路径 4）；
2. ~~【1 小时】让 `main` 支持 `--verbose` 参数：命中路径 2/3/4 时打印原始名和命中名~~（已完成，见上节；顺带加了 `pad()`/`clip()` 按显示宽度对齐、CSV 新增"命中键"列）；
3. ~~【半天】写 `tests.py`：7 个文件变成断言用例 + `classify_format` 补 6 个断言~~
   （已完成，见 `tests.py`；共 33 条，并因此修掉 3 个 bug。下一个练习：给
   `tests.py` 加"用例数不能减少"的守卫，防止以后误删断言）；
4. ~~【半天】给 `--pairs` 加 `--weight` 参数，按混淆矩阵第 3 列加权~~
   （已完成，见上节 `rank_pairs` / `print_weighted`；顺带发现镜像对必须先合并，
   否则同一对分歧会占两个名次。~~下一个练习：把 `--weight` 的分母做成 `--total`~~
   已于 2026-11-17 完成，见 v1.4 与 `examples/rate_demo.out`；
   下一个练习：让 `--total` 可以从 `TEanno.sum` / 混淆矩阵直接读出来，免去手填）；
5. ~~【半天】从 `aliases_issue542.csv` 筛出缺失行，自动生成 PR 用的新增行文本
   （对照 TE_SO 文件中同 SO 组的现有行格式）~~
   （已完成，见 `make_alias_pr.py` 与 `pr_issue542/`；实际产出不是"61 行新条目"
   而是"别名追加到已有行的第三列"：主补丁 21 条动 11 行，全量 55 条动 15 行。
   并顺带发现 6 条不该收——详见上一节。
   下一个练习：把 `DECISIONS` 字典改成外部 YAML/JSON，让分诊规则可配。）
