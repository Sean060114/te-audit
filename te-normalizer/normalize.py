"""
te-normalizer —— 模块 0：跨工具 TE 命名规范化器（v1）

它做什么：
    输入一批 TE 分类名（例如来自 EDTA 报告、TEsorter 输出、RepeatModeler2 库），
    逐个判断 → 三种出路：
        1. 命中 TE_SO 别名表   -> 给出 SO 标准名 + SO ID（直接可用）
        2. 命名格式不合规      -> 标记为 clade_style / malformed（建议规范化重写）
        3. 完全未知            -> 用编辑距离给出 top-3 相似名建议（而非直接丢弃）

    v1 新增（依据 2026-11-15 的 EDTA 注释阶段实测）：
        - describe_so()：给每个 SO 名补上「类别 / 角色 / 超家族」三级信息
        - --pairs：读两列名字，判定「等价 / 同目不同科 / 跨类冲突」
        - 自反性检查：输出的 SO 名必须能被自己再吃回去
    v1.3 新增（--pairs --weight）：
        - read_names 在配对模式下读第三列次数作为权重（缺省 1）
        - rank_pairs()：合并镜像对（Tc1→PIF 7 与 PIF→Tc1 1 是同一件事），
          权重相加、方向明细保留，再按权重排序
        - print_weighted()：输出「行数 / 对数 / 序列数」三个口径，
          并对比未加权与加权的判定占比
    v1.1 新增：
        - --verbose：路径 2/3/4/5 逐条打印「原始名 -> 库内实际写法 -> SO 名+ID」，
          因为真正需要复查的恰恰是"命中键 != 原始名"的那些（大小写、剥级、连字符）

怎么跑：
    python normalize.py examples/test_names.txt
    python normalize.py examples/test_names.txt --verbose  # 逐条看：原始名->命中键->SO名
    python normalize.py examples/test_names.txt --out report.csv
    python normalize.py examples/pairs.txt --pairs          # 两列名字，做冲突判定

代码结构刻意保持"一个文件、六个函数"，每个函数都有中文注释。
"""

import argparse
import csv
import os
import re
import sys
import unicodedata                  # 只为算"显示宽度"，让中文列在终端里对齐
from difflib import SequenceMatcher, get_close_matches   # 标准库自带，无需安装

# 五条查找路径的编号与说明（--verbose 打印时用）
# 之所以给每条路径编号：用户看到 "matched_ci" 不知道发生了什么，
# 看到 "路径2 大小写" 才知道"它命中了，但命中的是表里另一个写法"。
PATH_INFO = {
    "matched":             ("路径1 精确",   "与表里写法一字不差"),
    "matched_ci":          ("路径2 大小写", "表里的实际写法与输入大小写不同"),
    "matched_after_strip": ("路径3 剥级",   "三级命名，剥掉一层后命中"),
    "matched_norm":        ("路径4 连字符", "下划线与连字符写法不同"),
    "unknown":             ("路径5 未命中", "只给相似建议，不丢弃"),
}


def disp_width(s):
    """终端显示宽度：中日韩全角字符占 2 列，西文占 1 列。

    不能直接用 len(s) 补齐：'{:<20}'.format('路径2 大小写') 会按 8 个字符算，
    但它在终端里占 14 列，后面的列就全歪了。
    """
    return sum(2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1
               for ch in str(s))


def pad(s, width):
    """按显示宽度右补空格，用于对齐表格"""
    s = str(s)
    return s + " " * max(0, width - disp_width(s))


def clip(s, width):
    """按显示宽度截断，超长时以 … 结尾（不截断会让整列歪掉）"""
    s = str(s)
    if disp_width(s) <= width:
        return s
    out, w = "", 0
    for ch in s:
        cw = 2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1
        if w + cw > width - 1:          # 留 1 列给省略号
            break
        out += ch
        w += cw
    return out + "…"

# ---------------------------------------------------------------------------
# 层级知识：TE_SO 表只有「标准名 / SO ID / 别名」三列，没有层级。
# 而判断"两个名字是不是同一件事"，光有 SO ID 不够 —— 必须知道它们各自属于
# 哪一个"目"（Class I 逆转座子 / Class II DNA 转座子）。下面这几张表就是补这个。
# 名词对应关系：目 = order = Class I/II；科 = superfamily = Copia/Gypsy/hAT...
# ---------------------------------------------------------------------------

# 子特征：不是"一个转座子"，而是转座子的零件（LTR 两端、插入时复制的宿主序列等）
SUB_FEATURE = {
    "long_terminal_repeat", "target_site_duplication", "terminal_inverted_repeat",
    "primer_binding_site", "RR_tract", "U_box",
}
REGION   = {"repeat_region"}          # 装零件的容器
FRAGMENT = {"repeat_fragment"}        # 残缺到判不出类别的碎片
SIMPLE_REPEAT = {"centromeric_repeat", "knob", "satellite_DNA", "telomeric_repeat",
                 "subtelomere", "low_complexity"}
ORGANELLE = {"chloroplast_DNA", "mitochondrial_DNA"}
RNA_GENE  = {"snRNA", "scRNA", "rRNA_gene", "rDNA_intergenic_spacer_element",
             "cytosolic_2S_rRNA", "cytosolic_5S_rRNA", "cytosolic_16S_rRNA",
             "cytosolic_23S_rRNA", "cytosolic_28S_rRNA",
             # 下面 7 个是 TE_SO 里同样属于 rRNA 基因簇的条目，早先漏收，
             # 会让它们落进 "unknown" 从而污染未命中率
             "cytosolic_18S_rRNA", "cytosolic_25S_rRNA", "cytosolic_5_8S_rRNA",
             "rRNA_3_external_transcribed_sequence", "rRNA_5_external_transcribed_sequence",
             "rRNA_internal_transcribed_spacer1", "rRNA_internal_transcribed_spacer2"}


# ---------------------------------------------------------------------------
# 函数 1：加载 TE_SO 映射表
# 文件格式（制表符分隔，# 开头是注释）：
# 第 1 列 = SO 标准名    第 2 列 = SO ID    第 3 列起 = 别名（逗号分隔）
# ---------------------------------------------------------------------------
def load_te_so(path):
    """加载 TE_SO 映射表。

    返回 (alias_to_so, so_names, conflicts) 三元组。

    conflicts 记录「同一别名指向不同 SO 名」的情形。这类冲突不会报错——
    旧实现是后者静默覆盖前者，查表结果取决于文件里的行走顺序，出了问题
    也看不出痕迹（实测当前表 809 个别名中冲突为 0，但那是数据干净，不是
    代码保证了）。现在把它显式收集起来，由调用方报给人看。

    保留策略保持不变：仍是后者覆盖前者。改它会改变已有输出，而报告只是
    让覆盖可见，不是修正它。
    """
    alias_to_so = {}          # 别名 -> (SO标准名, SO ID)
    so_names = set()          # 所有 SO 标准名（建议功能要用）
    conflicts = {}            # 别名 -> [(SO标准名, SO ID), ...]，仅记录冲突
    # 表只有一百多行，一次性读入即可；用 with 是为了显式关闭文件——
    # 原写法 for line in open(...) 依赖 GC 回收，会抛 ResourceWarning，
    # 在 Windows 上还可能占着文件锁。
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    for line in lines:
        # 跳过无效行
        if line.startswith("#") or not line.strip():
            continue
        # 把一行文本变成干净的列清单
        cols = []
        for c in line.rstrip("\n").split("\t"):             # 按制表符分割，可能会形成空串
            if c.strip():                                   # 剥掉首尾空白后非空才保留
                cols.append(c.strip())
        # 跳过格式不完整的行
        if len(cols) < 3:
            # 用连续空格而非制表符分隔（作者手写不一致）。按 \t 只能切出 1 列，
            # 会被上面的 len<3 静默丢弃，损失 35 个别名。
            # 用"任意连续空白"重切一次即可救回。别名串内部用 ", " 分隔，
            # split() 会把它拆开——但 len<3 的兜底分支本来就只救空格分隔行，
            # 正常行不会走到这里
            cols = line.split()
        if len(cols) < 3:
            continue          # 两次解析都不完整，确认是坏行，跳过
        so_name, so_id, alias_part = cols[0], cols[1], cols[2:]
        so_names.add(so_name)                 # set集合可自动去重
        for group in alias_part:
            for alias in re.split(r"[,\s]+", group):        # 把这一组拆成单个别名
                if alias:                                   #  过滤空串
                    # 覆盖前先记账：同一个别名指向不同 SO 名时，后者会覆盖前者，
                    # 若不记录，这条覆盖就永远不可见。
                    prev = alias_to_so.get(alias)
                    if prev and prev[0] != so_name:
                        lst = conflicts.setdefault(alias, [prev])
                        if (so_name, so_id) not in lst:
                            lst.append((so_name, so_id))
                    alias_to_so[alias] = (so_name, so_id)
    return alias_to_so, so_names, conflicts


# ---------------------------------------------------------------------------
# 函数 2：给 SO 标准名补上「类别 / 角色 / 超家族」
# 为什么需要（2026-11-15 实测依据）：
#   1) EDTA 的 TEanno.gff3 里，同一个 classification=LTR/Copia 同时对应 4 个 SO ID
#      （2264 元件本体、0286 它的 LTR、0434 它的 TSD、0657 装它的容器）。
#      不看"角色"，就无法把 classification 映射到唯一的 SO ID。
#   2) 注释行 3666 行 vs 元素级 3461 行，差的 205 个正是 role=sub_feature/region。
#      不会区分角色，就会把零件当整机数。
# ---------------------------------------------------------------------------
def describe_so(so_name):
    n = so_name or ""
    if not n:
        return {"order": "unknown", "role": "unknown", "superfamily": ""}
    if n in SUB_FEATURE:
        return {"order": "-", "role": "sub_feature", "superfamily": n}
    if n in REGION:
        return {"order": "-", "role": "region", "superfamily": n}
    if n in FRAGMENT:
        return {"order": "-", "role": "fragment", "superfamily": n}
    if n in SIMPLE_REPEAT:
        return {"order": "-", "role": "simple_repeat", "superfamily": n}
    if n in ORGANELLE:
        return {"order": "-", "role": "organelle", "superfamily": n}
    if n in RNA_GENE:
        return {"order": "-", "role": "RNA_gene", "superfamily": n}

    # 必须放在 YR 规则之前：Crypton 名字里带 YR，但它属于 Class II（DNA 转座子）
    if n in ("Crypton_YR_transposon", "polinton"):
        return {"order": "Class II", "role": "element", "superfamily": n.split("_")[0]}

    # ---- Class I：逆转座子（RNA 中介）----
    # 必须放在 LTR 分支之前：non_LTR_retrotransposon 这个名字以 _LTR_retrotransposon
    # 结尾，若先走 LTR 分支会被截成 "non"（假超家族）。2026-11-15 写测试时才发现，
    # 之前放在 LINE 分支里的那段兜底是死代码，永远执行不到。
    if n == "non_LTR_retrotransposon":
        return {"order": "Class I", "role": "element", "superfamily": "non_LTR"}
    if n.endswith("_LTR_retrotransposon") or n in ("LTR_retrotransposon", "TRIM", "LARD"):
        return {"order": "Class I", "role": "element",
                "superfamily": n.replace("_LTR_retrotransposon", "")}
    if n.endswith("_YR_retrotransposon") or n == "YR_retrotransposon":
        return {"order": "Class I", "role": "element",
                "superfamily": n.replace("_YR_retrotransposon", "") or "YR"}
    if n.endswith("_LINE_retrotransposon") or n in ("LINE_element",
                                                    "Penelope_retrotransposon"):
        return {"order": "Class I", "role": "element",
                "superfamily": n.replace("_LINE_retrotransposon", "")
                                 .replace("_retrotransposon", "")}
    if n.endswith("_SINE_retrotransposon") or n == "SINE_element":
        return {"order": "Class I", "role": "element",
                "superfamily": n.replace("_SINE_retrotransposon", "")}
    if n in ("pararetrovirus", "ERTBV_retrotransposon"):
        return {"order": "Class I", "role": "element", "superfamily": "pararetrovirus"}
    if n == "retrotransposon":
        return {"order": "Class I", "role": "element", "superfamily": "unknown"}

    # ---- Class II：DNA 转座子（DNA 中介）----
    if n.endswith("_TIR_transposon") or n in ("DNA_transposon", "MITE",
                                              "terminal_inverted_repeat_element"):
        return {"order": "Class II", "role": "element",
                "superfamily": n.replace("_TIR_transposon", "")}
    if n == "helitron":
        return {"order": "Class II", "role": "element", "superfamily": "helitron"}

    return {"order": "unknown", "role": "unknown", "superfamily": n}


# ---------------------------------------------------------------------------
# 函数 3：命名格式合法性检查
# RepeatMasker 库的规范是"最多两级"：Class/Superfamily，如 DNA/hAT、LTR/Copia。
# 三级以上（如 LTR/Copia/SIRE）是 clade 级命名 —— 不是错误，但不合规范，
# ---------------------------------------------------------------------------
def classify_format(name):
    levels = name.split("/")
    # 空级检查必须放在最前面。写在 len==3 之后的话，Unknown//ClassI 会被
    # 当成"三级命名"判成 clade_style（它是坏的，不是 clade 级命名）。
    # 2026-11-15 写断言时发现：原来的顺序让这条检查成了死代码。
    if any(not lv for lv in levels):
        return "malformed"          # 出现空级，如 "Unknown//ClassI"
    if len(levels) > 3:
        return "malformed"          # 超过三级，疑似拼接错误
    if len(levels) == 3:
        return "clade_style"        # 三级：clade 级命名，建议剥离第三级
    return "standard"


# ---------------------------------------------------------------------------
# 函数 4：核心查找 —— 精确 -> 大小写不敏感 -> 规范化重写 -> 相似建议
# 返回 dict：{status, so_name, so_id, suggestion}
# ---------------------------------------------------------------------------
def lookup(name, alias_to_so, so_names):
    # 每条返回值都带 hit："表里实际命中的那个键"。
    # 路径 1 时 hit 就等于输入本身；路径 2/3/4 时 hit 与输入不同 ——
    # 这个差异正是要给人看的（输入写成 DNA/hAt，表里其实是 DNA/hAT）。
    # 路径 1：精确命中
    if name in alias_to_so:
        return {"status": "matched", "so_name": alias_to_so[name][0],
                "so_id": alias_to_so[name][1], "suggestion": "", "hit": name}

    # 路径 2：大小写不敏感命中（TE 命名大小写混乱是常态）
    lower_map = {}
    for k in alias_to_so:
        lower_map[k.lower()] = k
    if name.lower() in lower_map:
        real = lower_map[name.lower()]
        return {"status": "matched_ci", "so_name": alias_to_so[real][0],
                "so_id": alias_to_so[real][1], "suggestion": real, "hit": real}

    # 路径 3：三个层级的名字 -> 剥掉一层再查。两种层级顺序都存在，必须都试：
    #   RepeatMasker 风格 LTR/Copia/SIRE ：超家族在前，第三级是 clade -> 砍掉最后一级
    #   TEsorter 风格  ClassI/LTR/Gypsy  ：目/亚类在前，超家族在最后 -> 取最后两级
    # 只砍最后一级的话，ClassI/LTR/Gypsy 会变成 ClassI/LTR（不是超家族），永远查不到。
    if classify_format(name) == "clade_style":
        levels = name.split("/")
        for cand in ("/".join(levels[:-1]),          # LTR/Copia/SIRE -> LTR/Copia
                     "/".join(levels[-2:])):         # ClassI/LTR/Gypsy -> LTR/Gypsy
            if cand in alias_to_so:
                return {"status": "matched_after_strip", "so_name": alias_to_so[cand][0],
                        "so_id": alias_to_so[cand][1], "suggestion": cand, "hit": cand}

    # 路径 4：连字符/下划线差异 —— 两个方向都试，因为表里两种写法都存在
    # 例：输入 DNA/CMC_Transib，表里是 DNA/CMC-Transib（反向情形如 18S-rRNA / 18S_rRNA）
    for variant in (name.replace("-", "_"), name.replace("_", "-")):
        if variant != name and variant in alias_to_so:
            return {"status": "matched_norm", "so_name": alias_to_so[variant][0],
                    "so_id": alias_to_so[variant][1], "suggestion": variant, "hit": variant}

    # 路径 5：编辑距离建议 —— 给 top-3，而不是像 EDTA 那样直接降级丢弃
    scored = []
    for cand in set(list(alias_to_so) + list(so_names)):
        ratio = SequenceMatcher(None, name, cand).ratio()
        scored.append((ratio, cand))
    scored.sort(reverse=True)
    picked = [f"{cand}({ratio:.2f})" for ratio, cand in scored[:3] if ratio >= 0.6]
    if not picked and scored:
        # 一个都没到 60% 时也必须给点东西：空建议等于静默丢弃，与设计原则冲突。
        # 这里退回"最接近的那个"，并明确标注置信度不足，由人决定要不要采信。
        ratio, cand = scored[0]
        picked = [f"{cand}({ratio:.2f}) ← 低置信，最高相似度不足 0.60"]
    return {"status": "unknown", "so_name": "", "so_id": "",
            "suggestion": "; ".join(picked), "hit": ""}


# ---------------------------------------------------------------------------
# 函数 5：读入一批名字
# 模式 A（默认，每行一个名字）：支持 "名字 数量" 或 "数量 名字" 两种真实格式
# 模式 B（--pairs，每行两个名字）：支持 "名A  名B" 或 "名A  名B  次数"
# ---------------------------------------------------------------------------
def resolve_input(path):
    """把命令行给的路径解析成真实文件；失败时返回 (None, 人话提示)。

    为什么要这一层：直接 open() 抛 FileNotFoundError，用户看到的是八行
    traceback，最后一行才是 "No such file"。而这类错误九成只有两种：
      1. 记错文件名（pairs.txt 还是 pairs_matrix.txt）
      2. 站错目录（在仓库根目录跑，而 normalize.py 在 te-normalizer/ 下）
    这两种都不该让人去读栈。
    """
    if os.path.isfile(path):
        return path, ""
    if os.path.isdir(path):
        return None, (f"{path} 是目录，不是文件。\n"
                      f"请指向具体文件，例如 {path.rstrip('/\\')}/pairs_matrix.txt")

    def show(p):                       # 提示里的路径统一用正斜杠：与用户输入的习惯一致
        return p.replace("\\", "/")

    base = os.path.basename(path)
    here = os.path.dirname(os.path.abspath(__file__))
    probes = [os.path.dirname(path) or ".", os.path.join(here, "examples")]

    lines = [f"找不到输入文件：{path}"]
    seen, guesses, listings = set(), [], []
    for d in probes:
        if not os.path.isdir(d):
            continue
        files = sorted(n for n in os.listdir(d)
                       if os.path.isfile(os.path.join(d, n)) and not n.startswith("."))
        for g in get_close_matches(base, files, n=3, cutoff=0.5):
            p2 = os.path.join(d, g)
            if p2 not in seen:
                seen.add(p2)
                guesses.append(p2)
        txts = [n for n in files if n.endswith(".txt")]
        if txts:
            listings.append((d, txts))

    if guesses:
        lines.append("")
        lines.append("是不是想用：")
        lines += [f"  {show(g)}" for g in guesses[:3]]
    if listings:
        d, txts = listings[0]
        lines.append("")
        lines.append(f"{show(d)} 里可用的输入文件（约定：.txt 是输入，.out / .csv 是产物）：")
        lines += [f"  {n}" for n in txts[:10]]
        if len(txts) > 10:
            lines.append(f"  …… 其余 {len(txts) - 10} 个略")
    lines.append("")
    lines.append("提示：--pairs 需要两列（名字A <TAB> 名字B [次数]）；一列名字请去掉 --pairs。")
    return None, "\n".join(lines)


def read_names(path, pairs=False, issues=None):
    """读入一批名字。

    模式 A（默认，每行一个名字）：支持 "名字 数量" 或 "数量 名字" 两种真实格式
    模式 B（--pairs，每行两个名字）：支持 "名A  名B" 或 "名A  名B  次数"

    issues 是可选的输出列表；传入时把"读不对的行"追加进去，每项为
    (行号, 原文, 原因)。不传就保持旧行为（静默），传了就让调用方把问题报出去。

    为什么需要 issues：解析失败的代价是"清单里少了几行"，而下游拿到的
    只是一份更短的列表，它无从判断少了什么。不报出来，这些行就等于
    凭空消失了——比读错更难查。
    """
    names = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = re.split(r"\s+", line)
            # 只保留非数字 token：行内可能混有计数值，三种真实格式都存在：
            # "LTR/Copia/SIRE 104231"    （EDTA 报告摘录：名字在前）
            # "104231 LTR/Copia/SIRE"    （uniq -c 输出：数量在前）
            # "DNA/Tc1_Mariner  DNA/PIF_Harbinger  7"  （混淆矩阵：两个名字 + 次数）
            picked = [t for t in tokens if not t.isdigit()]
            if not picked:
                # 整行没有任何非数字 token —— 这一行根本没有名字。
                # 旧实现在这里用 tokens[-1:] 兜底，会把纯数字当成名字造出来：
                # 该兜底声称处理的是"数量在前的格式"（那种行里还有名字，
                # picked 不会为空），而 picked 为空的唯一原因恰恰是"没有名字"。
                # 前提与场景自相矛盾，所以这里跳过而不是猜。
                if issues is not None:
                    issues.append((lineno, line, "整行全是数字，没有名字可读，已跳过"))
                continue
            if pairs:
                if len(picked) >= 2:
                    # 行尾若带次数，就把它当权重。为什么不忽略它：
                    # EDTA 混淆矩阵里 Tc1→PIF 是 7 条、CACTA→Mutator 是 1 条，
                    # 两者都只是"一对分歧"，不加权就会当成同等严重。
                    # 没有次数行时权重记 1——一条分歧仍然算一条，只是最轻。
                    weight = next((int(t) for t in reversed(tokens) if t.isdigit()), 1)
                    names.append((picked[0], picked[1], weight))
                    if len(picked) > 2 and issues is not None:
                        issues.append((lineno, line,
                                       f"配对模式只取前两个名字，忽略 {picked[2:]}"))
                else:
                    # 配不成对。旧实现让 len<2 的行直接什么都不做，
                    # 整个行静默消失——文件看着有内容，统计里却没有它。
                    if issues is not None:
                        issues.append((lineno, line,
                                       "配对模式需要两个名字，本行只有 1 个，已跳过"))
            else:
                names.append(picked[0])
                if len(picked) > 1 and issues is not None:
                    issues.append((lineno, line,
                                   f"单列模式只取第一个名字，忽略 {picked[1:]}"))
    return names


# ---------------------------------------------------------------------------
# 函数 6：两个名字之间的关系判定（--pairs 模式的核心）
# 为什么要它（2026-11-15 实测依据）：
#   EDTA 已经算了一个 Inconsistent_rate，但它的判定条件只到"目"级
#   （cleanup_misclas.pl 把 Tc1_Mariner→DTT、PIF_Harbinger→DTH，再都归成 DNA）。
#   所以 14.89% 的 Tc1↔PIF 分歧被记录后没有任何处置。
#   要补上这一步，先得能机械地说清楚：这两个名字，到底是同一件事、还是同目不同科。
# ---------------------------------------------------------------------------
def relate(a, b):
    if a["status"] == "unknown" or b["status"] == "unknown":
        return "无法判定", "至少一方未命中 TE_SO"
    if a["so_id"] == b["so_id"]:
        return "等价", f"归一后同为 {a['so_name']}"
    da, db = describe_so(a["so_name"]), describe_so(b["so_name"])
    if da["role"] != db["role"]:
        return "角色不同", f"{da['role']} vs {db['role']}（不能直接比较）"
    if da["order"] != db["order"]:
        return "跨类冲突", f"{da['order']} vs {db['order']}（Class I / II 都不同，必然有一方错）"
    if da["superfamily"] != db["superfamily"]:
        return "同目不同科", (f"同属 {da['order']}，但 {da['superfamily']} vs "
                             f"{db['superfamily']}（EDTA 只记录、不处置的那一类）")
    return "同科异名", f"同一超家族 {da['superfamily']} 的两种 SO 写法"


# ---------------------------------------------------------------------------
# 函数 7：按权重排序分歧对，并把镜像对合并（--pairs --weight 的核心）
# 为什么必须合并镜像对：
#   混淆矩阵是有方向的（行＝EDTA 原始分类，列＝TEsorter 重判结果），所以
#   "Tc1→PIF 7 条" 和 "PIF→Tc1 1 条" 各占一行。它们说的是同一件事——
#   这两个超家族被互相判错——只是方向相反。不合并的话，同一对分歧会占据
#   排行榜两个名次，把真正的第三名挤到看不见。
#   合并后次数相加（7+1=8），但方向明细保留：方向能说明"谁在系统性地把谁
#   吸走"，这是判断该修哪一方分类器的依据，不能丢。
#
# 两个坑（2026-11-19 修）：
#   1) 身份必须分命名空间。未命中的一侧退回原始名，而 SO ID 也是一个字符串，
#      两套东西混在一个命名空间里，输入里出现 "SO:0002278" 这种名字时就会与
#      真正的 ID 撞成同一个键，把"等价"和"无法判定"两类记录并成一组。
#      退回时加 "name:" 前缀，两者永不相等。
#   2) 一个组里可能有不止一种判定（撞键、或将来 relate 引入方向性判定）。
#      旧实现 groups.setdefault(key, {"verdict": r["verdict"]}) 只取**最先进入
#      那一行**的判定，同一份数据换行序就换标签，且毫无痕迹。现在收集全部判定，
#      只有一种时用那一种，多于一种时标为"混合"，不替数据选边。
# ---------------------------------------------------------------------------
def _ident(side_id, side_name):
    """配对身份：优先用 SO ID；未命中时退回原始名，但加命名空间前缀"""
    return side_id or ("name:" + side_name)


def _short(name, resolved):
    """方向明细里的短名：已解析的取最后一级（LTR/Copia/Ale → Ale）。

    未命中的名字**不能截**：`DNA/CACTA/` 会截成空串（箭头两边空），
    `Unknown//ClassI` 会截成 `ClassI`——把一个"未知类"显示成一个类别名，
    比不显示更糟。
    """
    if not resolved.get(name):
        return name
    return name.split("/")[-1] or name


def rank_pairs(recs):
    groups, resolved = {}, {}
    for r in recs:
        resolved[r["a"]] = bool(r["a_id"])
        resolved[r["b"]] = bool(r["b_id"])
        key = tuple(sorted((_ident(r["a_id"], r["a"]),
                            _ident(r["b_id"], r["b"]))))
        g = groups.setdefault(key, {"weight": 0, "verdicts": {}, "dirs": {}})
        g["weight"] += r["w"]
        g["verdicts"][r["verdict"]] = g["verdicts"].get(r["verdict"], 0) + 1
        g["dirs"][(r["a"], r["b"])] = g["dirs"].get((r["a"], r["b"]), 0) + r["w"]

    out = []
    for g in groups.values():
        # 平局按名字兜底：只按权重排会退化成"插入顺序"（sorted 稳定），
        # 同一份数据换个行序，a/b 就互换，报告里的行标识跟着漂。
        dirs = sorted(g["dirs"].items(), key=lambda x: (-x[1], x[0]))
        # 单向也写箭头：原来的 "单向 2" 配上分歧对那列的 "↔"，读者看不出
        # 是谁被判成了谁，而这正是"该修哪一方"的依据。
        detail = " : ".join(f"{_short(a, resolved)}→{_short(b, resolved)} {w}"
                            for (a, b), w in dirs)
        vs = g["verdicts"]
        if len(vs) == 1:
            verdict, mixed = next(iter(vs)), ""
        else:
            verdict = "混合"
            mixed = " / ".join(f"{k}×{n}" for k, n in
                               sorted(vs.items(), key=lambda x: -x[1]))
        out.append({"a": dirs[0][0][0], "b": dirs[0][0][1],
                    "weight": g["weight"], "verdict": verdict,
                    "mixed": mixed, "dirs": detail})
    out.sort(key=lambda x: (-x["weight"], x["verdict"], x["a"]))
    return out


# ---------------------------------------------------------------------------
# 函数 8：把加权结果打印出来（--weight 的输出层）
# 三块：总量 -> 判定分布（未加权 vs 加权）-> 按权重排序的分歧对
# 为什么第一块要报"行数 / 对数 / 序列数"三个数：
#   8 行、7 对、17 条——同一份数据三个口径。这本身又是模块 2 那个教训的重演，
#   报数的人必须说清用的是哪一个。
# ---------------------------------------------------------------------------
def print_weighted(recs):
    """--weight 的输出层。

    口径只有一处：**一切以 rank_pairs 合并后的"对"为准**。
    旧实现一边报"合并镜像后 6 对"、一边在表里把行数当对数（合计 8），
    同一屏两个数打架——这是本项目第四次栽在"分母是什么"上。
    """
    allrank = rank_pairs(recs)
    ranked = [r for r in allrank if r["verdict"] != "等价"]   # 等价对不是分歧，不进排序
    n_row, n_pair = len(recs), len(allrank)

    gtally, gwtally = {}, {}
    for g in allrank:
        k = g["verdict"]
        gtally[k] = gtally.get(k, 0) + 1
        gwtally[k] = gwtally.get(k, 0) + g["weight"]
    total_w = sum(gwtally.values())

    print(f"共 {n_row} 行记录，合并镜像后 {n_pair} 对"
          f"（其中等价 {gtally.get('等价', 0)} 对），涉及序列 {total_w} 条")
    print()
    print("判定分布：未加权（按对数） vs 加权（按次数）")
    print("  " + pad("判定", 16) + pad("对数", 6) + pad("占比", 9)
          + pad("序列数", 8) + pad("占比", 9) + "变化")
    print("  " + "-" * 54)
    keys = sorted(gtally, key=lambda k: (-gwtally[k], -gtally[k]))
    shifts = []
    for k in keys:
        pa, pw = gtally[k] / n_pair * 100, gwtally[k] / total_w * 100
        shifts.append((pw - pa, k))
        flag = f"  {pw - pa:+.1f}pp" if abs(pw - pa) >= 0.05 else ""
        print("  " + pad(k, 16) + pad(gtally[k], 6) + pad(f"{pa:.1f}%", 9)
              + pad(gwtally[k], 8) + pad(f"{pw:.1f}%", 9) + flag)
    print("  " + "-" * 54)
    print("  " + pad("合计", 16) + pad(n_pair, 6) + pad("100.0%", 9)
          + pad(total_w, 8) + "100.0%")

    print()
    print("按涉及序列数排序（只列非等价对；镜像对已合并，末列是方向明细）")
    print("  " + pad("序列数", 8) + pad("判定", 16) + pad("分歧对", 46) + "方向")
    print("  " + "-" * 76)
    if not ranked:
        print("  （没有非等价对）")
    for r in ranked:
        pair = f"{r['a']} ↔ {r['b']}"
        tail = r["dirs"]
        if r.get("mixed"):
            # 组内判定不一致：不隐藏，直接写在行尾。这类行意味着"同一个键下
            # 挤进了两种关系"，比单纯的方向明细更值得怀疑。
            tail = f"{tail}  ｜ 判定混合：{r['mixed']}"
        print("  " + pad(r["weight"], 8) + pad(r["verdict"], 16)
              + pad(clip(pair, 44), 46) + tail)

    # 只在加权确实改变了结论时才多说一句，避免刷无意义的提示
    d, k = max(shifts, key=lambda x: abs(x[0]))
    if abs(d) >= 5:
        print()
        print(f"读法：{k} 未加权时占 {gtally[k] / n_pair * 100:.1f}%（{gtally[k]} 对），"
              f"加权后占 {gwtally[k] / total_w * 100:.1f}%（{gwtally[k]} 条），差 {d:+.1f}pp。")
        print("      差异来自少数几个大格子——把每一对分歧都当成同等严重，会低估它们。")


# ---------------------------------------------------------------------------
# 函数 8b：把"加权不一致率"算出来（--total 的输出层）
#
# 为什么需要它：EDTA 自己会打印一个 Inconsistent_rate，分母是**库里全部序列数**
# （本次 Chr1 实测 17 / 1593 = 1.07%）。我们的 --weight 只报"分歧对占多少对 /
# 涉及多少条"，分母是"出现分歧的那部分"，与 EDTA 的数不可比——两个都在 1% 附近
# 的百分数，分母差一个数量级。给了 --total 之后两边才落在同一把尺子上。
#
# 分子分两档，因为"不一致"有两种含义：
#   全部非等价对     ——含"角色不同"（element vs sub_feature），可比 EDTA 的口径
#   只算真冲突       ——去掉"角色不同"和"等价"，才是真正的判断分歧
# 把两个都报出来，读者自己选，不要替他决定哪个才算"错"。
# ---------------------------------------------------------------------------
def print_rate(wtally, total, baseline=None):
    total_w = sum(wtally.values())
    agree = wtally.get("等价", 0)
    all_bad = total_w - agree
    role_only = wtally.get("角色不同", 0)
    real_bad = all_bad - role_only

    print()
    print("不一致率（与 EDTA 同一分母）")
    print("  " + pad("库序列总数 N", 20) + pad(total, 10))
    print("  " + pad("其中等价", 20) + pad(agree, 10))
    print("  " + pad("非等价", 20) + pad(all_bad, 10)
          + f"  → {all_bad / total * 100:.2f}%")
    print("  " + pad("只算真冲突", 20) + pad(real_bad, 10)
          + f"  → {real_bad / total * 100:.2f}%"
          + (f"（另有 {role_only} 条仅角色不同，不计入）" if role_only else ""))
    if baseline is not None:
        diff = all_bad / total * 100 - baseline
        dtxt = "一致" if abs(diff) < 0.005 else f"差 {diff:+.2f} pp"
        print("  " + pad("EDTA 自报", 20) + pad(f"{baseline:.2f}%", 10) + f"  {dtxt}")
        print("  读法：两边同分母可比；差异只能来自分桶规则，不是数据。")


# ---------------------------------------------------------------------------
# 函数 9：主流程 —— 解析参数 -> 加载 -> 逐个查 -> 打印/写 CSV
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="TE 命名规范化器 v1")
    ap.add_argument("input", help="待检查名称文件（每行一个，可带数量列）")
    ap.add_argument("--te-so", default="data/TE_Sequence_Ontology.txt",
                    help="TE_SO 映射表路径")
    ap.add_argument("--out", default="", help="结果 CSV 输出路径（可选）")
    ap.add_argument("--pairs", action="store_true",
                    help="两列名字模式：判定每一对的等价/冲突关系")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="逐条打印命中明细：原始名 -> 命中的库内写法 -> SO 名与 SO ID"
                         "（路径 2/3/4/5；路径 1 精确命中只给计数）")
    ap.add_argument("--weight", action="store_true",
                    help="配合 --pairs：把第三列的次数当权重，合并镜像对后按权重排序，"
                         "并对比未加权/加权的判定分布")
    ap.add_argument("--total", type=int, default=None,
                    help="配合 --pairs：库/混淆矩阵的总序列数（EDTA 报 Inconsistent_rate 的"
                         "分母）。给了才输出「加权不一致率」，与 EDTA 自报值同分母可比")
    ap.add_argument("--baseline", type=float, default=None,
                    help="配合 --total：EDTA 自报的不一致率（%%），用于算差值")
    args = ap.parse_args()
    if args.weight and not args.pairs:
        ap.error("--weight 只能配合 --pairs 使用（权重来自分歧对的第三列次数）")
    if (args.total is not None or args.baseline is not None) and not args.pairs:
        ap.error("--total / --baseline 只能配合 --pairs 使用（分子来自配对判定）")
    if args.baseline is not None and args.total is None:
        ap.error("--baseline 需要同时给 --total：没有分母就算不出可比的不一致率")
    if args.total is not None and args.total <= 0:
        ap.error("--total 必须是正整数（它是库序列总数，不是分歧条数）")

    alias_to_so, so_names, conflicts = load_te_so(args.te_so)
    if conflicts:
        print(f"警告：{len(conflicts)} 个别名指向多个 SO 名，查表取文件中最后出现的那个，"
              f"结果依赖行走顺序：", file=sys.stderr)
        for a in sorted(conflicts)[:10]:
            forms = " / ".join(n for n, _ in conflicts[a])
            print(f"  {a!r} -> {forms}（最终取 {alias_to_so[a][0]}）", file=sys.stderr)
        if len(conflicts) > 10:
            print(f"  …… 其余 {len(conflicts) - 10} 个略", file=sys.stderr)
        print("  建议在 TE_Sequence_Ontology.txt 里消除重复：否则同一条注释在两个版本间"
              "可能归到不同的 SO 名。", file=sys.stderr)
    real_input, why = resolve_input(args.input)
    if real_input is None:
        print(why, file=sys.stderr)
        sys.exit(2)

    read_issues = []
    rows = read_names(real_input, pairs=args.pairs, issues=read_issues)
    if read_issues:
        mode = "配对" if args.pairs else "单列"
        print(f"警告：{len(read_issues)} 行无法按{mode}模式完整解析：", file=sys.stderr)
        for ln, txt, why in read_issues[:10]:
            print(f"  行{ln}: {why}", file=sys.stderr)
            print(f"        原文: {clip(txt, 68)}", file=sys.stderr)
        if len(read_issues) > 10:
            print(f"  …… 其余 {len(read_issues) - 10} 行略", file=sys.stderr)
        print("  这些行没有进入统计，上面的计数可能偏低。"
              + ("（模式用错了？本文件像两列名字，试试 --pairs）" if not args.pairs else ""),
              file=sys.stderr)

    if args.pairs:
        recs, non_exact = [], []        # non_exact 供 --verbose 用
        for a_name, b_name, w in rows:
            ra = lookup(a_name, alias_to_so, so_names)
            rb = lookup(b_name, alias_to_so, so_names)
            for nm, r in ((a_name, ra), (b_name, rb)):
                if r["status"] != "matched":
                    non_exact.append((nm, r))
            verdict, why = relate(ra, rb)
            recs.append({"a": a_name, "a_so": ra["so_name"], "a_id": ra["so_id"],
                         "a_order": describe_so(ra["so_name"])["order"],
                         "b": b_name, "b_so": rb["so_name"], "b_id": rb["so_id"],
                         "b_order": describe_so(rb["so_name"])["order"],
                         "verdict": verdict, "why": why, "w": w})

        tally, wtally = {}, {}
        for r in recs:
            tally[r["verdict"]] = tally.get(r["verdict"], 0) + 1
            wtally[r["verdict"]] = wtally.get(r["verdict"], 0) + r["w"]

        if args.weight:
            print_weighted(recs)
        else:
            print(f"共 {len(recs)} 对名字：" +
                  " ｜ ".join(f"{k} {v}" for k, v in sorted(tally.items(), key=lambda x: -x[1])))
            for r in recs:
                if r["verdict"] != "等价":
                    print(f"  [{r['verdict']:<6}] {r['a']:<20} vs {r['b']:<20} -> {r['why']}")
        if args.total is not None:
            involved = sum(wtally.values())
            if involved > args.total:
                print(f"警告：配对记录共涉及 {involved} 条序列，超过 --total {args.total}。"
                      "分母给错了？不一致率未输出。", file=sys.stderr)
            else:
                print_rate(wtally, args.total, args.baseline)
        if args.verbose and non_exact:
            print("\n" + "-" * 70)
            print("配对模式下的归一明细（--verbose）")
            print("-" * 70)
            print("  " + pad("路径", 14) + pad("原始名", 26) + pad("命中键", 26) + "SO 标准名 (SO ID)")
            for nm, r in non_exact:
                label = PATH_INFO[r["status"]][0]
                hit = r["hit"] or "—"
                target = (f"{r['so_name']} ({r['so_id']})" if r["so_name"]
                          else "建议：" + (r["suggestion"] or "无相似名"))
                print("  " + pad(label, 14) + pad(clip(nm, 25), 26) + pad(clip(hit, 25), 26) + target)
        if args.out:
            with open(args.out, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                # 次数列放在最后：不加 --weight 时全为 1，加了的才有意义。
                # 保留原始行而不是合并后的结果——CSV 是审计用的流水，
                # 合并在控制台看即可，别把原始记录改掉。
                w.writerow(["名字A", "A_SO名", "A_SO_ID", "A类别",
                            "名字B", "B_SO名", "B_SO_ID", "B类别", "判定", "说明", "次数"])
                for r in recs:
                    w.writerow([r["a"], r["a_so"], r["a_id"], r["a_order"],
                                r["b"], r["b_so"], r["b_id"], r["b_order"],
                                r["verdict"], r["why"], r["w"]])
            print(f"报告已写入 {args.out}")
        return 0

    results = []
    for name in rows:
        r = lookup(name, alias_to_so, so_names)
        r["name"] = name
        r["format"] = classify_format(name)
        r.update(describe_so(r["so_name"]))     # 补 order / role / superfamily
        results.append(r)

    # 控制台摘要
    n = {"matched": 0, "matched_ci": 0, "matched_after_strip": 0,
         "matched_norm": 0, "unknown": 0}
    for r in results:
        n[r["status"]] += 1
    print(f"共 {len(results)} 个名称："
          f"精确命中 {n['matched']} ｜ 大小写命中 {n['matched_ci']} ｜ "
          f"剥级后命中 {n['matched_after_strip']} ｜ 连字符归一 {n['matched_norm']} ｜ "
          f"未知 {n['unknown']}")

    # 角色分布：把"元件"和"零件"分开数（注释行数 vs 元素级行数的差别就在这里）
    roles = {}
    for r in results:
        roles[r["role"]] = roles.get(r["role"], 0) + 1
    print(f"角色分布：{' ｜ '.join(f'{k} {v}' for k, v in sorted(roles.items(), key=lambda x: -x[1]))}")

    # 自反性检查：输出的 SO 标准名必须能被自己再吃回去，否则报告不可信
    produced = [r["so_name"] for r in results if r["so_name"]]
    bad = [s for s in set(produced) if s not in alias_to_so]
    if not produced:
        # 没有输出时说"全部可回查 ✓"是空话，等于把真空当成通过
        print("自反性检查：本次没有输出 SO 名，跳过")
    else:
        print(f"自反性检查：{len(set(produced))} 个输出 SO 名，"
              + ("全部可回查 ✓" if not bad else f"{len(bad)} 个无法回查 ✗ -> {bad[:5]}"))

    detail = [r for r in results if r["status"] != "matched"]

    if args.verbose:
        # 为什么要"原始名 -> 命中键 -> SO 名"三段：
        # 只看到 "DNA/hAt -> hAT_TIR_transposon" 时，会以为表里有 DNA/hAt；
        # 看到命中键是 DNA/hAT，才知道差在一个字母的大小写。
        # 这一列是后续判断"要不要修上游库"的唯一依据，所以必须显示出来。
        print("\n" + "=" * 70)
        print("逐条命中明细（--verbose）：原始名 -> 库内实际写法 -> SO 标准名")
        print("=" * 70)
        if detail:
            print("  " + pad("路径", 14) + pad("原始名", 26) + pad("命中键", 26)
                  + "SO 标准名 (SO ID)")
            print("  " + "-" * 68)
            for r in detail:
                label = PATH_INFO[r["status"]][0]
                hit = r["hit"] or "—"
                target = (f"{r['so_name']} ({r['so_id']})" if r["so_name"]
                          else "建议：" + (r["suggestion"] or "无相似名"))
                print("  " + pad(label, 14) + pad(clip(r["name"], 25), 26)
                      + pad(clip(hit, 25), 26) + target)
            print("  " + "-" * 68)
        for st, (label, why) in PATH_INFO.items():
            if n[st]:
                print(f"  {pad(label, 14)} {n[st]:>2} 条 —— {why}")
        if not detail:
            print("\n  全部走路径1（原始名＝命中键），没有需要归一的名称。")
    else:
        for r in detail:
            print(f"  [{r['status']:<19}] {r['name']:<30} -> {r['so_name'] or r['suggestion']}")

    # 可选：写 CSV 报告
    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(["输入名", "格式", "状态", "命中键", "SO 标准名", "SO ID",
                        "类别", "角色", "超家族", "建议"])
            for r in results:
                w.writerow([r["name"], r["format"], r["status"], r["hit"], r["so_name"], r["so_id"],
                            r["order"], r["role"], r["superfamily"], r["suggestion"]])
        print(f"报告已写入 {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
