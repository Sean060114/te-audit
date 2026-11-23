#!/usr/bin/env python3
"""
te_check — EDTA 注释输出的完整性事实检查器（模块 1）

不做重新注释，只回答一个问题：**这份注释的哪些结论缺少证据支撑，或证据自相矛盾。**

文件级六项（只需 EDTA 输出目录）：
  C1 字段可得性      按类别统计哪些质控字段有值；整类缺字段的，其结论无法自校验
  C2 结构长度异常    内部区过短的"完整元件"（更可能是 MITE 或假阳性）
  C3 反向重复质量    末端反向重复的匹配率、臂长、是否完美；短臂需与随机期望对照
  C4 注释↔分类器覆盖 注释里存在的元件，是否都进入了分类器的输入表
  C5 库层面跨类      分类表里被判为其他类别的条目（建库污染或命名错位）
  C6 一致性与域证据  LTR 一致性偏高但无蛋白域证据（年龄推断可能失真）

结构与语义七项（C7-C11 另需 --fasta；C12/C13 只需输出目录）：
  C7 两套坐标        LTR 的 repeat_region（含 TSD）与分类行（不含）差值是否等于 TSD 长度
  C8 motif 位置      motif 前后 2 nt 是否在 LTR 末端（负链需 rc 且交叉）
  C9 TIR 坐标        TIR 的 TSD 是否在元件**外侧**（与 LTR 约定相反）
  C10 identity 虚高  TIR 的 identity= 与两臂实测相似度差多少（严格对位，不滑移）
  C11 长臂内侧外推   长臂内侧 1/3 的错配率是否显著高于外侧 2/3
  C12 近重复记录     起止相差 ≤2 bp 的同类别记录（同一元件被输出两次）
  C13 家族标签分散   同一 Name 在不同位点拿到不同超家族；分真冲突/含未定/仅自主性

C1-C6 查的是"文件内部矛盾"，C7-C13 查的是"**你会把哪个数字读反**"。后者每条都
有明确的失效场景：用错一套坐标、用错链方向、把字段当实测、只看总相似度不看位置
分布。这些错误都不会让程序报错，只会让结论整体偏掉。

用法：
    python te_check.py --edta-dir /path/to/edta_run
    python te_check.py --edta-dir /path/to/edta_run --json report.json
    python te_check.py --edta-dir ... --strict     # 有严重问题则退出码 1
    python te_check.py --edta-dir ... --fasta genome.fa [more.fa]   # 启用 C7-C11

只依赖标准库。
"""
import argparse
import glob
import json
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

ELEMENT_TYPES = {
    "CACTA_TIR_transposon", "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon",
    "LTR_retrotransposon", "Mutator_TIR_transposon", "PIF_Harbinger_TIR_transposon",
    "Tc1_Mariner_TIR_transposon", "hAT_TIR_transposon", "helitron",
}
TIR_TYPES = {t for t in ELEMENT_TYPES if "TIR" in t}
LTR_TYPES = {t for t in ELEMENT_TYPES if "LTR" in t}
HELITRON = {"helitron"}

# intact.gff3 里"是零件的行"，不是元件本身：解析时按设计跳过，不算漏检。
# 与"未建模的元件类型"分开——后者要计数并报出来，否则就是静默丢弃。
SUBFEATURE_TYPES = {
    "repeat_region",            # 元件容器，坐标含 TSD；另挂在 pstart/pend 上
    "long_terminal_repeat",     # LTR 本身
    "target_site_duplication",  # TSD 本身
    "primer_binding_site",
    "repeat_fragment",          # 碎片注释，不是 intact 元件
}

# 这个 SO 类型本身就是"LTR 但超家族未知"。TEsorter 无法为它定超家族，
# 因此它不进分类表属预期行为，不能与已定超家族的缺失混为一谈。
# 水稻 218 个 LTR 里 20 个属此类、且 100% 未分类；若混算，会让总覆盖率
# 从 94.4% 掉到 85.8%，把 11 个真正的静默丢弃稀释掉。
UNASSIGNED_LTR = "LTR_retrotransposon"

MIN_INTERNAL_BP = 100       # C2 阈值：内部区低于此值视为可疑
MIN_ARM_BP = 10             # C3：短于此的臂不做统计
SHORT_ARM_WARN = 20         # C3：臂长短于此值提示"分辨率粗"

# C7-C13 阈值（C7-C11 需要 --fasta）
TIR_INFLATION_TOL = 0.20    # C10：字段 identity 高出实测多少算虚高
TIR_LONG_ARM = 100          # C11：多长的臂算"长臂"
TIR_INNER_MISMATCH = 0.30   # C11：内侧 1/3 错配率超过此值算异常
TIR_OUTER_GAP = 0.20        # C11：内侧比外侧高出多少才算"边界被外推"
DUP_TOL = 2                 # C12：起止相差多少 bp 算近重复记录
FAMILY_MIN_OCC = 3          # C13：家族至少几个位点才检查标签一致性


def rc(s):
    return s.translate(COMP)[::-1]


def attrs(a):
    d = {}
    for kv in a.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


# ---------------------------------------------------------------- 解析

def find_files(edta_dir, prefix=None):
    """发现 EDTA 输出文件。

    prefix 用于同一目录下存在多份运行结果时（例如把多个基因组的输出放在一起），
    按文件名前缀挑出目标那一份。不指定且存在多份时，会警告并取文件名最靠前的一份
    ——静默取错会让报告张冠李戴，所以必须出声。
    """
    out = {}
    for pat, key in (
        ("*.EDTA.intact.gff3", "intact"),
        ("*.EDTA.TEanno.gff3", "anno"),
        ("*.EDTA.TEanno.sum", "sum"),
        ("*.EDTA.TElib.fa", "tElib"),
    ):
        hits = sorted(glob.glob(os.path.join(edta_dir, pat)))
        if prefix:
            hits = [h for h in hits if os.path.basename(h).startswith(prefix)]
        elif len(hits) > 1 and key == "intact":
            print(f"警告：目录中有 {len(hits)} 份 intact.gff3，仅分析 "
                  f"{os.path.basename(hits[0])}；如需指定请用 --prefix", file=sys.stderr)
        out[key] = hits[0] if hits else None
    cls = sorted(glob.glob(os.path.join(edta_dir, "*.EDTA.raw", "*", "*.rexdb.cls.tsv")))
    if prefix:
        cls = [c for c in cls if os.path.basename(os.path.dirname(os.path.dirname(c))).startswith(prefix)]
    out["cls"] = cls
    return out


def parse_intact(path, unmodeled=None):
    """返回元素级记录列表（排除 LTR/TSD/repeat_region 等子特征）。

    `unmodeled` 传入一个 dict 时，会把"既不是子特征、也不在 ELEMENT_TYPES 里"的
    类型计数写进去（类型 -> 条数）。这类行没有判据可用，只能跳过；但跳过必须被
    报出来——否则 EDTA 哪天多输出一种合法 SO 类型，元件会无声消失，这正是本
    工具在别处反复抓的失效模式。传 None 时行为与旧版完全一致。

    **同时收集 repeat_region 行的坐标**，挂到 pstart/pend。原因：同一元件在这里
    有**两套坐标**，差值恰好等于 TSD 长度（水稻 410/420 条为 5,5）。repeat_region
    含两侧 TSD，分类行（Gypsy_LTR_retrotransposon 等）不含。于是：
      - TSD 只能在父行上查（含），
      - motif 只能在分类行上查（不含）。
    用错一套的后果是 0% 与 1%，而两个结果都不报错——只会让人以为数据有问题。
    """
    recs = []
    parents = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "repeat_region":
                p = attrs(f[8])
                parents[p.get("ID", "")] = (int(f[3]), int(f[4]))
                continue
            if f[2] not in ELEMENT_TYPES:
                if unmodeled is not None and f[2] not in SUBFEATURE_TYPES:
                    unmodeled[f[2]] = unmodeled.get(f[2], 0) + 1
                continue
            d = attrs(f[8])
            start, end = int(f[3]), int(f[4])
            r = dict(
                chrom=f[0],
                type=f[2], start=start, end=end, length=end - start + 1,
                strand=f[6], cls=d.get("classification", ""),
                identity=None, ltr_identity=None,
                tsd=None, tsd_pct=None, tir=("", ""),
                ir_ident=None, ir_shift=None, armlen=0,
                motif=d.get("motif"), strand_ev=d.get("strand_evidence"),
                method=d.get("method"),
                # C13 需要家族名。Name 是**家族**不是元件 ID（同一家族多位点共用），
                # 个别还是反写坐标串（OsChr1:19588987..19579601），用前必须过滤。
                sid=d.get("ID", ""), name=d.get("Name", ""),
                parent=d.get("Parent", ""),
            )
            v = d.get("identity")
            if v and v != "NA":
                try:
                    r["identity"] = float(v)
                except ValueError:
                    pass
            v = d.get("ltr_identity")
            if v:
                try:
                    r["ltr_identity"] = float(v)
                except ValueError:
                    pass
            tsd = d.get("TSD") or d.get("tsd")
            if tsd:
                parts = tsd.split("_")
                r["tsd"] = parts[0]
                if len(parts) >= 3:
                    try:
                        r["tsd_pct"] = float(parts[2])
                    except ValueError:
                        pass
            tir = d.get("TIR")
            if tir and "_" in tir:
                p = tir.split("_")
                if len(p) == 2:
                    r["tir"] = (p[0], p[1])
            recs.append(r)
    for r in recs:
        pc = parents.get(r.get("parent", ""))
        if pc:
            r["pstart"], r["pend"] = pc
    return recs


def ir_identity(left, right, max_shift=10):
    """右臂对 revcomp(左臂) 滑动取最佳匹配率。"""
    if not left or not right:
        return None, None
    target = rc(left)
    n = len(target)
    best, best_shift = -1.0, None
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            a, b = target[shift:], right[:n - shift] if shift < n else ""
        else:
            a, b = target[:len(right) + shift], right[-shift:]
        m = min(len(a), len(b))
        if m < MIN_ARM_BP:
            continue
        ident = sum(1 for i in range(m) if a[i] == b[i]) / m
        if ident > best:
            best, best_shift = ident, shift
    return (best if best >= 0 else None), best_shift


def parse_cls(path):
    """解析 *.rexdb.cls.tsv，返回条目列表。"""
    rows = []
    lib = os.path.basename(path).split(".")[0:1]
    libname = os.path.basename(os.path.dirname(path))
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            m = re.match(r"(.+?):(\d+)\.\.(\d+)#(.+)", f[0])
            if not m:
                continue
            a, b = int(m.group(2)), int(m.group(3))
            rows.append(dict(
                seqid=m.group(1), start=min(a, b), end=max(a, b),
                orig_cls=m.group(4), order=f[1], superfamily=f[2],
                clade=f[3], complete=f[4], strand=f[5], domains=f[6],
                lib=libname,
            ))
    return rows


# ---------------------------------------------------------------- 检查

def check_field_availability(recs):
    """C1：每类元件有哪些质控字段真正可用。"""
    out = []
    by = defaultdict(list)
    for r in recs:
        by[r["type"]].append(r)
    for t in sorted(by, key=lambda x: -len(by[x])):
        g = by[t]
        out.append(dict(
            type=t, n=len(g),
            identity=sum(1 for r in g if r["identity"] is not None or r["ltr_identity"] is not None),
            tsd=sum(1 for r in g if r["tsd"]),
            tir=sum(1 for r in g if r["tir"][0]),
            motif=sum(1 for r in g if r["motif"]),
            strand_ev=sum(1 for r in g if r["strand_ev"]),
            median_len=statistics.median(r["length"] for r in g),
        ))
    return out


def check_short_elements(recs):
    """C2：内部区过短。只对 TIR 类有意义（臂长已知）。"""
    hits = []
    for r in recs:
        if r["type"] not in TIR_TYPES or not r["tir"][0]:
            continue
        internal = r["length"] - 2 * len(r["tir"][0])
        r["internal"] = internal
        if internal < MIN_INTERNAL_BP:
            hits.append(dict(r, internal=internal))
    return hits


def check_ir_quality(recs):
    """C3：反向重复匹配率、臂长、完美性。"""
    rows = []
    for r in recs:
        if r["type"] not in TIR_TYPES or not r["tir"][0]:
            continue
        ident, shift = ir_identity(r["tir"][0], r["tir"][1])
        r["ir_ident"], r["ir_shift"] = ident, shift
        r["armlen"] = len(r["tir"][0])
        if ident is not None:
            rows.append(r)
    return rows


def check_coverage(recs, cls_rows):
    """C4：注释中存在的元件，是否都进入了某个分类表。

    只有 LTR 类可以要求全覆盖：EDTA 的 LTR intact 元件全部经过 TEsorter，
    缺失即意味着被静默丢弃。TIR / helitron 的注释并不依赖分类器（identity
    多来自库同源匹配），分类表只覆盖建库时的一个子集，故仅作参考。
    """
    # 按染色体分组：多序列输入（如 Chr1+Chr12）时，不同染色体可能有相同坐标，
    # 只用 (start,end) 配对会把 A 染色体的分类结果当成 B 染色体的，静默错配。
    by_chrom = defaultdict(list)
    for c in cls_rows:
        by_chrom[c.get("seqid", "")].append((c["start"], c["end"]))
    every = [iv for v in by_chrom.values() for iv in v]
    tol = 5

    def covered(r):
        # 记录没有染色体信息时退化为全表匹配（保持旧行为，也兼容手工构造的记录）
        ch = r.get("chrom")
        pool = by_chrom.get(ch, ()) if ch else every
        return any(abs(s - r["start"]) <= tol and abs(e - r["end"]) <= tol
                   for s, e in pool)

    per_class = defaultdict(lambda: [0, 0])   # type -> [covered, total]
    missing_ltr, missing_other = [], []
    for r in recs:
        per_class[r["type"]][1] += 1
        if covered(r):
            per_class[r["type"]][0] += 1
        elif r["type"] in LTR_TYPES:
            missing_ltr.append(r)
        else:
            missing_other.append(r)
    return per_class, missing_ltr, missing_other


def check_cross_class(cls_rows):
    """C5：分类表里被判为其他类别的条目。"""
    issues = []
    for c in cls_rows:
        lib = c["lib"].lower()
        expect = None
        if "tir" in lib:
            expect = {"TIR"}
        elif "ltr" in lib:
            expect = {"LTR"}
        elif "helitron" in lib:
            expect = {"Helitron", "DNA"}      # TEsorter 对 helitron 输出 order=Helitron
        elif "line" in lib:
            expect = {"LINE"}
        elif "sine" in lib:
            expect = {"SINE"}
        else:
            expect = None
        if expect and c["order"] not in expect:
            issues.append(c)
    return issues


def check_ltr_conflict(recs, cls_rows):
    """C6：LTR 一致性偏高，但缺结构域证据。

    域证据有两个独立来源，任一存在即算有证据：
      1. 分类表（*.rexdb.cls.tsv）的 domains 列 —— TEsorter 的结果
      2. intact.gff3 的 strand_evidence 字段 —— LTR_retriever 自己的判定依据，
         其中 "domain" 表示检测到内部蛋白域

    只认第 1 个来源会误报：一个元件可以有域、却没进入分类器（C4 已单独报告）。
    水稻 880801..891653 就是此情形（strand_evidence=domain,PBS,PPT，
    ltr_identity=0.9976，但不在任何 cls.tsv 中）。
    """
    cls_by_coord = {(c.get("seqid", ""), c["start"], c["end"]): c for c in cls_rows}
    cls_by_pos = {(c["start"], c["end"]): c for c in cls_rows}   # 无染色体时的退化路径
    out = []
    for r in recs:
        if r["type"] not in LTR_TYPES:
            continue
        if r["ltr_identity"] is None or r["ltr_identity"] < 0.95:
            continue
        ch = r.get("chrom")
        c = (cls_by_coord.get((ch, r["start"], r["end"])) if ch
             else cls_by_pos.get((r["start"], r["end"])))
        cls_domain = bool(c and c["domains"] and c["domains"] != "NA")
        gff_domain = "domain" in (r.get("strand_ev") or "").lower()
        if not (cls_domain or gff_domain):
            out.append(dict(element=r, cls=c))
    return out


# ------------------------------------------------- C7-C13 字段语义与结构自证

def read_fasta(paths):
    """读 FASTA，返回 {染色体名: 大写序列}。

    不建 Os 前缀别名：同一序列挂两个键会在遍历时被算两次，凡是靠染色体名做的
    自身排除都会失效（pack_mule_test.py 曾因此把基线虚高到 100%）。gff 与 FASTA
    的命名差异由 _chrom_key 就地解析，不改表。
    """
    seq = {}
    for p in paths:
        name, buf = None, []
        with open(p, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        seq[name] = "".join(buf).upper()
                    name = line[1:].split()[0]
                    buf = []
                elif name:
                    buf.append(line.strip())
        if name:
            seq[name] = "".join(buf).upper()
    return seq


def _chrom_key(chrom, seq):
    """把 gff 的染色体名对到 FASTA 的键上（OsChr1 ↔ Chr1）。对不上返回 None。"""
    if chrom in seq:
        return chrom
    alt = chrom[2:] if chrom.startswith("Os") else "Os" + chrom
    return alt if alt in seq else None


def check_ltr_two_coords(recs, seq):
    """C7：LTR 元件的**两套坐标**是否自洽。

    同一元件在 intact.gff3 里出现两次，坐标不同：
      repeat_region 行   106517..118085  ← 含两侧 TSD
      分类行（Gypsy…）   106522..118080  ← 不含 TSD
    差值应当恰好等于 len(tsd)。水稻 410/420 条是 (5,5) 对 5 bp TSD，自洽；余下 10
    条里 5 条是 4/4 对 4 bp（也自洽），5 条不自洽（左端多差 2–3 bp）。

    为什么必须查：两套坐标来自**同一个文件**，没人会怀疑它们不一致；而只要下游
    一边用 repeat_region、另一边用分类行做边界比较，就会整体错位一个 TSD 长度。

    提供 FASTA 时再验一步：父行最左/最右 L bp 是否等于 tsd（负链要先 rc——tsd
    按元件方向书写）。
    """
    rows = []
    for r in recs:
        if r["type"] not in LTR_TYPES or not r.get("pstart"):
            continue
        L = len(r["tsd"]) if r["tsd"] else None
        d5, d3 = r["start"] - r["pstart"], r["pend"] - r["end"]
        rec = dict(r, d5=d5, d3=d3,
                   coord_ok=(L is not None and d5 == L and d3 == L))
        if seq and L:
            k = _chrom_key(r["chrom"], seq)
            if k:
                left = seq[k][r["pstart"] - 1:r["pstart"] - 1 + L]
                right = seq[k][r["pend"] - L:r["pend"]]
                if r["strand"] == "-":
                    left, right = rc(left), rc(right)
                rec["seq5"] = (left == r["tsd"])
                rec["seq3"] = (right == r["tsd"])
        rows.append(rec)
    return rows


def check_ltr_motif(recs, seq):
    """C8：motif= 的前 2 nt / 后 2 nt 是否落在 LTR 末端，按**元件方向**。

    反直觉之处：两翼 LTR 的**序列**按基因组正向书写（所以相似度不取 rc），
    但 tsd= 与 motif= 两个**字段**按元件方向书写（负链要 rc）。同一份文件里
    两种方向并存，混用会让一侧永远差 2 nt，而总体通过率可能高到看不出问题
    ——水稻 420 条里有 406 条的正读与 rc 读同时成立，把 5 条负链的错读掩盖了。
    """
    rows = []
    for r in recs:
        if r["type"] not in LTR_TYPES or not r["motif"] or len(r["motif"]) != 4:
            continue
        if not r["tsd"]:
            continue
        k = _chrom_key(r["chrom"], seq)
        if not k:
            continue
        s, e = r["start"] - 1, r["end"]
        if e - s < 40:
            continue
        m = r["motif"]
        # 用**分类行**坐标：它不含 TSD，所以元件两端就是两侧 LTR 的末端。
        # 若改用 repeat_region 坐标（含 TSD），5'LTR 起点要再 +len(tsd)。这里
        # 多算一个 L 会让通过率从 98.8% 掉到 1.0%，而结果只是"看起来数据坏了"。
        h = seq[k][s:s + 2]              # 5'LTR 最靠内 2 nt
        t = seq[k][e - 2:e]              # 3'LTR 最靠内 2 nt
        if r["strand"] == "-":
            # 基因组小坐标端是负链元件的 3' 端，所以要**交叉**
            h, t = rc(t), rc(h)
        rows.append(dict(r, motif_ok=(h == m[:2] and t == m[2:]),
                         head=h, tail=t))
    return rows


def check_tir_tsd_outside(recs, seq):
    """C9：TIR 元件的坐标是否**不含** TSD（与 LTR 相反）。

    这不是笔误。TIR-Learner 输出的元件两端就是臂本身，TSD 在臂的外侧；EDTA 的
    LTR 则把 TSD 算进元件，两者恰好差一个 TSD 长度（水稻实测恒差 5 bp）。按
    LTR 的方式去 TIR 上找 TSD，1,824 条会全部失败——反之亦然。跨工具比较时
    这是最先踩、也最容易被当成"数据错了"的一处。
    """
    rows = []
    for r in recs:
        if r["type"] not in TIR_TYPES or not r["tsd"]:
            continue
        k = _chrom_key(r["chrom"], seq)
        if not k:
            continue
        L = len(r["tsd"])
        s, e = r["start"] - 1, r["end"]
        if s - L < 0 or e + L > len(seq[k]):
            continue
        got5, got3 = seq[k][s - L:s], seq[k][e:e + L]
        rows.append(dict(r, out5=(got5 == r["tsd"]), out3=(got3 == r["tsd"]),
                         got5=got5, got3=got3))
    return rows


def check_tir_identity_inflation(recs):
    """C10：identity= 与实测反向重复相似度差得太多。

    这个字段**不是** TIR 的两臂相似度。EDTA 把 bed 第 6 列原样写进来（TIR-Learner
    的 score），实测部分元件的字段值比两臂真实相似度高出 0.20 以上。

    对齐方式与 C3 刻意不同：**这里不滑移**，逐位对位。
      C3 要判"是不是反向重复"，臂边界可能差几 bp，滑移 ±10 合理；
      C10 要判"字段是否虚高"，滑移会把短臂相似度抬高、正好把虚高掩盖掉。
    同一批数据、两个判据、两种对齐，是刻意的，不是不一致。
    """
    rows = []
    for r in recs:
        if r["type"] not in TIR_TYPES or r["identity"] is None or not r["tir"][0]:
            continue
        a, b = r["tir"]
        n = min(len(a), len(b))
        if n < MIN_ARM_BP:
            continue
        rb = rc(b)
        miss = sum(1 for i in range(n) if a[i] != rb[i])
        meas = 1 - miss / n
        rows.append(dict(r, measured=meas, inflation=r["identity"] - meas))
    return rows


def check_tir_arm_extrapolation(recs):
    """C11：长臂的**内侧** 1/3 错配率是否显著高于外侧 2/3。

    真臂的证据在外端，内侧边界是流程为了凑长度往外推出来的。实测：外侧错配
    0.09、内端 0.42。只看总相似度看不出这件事——1.5 kb 上散 150 个错配仍然
    "84% 相似"，但那些错配全挤在内侧，说明内侧根本不是重复序列。

    对位方式：a[0] 是左臂最外端，rc(b)[0] 是右臂最外端的互补，两者同向。
    """
    rows = []
    for r in recs:
        if r["type"] not in TIR_TYPES or not r["tir"][0]:
            continue
        a, b = r["tir"]
        L = min(len(a), len(b))
        if L < TIR_LONG_ARM:
            continue
        rb = rc(b)[:L]
        cut = (2 * L) // 3
        m_out = sum(1 for i in range(cut) if a[i] != rb[i]) / cut
        m_in = sum(1 for i in range(cut, L) if a[i] != rb[i]) / (L - cut)
        if m_in >= TIR_INNER_MISMATCH and m_in - m_out >= TIR_OUTER_GAP:
            rows.append(dict(r, arm=L, mismatch_outer=m_out, mismatch_inner=m_in))
    return rows


def check_near_duplicates(recs):
    """C12：起止相差 ≤DUP_TOL bp 的近重复记录。

    同一元件被写成两条（TSD 计入与否被当成两个元件）会让"有几个元件"整体
    虚高，所有按元件计数的比例随之偏。排序后相邻比较即可，不需要 O(n²)。
    """
    by = defaultdict(list)
    for r in recs:
        by[r["chrom"]].append(r)
    pairs = []
    for ch, v in by.items():
        v.sort(key=lambda x: (x["start"], x["end"]))
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                if v[j]["start"] - v[i]["start"] > DUP_TOL:
                    break
                if (abs(v[j]["end"] - v[i]["end"]) <= DUP_TOL
                        and v[i]["type"] == v[j]["type"]):
                    pairs.append((v[i], v[j]))
    return pairs


COORD_NAME = re.compile(r"^[\w.]+:\d+\.\.\d+$")


def _superfamily(cls):
    """取 classification 的超家族段。MITE/DTT 与 DNA/DTT 的超家族都是 DTT。"""
    if not cls:
        return ""
    parts = cls.split("/")
    return parts[-1] if len(parts) >= 2 else cls


def check_family_label_spread(recs):
    """C13：同一 Name（家族）在不同位点被标成不同**超家族**。

    三层必须分开，否则会误报成一片：
      真冲突   ≥2 个已定超家族（DTT vs DTH）—— 分类确实不稳
      含未定   1 个已定 + unknown —— unknown 是"没定"，不是"定成别的"。与 C4 的
                LTR 分层同理：把它算进冲突，会把预期的未决当成错误
      仅自主性 MITE/DTT vs DNA/DTT —— 超家族相同，差的只是自主性标签。不算冲突，
                但会直接影响"按 MITE 计数"的下游，所以单列

    Name 是家族名不是元件 ID；必须排除 `OsChr1:19588987..19579601` 这类反写坐标串，
    那是缺家族名时的占位，不排除会把每条都算成独立家族。
    """
    by = defaultdict(list)
    for r in recs:
        n = r.get("name") or ""
        if not n or COORD_NAME.match(n):
            continue
        by[n].append(r)
    out = []
    for n, g in by.items():
        if len(g) < FAMILY_MIN_OCC:
            continue
        labels = Counter(r["cls"] for r in g if r["cls"])
        sfs = Counter(_superfamily(c) for c in labels.elements())
        known = {k for k in sfs if k and k.lower() != "unknown"}
        n_unk = sum(v for k, v in sfs.items() if k.lower() == "unknown")
        if len(known) >= 2:
            kind = "真冲突"
        elif len(known) == 1 and n_unk:
            kind = "含未定"
        elif len(labels) >= 2:
            kind = "仅自主性"
        else:
            continue
        out.append(dict(name=n, n=len(g), kind=kind, labels=labels, sfs=sfs))
    order = {"真冲突": 0, "含未定": 1, "仅自主性": 2}
    return sorted(out, key=lambda d: (order[d["kind"]], -d["n"]))


# ---------------------------------------------------------------- 报告

def main():
    ap = argparse.ArgumentParser(description="EDTA 注释完整性事实检查器")
    ap.add_argument("--edta-dir", required=True, help="EDTA 输出目录（含 *.EDTA.intact.gff3）")
    ap.add_argument("--json", help="把结果写成 JSON")
    ap.add_argument("--prefix", help="同目录有多份运行结果时，按文件名前缀指定目标（如 col1m.fa）")
    ap.add_argument("--strict", action="store_true", help="存在严重问题时以退出码 1 结束")
    ap.add_argument("--fasta", nargs="+", metavar="FA",
                    help="基因组 FASTA（可给多个）。提供后才做 C7-C11 的序列级复算；"
                         "C12/C13 只读文件内部信息，无需 FASTA")
    args = ap.parse_args()

    files = find_files(args.edta_dir, args.prefix)
    if not files["intact"]:
        sys.exit(f"未在 {args.edta_dir} 找到 *.EDTA.intact.gff3")

    print("=" * 78)
    print("te_check — EDTA 注释完整性检查")
    print("=" * 78)
    print(f"目录: {args.edta_dir}")
    for k in ("intact", "anno", "sum", "tElib"):
        print(f"  {k:<7}: {os.path.basename(files[k]) if files[k] else '（缺失）'}")
    print(f"  分类表 : {len(files['cls'])} 个")
    for p in files["cls"]:
        print(f"           {os.path.basename(p)}")

    unmodeled = {}
    recs = parse_intact(files["intact"], unmodeled=unmodeled)
    cls_rows = []
    for p in files["cls"]:
        cls_rows.extend(parse_cls(p))
    print(f"\n元素级记录: {len(recs)}   分类器条目: {len(cls_rows)}")
    if unmodeled:
        detail = "、".join(f"{t} {n} 条" for t, n in
                           sorted(unmodeled.items(), key=lambda kv: (-kv[1], kv[0])))
        print(f"  ⚠ 另有 {sum(unmodeled.values())} 条元件属于未建模类型，已跳过：{detail}")
        print(f"     这些类型不在 ELEMENT_TYPES 的九类里，十三项检查不覆盖它们；"
              f"要纳入就得先补判据与测试（清单见 test_te_check.UNMODELED_SO_TYPES）。")
    if args.fasta:
        # 路径写错时给一句人话，不要 traceback。公开仓库的第一个访客最容易撞这里。
        missing = [p for p in args.fasta if not os.path.exists(p)]
        if missing:
            sys.exit("找不到 FASTA 文件：" + ", ".join(missing))
    seq = read_fasta(args.fasta) if args.fasta else {}
    if args.fasta:
        print(f"基因组序列: {len(seq)} 条（{', '.join(sorted(seq))[:60]}）")

    report = {"dir": args.edta_dir, "n_elements": len(recs), "n_classified": len(cls_rows)}
    if unmodeled:
        # 只在真的出现未建模类型时写这个键，保证老数据的 JSON 逐字节不变
        report["unmodeled_types"] = dict(sorted(unmodeled.items()))
    severe = 0
    notice = 0                     # C7-C13 的提示数，不进退出码

    # ---- C1
    print("\n" + "-" * 78)
    print("C1 字段可得性（按类别）")
    print("-" * 78)
    fa = check_field_availability(recs)
    print(f"{'类别':<30}{'n':>4}{'identity':>10}{'TSD':>6}{'TIR':>6}{'motif':>7}{'strand_ev':>11}")
    for d in fa:
        print(f"{d['type']:<30}{d['n']:>4}{d['identity']:>10}{d['tsd']:>6}{d['tir']:>6}{d['motif']:>7}{d['strand_ev']:>11}")
    no_id = [d for d in fa if d["identity"] == 0]
    for d in no_id:
        print(f"  ⚠ {d['type']}：{d['n']} 个元件全部无 identity 数值，其边界无法用文件内信息自校验")
    report["C1_field_availability"] = fa
    report["C1_no_identity_classes"] = [d["type"] for d in no_id]

    # ---- C2
    print("\n" + "-" * 78)
    print(f"C2 结构长度异常（TIR 类内部区 < {MIN_INTERNAL_BP} bp）")
    print("-" * 78)
    short = check_short_elements(recs)
    if not short:
        print("  未发现")
    else:
        severe += len(short)
        print(f"{'坐标':<24}{'类型':<26}{'总长':>6}{'臂长':>5}{'内部区':>7}")
        for r in sorted(short, key=lambda x: x["internal"]):
            print(f"{r['start']}..{r['end']:<10}{r['type'][:24]:<26}{r['length']:>6}"
                  f"{len(r['tir'][0]):>5}{r['internal']:>7}")
        print(f"  合计 {len(short)} 个。内部区过短者更可能是 MITE 或假阳性，不应与完整元件同列。")
    report["C2_short_elements"] = [
        dict(start=r["start"], end=r["end"], type=r["type"], length=r["length"],
             armlen=len(r["tir"][0]), internal=r["internal"]) for r in short]

    # ---- C3
    print("\n" + "-" * 78)
    print("C3 反向重复质量（±10 bp 滑移）")
    print("-" * 78)
    ir_rows = check_ir_quality(recs)
    if ir_rows:
        idents = [r["ir_ident"] for r in ir_rows]
        arms = [r["armlen"] for r in ir_rows]
        perfect = sum(1 for i in idents if i >= 0.999)
        print(f"  元件数 {len(ir_rows)}   完美反向重复 {perfect}")
        print(f"  匹配率: 中位 {statistics.median(idents):.3f}  范围 {min(idents):.3f}–{max(idents):.3f}")
        print(f"  臂长  : 中位 {statistics.median(arms):.0f} bp  范围 {min(arms)}–{max(arms)}")
        short_arm = sum(1 for a in arms if a < SHORT_ARM_WARN)
        if short_arm:
            print(f"  ⚠ {short_arm} 个元件的臂长 < {SHORT_ARM_WARN} bp：该长度下 2–3 个错配就足以把"
                  f"匹配率压到 0.8 左右，单看匹配率会误判为'高度相似'")
        b = Counter()
        for i in idents:
            b["完美(=1)"] += i >= 0.999
            b["0.95–0.999"] += 0.95 <= i < 0.999
            b["0.85–0.95"] += 0.85 <= i < 0.95
            b["0.70–0.85"] += 0.70 <= i < 0.85
            b["<0.70"] += i < 0.70
        for k, v in b.items():
            if v:
                print(f"    {k:<12} {v:>3}")

        # 匹配率的下限不是本工具算出来的，而是上游 TIR-Learner 的 candidate 筛选阈值。
        # TIR-Learner 的 minimum similarity 为 80%，低于此值的候选根本不会出现在
        # intact.gff3 里，所以此处只能在 [0.80, 1.0] 内取值。不点明这一点，
        # "中位 0.81" 会被读成"两臂只有 81% 相似"，进而被当成物种间的保守性差异。
        floor = min(idents)
        at_floor = sum(1 for i in idents if i < floor + 0.01)
        if floor >= 0.795:
            print(f"  ⚠ 匹配率下限 {floor:.3f}，{at_floor}/{len(idents)} 个元件贴在此下限上。")
            print(f"     该下限来自上游 TIR-Learner 的 minimum similarity（80%），不是本工具的计算边界——")
            print(f"     低于 80% 的候选不会进入 intact.gff3。因此中位匹配率测的是上游筛选阈值，")
            print(f"     不能当作 TIR 保守性度量，也不能用于跨物种/跨数据集比较。")
            print(f"     臂长同理：TIR-Learner 最小 TIR 长度 10 bp，故臂长分布也在此处截断。")
    _floor = min((r["ir_ident"] for r in ir_rows), default=None)
    report["C3_ir_quality"] = {
        "n": len(ir_rows),
        "perfect": sum(1 for r in ir_rows if r["ir_ident"] >= 0.999),
        "median_ident": round(statistics.median([r["ir_ident"] for r in ir_rows]), 4) if ir_rows else None,
        "median_armlen": statistics.median([r["armlen"] for r in ir_rows]) if ir_rows else None,
        # 匹配率的经验下界与贴界元件数：用于提示该指标被上游阈值截断
        "ident_floor": round(_floor, 4) if _floor is not None else None,
        "at_floor_n": (sum(1 for r in ir_rows if r["ir_ident"] < _floor + 0.01)
                       if _floor is not None else None),
    }

    # ---- C4
    print("\n" + "-" * 78)
    print("C4 注释 ↔ 分类器覆盖")
    print("-" * 78)
    per_class, missing_ltr, missing_other = check_coverage(recs, cls_rows)
    print(f"  {'类别':<30}{'总数':>6}{'已分类':>8}{'覆盖率':>9}")
    for t in sorted(per_class, key=lambda x: -per_class[x][1]):
        c, n = per_class[t]
        flag = "  ← 应全覆盖" if t in LTR_TYPES else ""
        print(f"  {t:<30}{n:>6}{c:>8}{c / n * 100:>8.1f}%{flag}")

    # LTR 必须分层报告。LTR_retrotransposon 是"LTR 但超家族未知"的 SO 类型，
    # TEsorter 无法给它定超家族，它不进分类表是预期行为；而已定超家族（Copia/Gypsy）
    # 的缺失才是静默丢弃。两者混在一个百分比里，20 个预期内缺失会把 11 个真问题稀释掉。
    asg_all = [r for r in recs if r["type"] in LTR_TYPES and r["type"] != UNASSIGNED_LTR]
    asg_miss = [r for r in missing_ltr if r["type"] != UNASSIGNED_LTR]
    una_all = [r for r in recs if r["type"] == UNASSIGNED_LTR]
    una_miss = [r for r in missing_ltr if r["type"] == UNASSIGNED_LTR]
    if asg_all or una_all:
        print(f"\n  LTR 分层（{len(asg_all) + len(una_all)} 个 LTR 元件）")
        for name, tot, miss, note in (
            ("已定超家族", asg_all, asg_miss, "← 缺失即静默丢弃，须逐个核查"),
            ("超家族未定", una_all, una_miss, "← 缺失属预期（TEsorter 无法定超家族）"),
        ):
            if not tot:
                continue
            print(f"    {name}  {len(tot) - len(miss)}/{len(tot)} 已分类"
                  f"（缺 {len(miss)}，{len(miss) / len(tot) * 100:.1f}%）   {note}")
        report["C4_ltr_tiers"] = {
            "assigned_superfamily": {"total": len(asg_all), "missing": len(asg_miss)},
            "unassigned_superfamily": {"total": len(una_all), "missing": len(una_miss)},
        }

    if asg_miss:
        severe += len(asg_miss)
        print(f"\n  ⚠ {len(asg_miss)} 个已定超家族的 LTR 元件未进分类表（真问题）：")
        for r in asg_miss:
            print(f"      {r['start']}..{r['end']}  {r['type']:<26} {r['cls']:<14} "
                  f"strand_ev={r['strand_ev'] or '(无)'}")
        print("  含义：这些元件结构完整、超家族已定，却从未进入分类器。")
        print("        其分类标签不是分类器给出的结果，LTR 一致性、结构域等判断可能因此失真。")
    elif asg_all:
        print("\n  已定超家族的 LTR 全部覆盖，无静默丢弃")

    if una_miss:
        print(f"\n  {len(una_miss)}/{len(una_all)} 个超家族未定的 LTR 未进分类表 —— 预期行为，"
              f"单独计数，不计入严重项。")
        for r in una_miss[:5]:
            print(f"      {r['start']}..{r['end']}  {r['cls']:<14} "
                  f"strand_ev={r['strand_ev'] or '(无)'}")
        if len(una_miss) > 5:
            print(f"      …… 其余 {len(una_miss) - 5} 个见 JSON 的 C4_not_classified")
    if missing_other:
        print(f"\n  其他类别未覆盖 {len(missing_other)} 个（TIR/helitron 的注释不依赖分类器，"
              f"仅作参考）")
    report["C4_coverage"] = {t: {"covered": v[0], "total": v[1]} for t, v in per_class.items()}
    report["C4_not_classified"] = [
        dict(start=r["start"], end=r["end"], type=r["type"], cls=r["cls"]) for r in missing_ltr]

    # ---- C5
    print("\n" + "-" * 78)
    print("C5 库层面跨类（分类表里被判为其他类别）")
    print("-" * 78)
    cross = check_cross_class(cls_rows)
    if not cross:
        print("  未发现")
    else:
        by_lib = Counter((c["lib"], c["order"]) for c in cross)
        for (lib, order), n in by_lib.most_common():
            print(f"  {lib:<12} 被判为 {order:<8} {n} 条")
        print("  含义：建库时混入了其他类别的序列。分类器判定正确、后续被排除属预期；")
        print("        但若这些条目最终仍出现在注释里，则说明过滤环节失效。")
    report["C5_cross_class"] = [
        dict(start=c["start"], end=c["end"], lib=c["lib"], order=c["order"],
             orig=c["orig_cls"], domains=c["domains"]) for c in cross]

    # ---- C6
    print("\n" + "-" * 78)
    print("C6 LTR 一致性偏高但缺结构域证据")
    print("-" * 78)
    conflicts = check_ltr_conflict(recs, cls_rows)
    if not conflicts:
        print("  未发现")
    else:
        for d in conflicts:
            r = d["element"]
            print(f"  {r['start']}..{r['end']}  {r['type'][:24]:<26} "
                  f"ltr_identity={r['ltr_identity']:.4f}  strand_ev={r['strand_ev'] or '(无)'}")
            print(f"    分类表: {'无对应条目' if d['cls'] is None else d['cls']['order'] + '/' + d['cls']['superfamily'] + ' 域=' + d['cls']['domains']}")
        print("  含义：高一致性常被读作'近期插入'，但缺域证据时该推断不可靠——")
        print("        LTR 间基因转换会重置一致性，使老元件看起来年轻。")
    report["C6_ltr_conflict"] = [
        dict(start=d["element"]["start"], end=d["element"]["end"],
             type=d["element"]["type"], ltr_identity=d["element"]["ltr_identity"],
             strand_evidence=d["element"]["strand_ev"],
             classified=d["cls"] is not None) for d in conflicts]

    # ---- C7-C13
    print("\n" + "=" * 78)
    print("C7-C13 字段语义与结构自证")
    print("=" * 78)
    if not seq:
        print("  未提供 --fasta：C7-C11（序列级五项）跳过，C12/C13 照常执行。")
        print("  这五项无法只靠文件内部信息判定——它们的判据就是序列本身。")

    # ---- C7
    if seq:
        print("\n" + "-" * 78)
        print("C7 LTR 两套坐标是否自洽（repeat_region 含 TSD，分类行不含）")
        print("-" * 78)
        c7 = check_ltr_two_coords(recs, seq)
        bad = []
        if not c7:
            print("  无可判元件（无带 repeat_region 父行的 LTR）")
        else:
            okc = [r for r in c7 if r["coord_ok"]]
            print(f"  可判 {len(c7)} 条   自洽 {len(okc)} = {len(okc) / len(c7) * 100:.1f}%"
                  f"（判据 Δstart = Δend = len(tsd)）")
            for k, v in Counter((r["d5"], r["d3"], len(r["tsd"]) if r["tsd"] else 0)
                                for r in c7).most_common(3):
                print(f"    (Δstart, Δend, len(tsd)) = {k} → {v} 条")
            print(f"  父行序列验 TSD：5′ {sum(1 for r in c7 if r.get('seq5'))}/{len(c7)}"
                  f"   3′ {sum(1 for r in c7 if r.get('seq3'))}/{len(c7)}")
            bad = [r for r in c7 if not r["coord_ok"]]
            notice += len(bad)
            for r in bad[:5]:
                print(f"    {r['chrom']}:{r['start']}..{r['end']}  父行 {r['pstart']}..{r['pend']}"
                      f"  Δ={r['d5']},{r['d3']}  tsd={r['tsd']}")
            if len(bad) > 5:
                print(f"    …… 其余 {len(bad) - 5} 条见 JSON 的 C7")
            print("  含义：两套坐标来自同一个文件，差值应恰为一侧 TSD 长度。")
            print("        不自洽时，混用两套坐标做边界比较会整体错位一个 TSD，而不会报错。")
        report["C7_ltr_two_coords"] = dict(
            checked=len(c7), ok=sum(1 for r in c7 if r["coord_ok"]),
            seq5=sum(1 for r in c7 if r.get("seq5")),
            seq3=sum(1 for r in c7 if r.get("seq3")),
            bad=[dict(chrom=r["chrom"], start=r["start"], end=r["end"],
                      parent=[r["pstart"], r["pend"]], d5=r["d5"], d3=r["d3"],
                      tsd=r["tsd"]) for r in bad[:50]])

        # ---- C8
        print("\n" + "-" * 78)
        print("C8 LTR motif 位置（负链需 rc 且交叉）")
        print("-" * 78)
        c8 = check_ltr_motif(recs, seq)
        bad = []          # 空结果时下面的 report 也要能构造
        if not c8:
            print("  无可判元件（无 length=4 的 motif 字段）")
        else:
            ok = [r for r in c8 if r["motif_ok"]]
            print(f"  可判 {len(c8)} 条   通过 {len(ok)} = {len(ok) / len(c8) * 100:.1f}%")
            neg = [r for r in c8 if r["strand"] == "-"]
            if neg:
                print(f"    其中负链 {len(neg)} 条，通过 {sum(r['motif_ok'] for r in neg)}")
            bad = [r for r in c8 if not r["motif_ok"]]
            notice += len(bad)
            for r in bad[:5]:
                print(f"    {r['chrom']}:{r['start']}..{r['end']} {r['strand']}  "
                      f"motif={r['motif']}  末端取到 {r['head']}/{r['tail']}  tsd={r['tsd']}")
            if len(bad) > 5:
                print(f"    …… 其余 {len(bad) - 5} 条见 JSON 的 C8")
            if bad:
                print("  含义：负链元件要按元件方向取 rc 并交叉（小坐标端是它的 3′ 端）。")
                print("        不这样做时总体通过率可能仍在 97% 以上——因为末端 2 nt 常相同或回文，")
                print("        会把少数负链元件的错读整个掩盖。查错要看负链子集，不看总体。")
        report["C8_ltr_motif"] = dict(
            checked=len(c8), ok=sum(1 for r in c8 if r["motif_ok"]),
            bad=[dict(chrom=r["chrom"], start=r["start"], strand=r["strand"],
                      motif=r["motif"], head=r["head"], tail=r["tail"]) for r in bad[:50]])

        # ---- C9
        print("\n" + "-" * 78)
        print("C9 TIR 坐标是否不含 TSD（TIR-Learner 约定：不含）")
        print("-" * 78)
        c9 = check_tir_tsd_outside(recs, seq)
        bad = []          # 空结果时下面的 report 也要能构造
        if not c9:
            print("  无可判元件（无带 TSD 的 TIR）")
        else:
            both = [r for r in c9 if r["out5"] and r["out3"]]
            print(f"  可判 {len(c9)} 条   左侧外一致 {sum(r['out5'] for r in c9)}"
                  f"   右侧外一致 {sum(r['out3'] for r in c9)}"
                  f"   两侧均一致 {len(both)} = {len(both) / len(c9) * 100:.1f}%")
            bad = [r for r in c9 if not (r["out5"] and r["out3"])]
            notice += len(bad)
            for r in bad[:5]:
                print(f"    {r['chrom']}:{r['start']}..{r['end']}  tsd={r['tsd']}"
                      f"  左外取到 {r['got5']}  右外取到 {r['got3']}")
            if len(bad) > 5:
                print(f"    …… 其余 {len(bad) - 5} 条见 JSON 的 C9")
            print("  含义：TIR 的元件两端就是臂，TSD 在臂外侧；LTR 则把 TSD 算进元件。")
            print("        两者差恰好一个 TSD 长度（水稻实测恒差 5 bp）。按对方约定去查会全军覆没。")
        report["C9_tir_tsd_outside"] = dict(
            checked=len(c9),
            ok5=sum(1 for r in c9 if r["out5"]), ok3=sum(1 for r in c9 if r["out3"]),
            bad=[dict(chrom=r["chrom"], start=r["start"], end=r["end"],
                      tsd=r["tsd"], got5=r["got5"], got3=r["got3"]) for r in bad[:50]])

        # ---- C10
        print("\n" + "-" * 78)
        print(f"C10 TIR identity 字段 vs 实测（差 ≥{TIR_INFLATION_TOL:.2f} 算虚高）")
        print("-" * 78)
        c10 = check_tir_identity_inflation(recs)
        if not c10:
            print("  无可判元件（无 identity 数值的 TIR）")
        else:
            hits = [r for r in c10 if r["inflation"] >= TIR_INFLATION_TOL]
            med = statistics.median(r["inflation"] for r in c10)
            print(f"  可判 {len(c10)} 条   |字段 − 实测| 中位 {med:+.3f}"
                  f"   虚高 ≥{TIR_INFLATION_TOL:.2f} 的 {len(hits)} 条 = "
                  f"{len(hits) / len(c10) * 100:.1f}%")
            notice += len(hits)
            for r in sorted(hits, key=lambda x: -x["inflation"])[:5]:
                print(f"    {r['chrom']}:{r['start']}..{r['end']}  臂 {len(r['tir'][0])}/"
                      f"{len(r['tir'][1])}  字段 {r['identity']:.3f}  实测 {r['measured']:.3f}"
                      f"  差 {r['inflation']:+.3f}")
            if len(hits) > 5:
                print(f"    …… 其余 {len(hits) - 5} 条见 JSON 的 C10")
            print("  含义：identity= 不是 TIR 的两臂相似度，是 TIR-Learner 的 score 原样写入。")
            print("        虚高的那批不能拿它当'臂有多像'用；判臂质量请看 C3 的实测值。")
        report["C10_tir_identity_inflation"] = dict(
            checked=len(c10),
            median_delta=round(statistics.median(r["inflation"] for r in c10), 4) if c10 else None,
            inflated=sum(1 for r in c10 if r["inflation"] >= TIR_INFLATION_TOL),
            worst=[dict(chrom=r["chrom"], start=r["start"], end=r["end"],
                        field=r["identity"], measured=round(r["measured"], 4),
                        delta=round(r["inflation"], 4))
                   for r in sorted(c10, key=lambda x: -x["inflation"])[:50]])

        # ---- C11
        print("\n" + "-" * 78)
        print(f"C11 TIR 长臂内侧外推（臂 ≥{TIR_LONG_ARM}，内侧错配 ≥{TIR_INNER_MISMATCH:.2f}）")
        print("-" * 78)
        c11 = check_tir_arm_extrapolation(recs)
        n_long = sum(1 for r in recs
                     if r["type"] in TIR_TYPES and r["tir"][0]
                     and min(len(r["tir"][0]), len(r["tir"][1])) >= TIR_LONG_ARM)
        if not n_long:
            print(f"  无臂长 ≥{TIR_LONG_ARM} bp 的 TIR，不适用")
        elif not c11:
            print(f"  长臂 {n_long} 条，未发现内侧错配显著高于外侧者")
        else:
            print(f"  长臂 {n_long} 条   可疑 {len(c11)} 条 = {len(c11) / n_long * 100:.1f}%")
            notice += len(c11)
            for r in sorted(c11, key=lambda x: -x["mismatch_inner"])[:5]:
                print(f"    {r['chrom']}:{r['start']}..{r['end']}  臂 {r['arm']} bp"
                      f"  外侧错配 {r['mismatch_outer']:.3f}  内侧 {r['mismatch_inner']:.3f}")
            if len(c11) > 5:
                print(f"    …… 其余 {len(c11) - 5} 条见 JSON 的 C11")
            print("  含义：真臂的证据在外端，内侧边界是流程为凑长度外推出来的。")
            print("        总相似度看不出来（1.5 kb 散 150 个错配仍'84% 相似'），要看错配的**位置分布**。")
        report["C11_tir_arm_extrapolation"] = dict(
            long_arms=n_long, suspicious=len(c11),
            rows=[dict(chrom=r["chrom"], start=r["start"], end=r["end"], arm=r["arm"],
                       outer=round(r["mismatch_outer"], 4),
                       inner=round(r["mismatch_inner"], 4)) for r in c11[:50]])

    # ---- C12
    print("\n" + "-" * 78)
    print(f"C12 近重复记录（起止相差 ≤{DUP_TOL} bp、同类别）")
    print("-" * 78)
    c12 = check_near_duplicates(recs)
    if not c12:
        print("  未发现")
    else:
        print(f"  {len(c12)} 对，涉及 {len({id(x) for p in c12 for x in p})} 条记录")
        notice += len(c12)
        for a, b in c12[:5]:
            print(f"    {a['chrom']}:{a['start']}..{a['end']} 与 "
                  f"{b['start']}..{b['end']}   {a['type']}")
        if len(c12) > 5:
            print(f"    …… 其余 {len(c12) - 5} 对见 JSON 的 C12")
        print("  含义：同一元件被输出两次（TSD 计入与否被当成两个元件）会让元件数虚高，")
        print("        所有按元件计数的比例随之偏。去重后再算任何率。")
    report["C12_near_duplicates"] = [
        dict(chrom=a["chrom"], a=[a["start"], a["end"]], b=[b["start"], b["end"]],
             type=a["type"]) for a, b in c12]

    # ---- C13
    print("\n" + "-" * 78)
    print(f"C13 同家族多标签（位点 ≥{FAMILY_MIN_OCC} 且标签 ≥2 种）")
    print("-" * 78)
    c13 = check_family_label_spread(recs)
    if not c13:
        print("  未发现")
    else:
        kinds = Counter(d["kind"] for d in c13)
        print(f"  涉及 {len(c13)} 个家族：" + "  ".join(
            f"{k} {kinds[k]}" for k in ("真冲突", "含未定", "仅自主性") if kinds[k]))
        notice += kinds["真冲突"]
        for d in c13[:8]:
            print(f"    [{d['kind']}] {d['name']:<14}{d['n']:>3} 位点  "
                  + "  ".join(f"{k}×{v}" for k, v in d["labels"].most_common(3)))
        if len(c13) > 8:
            print(f"    …… 其余 {len(c13) - 8} 个家族见 JSON 的 C13")
        print("  含义：三层要分开读——真冲突是分类不稳；含未定是预期内（unknown 不是'定成别的'）；")
        print("        仅自主性（MITE/DTT vs DNA/DTT）超家族相同，但会让按 MITE 计数的下游偏。")
    report["C13_family_label_spread"] = [
        dict(name=d["name"], n=d["n"], kind=d["kind"], labels=dict(d["labels"]))
        for d in c13]

    # ---- 汇总
    print("\n" + "=" * 78)
    print(f"汇总：严重问题 {severe} 项（C2 长度异常 + C4 未分类元件）")
    if notice:
        print(f"      结构级提示 {notice} 项（C7-C13；不改变退出码，需人工判读）")
    print("=" * 78)
    report["severe_total"] = severe
    report["note_total"] = notice

    if args.json:
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print(f"[写出] {args.json}")

    if args.strict and severe:
        sys.exit(1)


if __name__ == "__main__":
    main()
