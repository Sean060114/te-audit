#!/usr/bin/env python3
"""
te_check — EDTA 注释输出的完整性事实检查器（模块 1）

不做重新注释，只回答一个问题：**这份注释的哪些结论缺少证据支撑，或证据自相矛盾。**

六项检查：
  C1 字段可得性      按类别统计哪些质控字段有值；整类缺字段的，其结论无法自校验
  C2 结构长度异常    内部区过短的"完整元件"（更可能是 MITE 或假阳性）
  C3 反向重复质量    末端反向重复的匹配率、臂长、是否完美；短臂需与随机期望对照
  C4 注释↔分类器覆盖 注释里存在的元件，是否都进入了分类器的输入表
  C5 库层面跨类      分类表里被判为其他类别的条目（建库污染或命名错位）
  C6 一致性与域证据  LTR 一致性偏高但无蛋白域证据（年龄推断可能失真）

用法：
    python te_check.py --edta-dir /path/to/edta_run
    python te_check.py --edta-dir /path/to/edta_run --json report.json
    python te_check.py --edta-dir ... --strict     # 有严重问题则退出码 1

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

# 这个 SO 类型本身就是"LTR 但超家族未知"。TEsorter 无法为它定超家族，
# 因此它不进分类表属预期行为，不能与已定超家族的缺失混为一谈。
# 水稻 218 个 LTR 里 20 个属此类、且 100% 未分类；若混算，会让总覆盖率
# 从 94.4% 掉到 85.8%，把 11 个真正的静默丢弃稀释掉。
UNASSIGNED_LTR = "LTR_retrotransposon"

MIN_INTERNAL_BP = 100       # C2 阈值：内部区低于此值视为可疑
MIN_ARM_BP = 10             # C3：短于此的臂不做统计
SHORT_ARM_WARN = 20         # C3：臂长短于此值提示"分辨率粗"


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


def parse_intact(path):
    """返回元素级记录列表（排除 LTR/TSD/repeat_region 等子特征）。"""
    recs = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ELEMENT_TYPES:
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


# ---------------------------------------------------------------- 报告

def main():
    ap = argparse.ArgumentParser(description="EDTA 注释完整性事实检查器")
    ap.add_argument("--edta-dir", required=True, help="EDTA 输出目录（含 *.EDTA.intact.gff3）")
    ap.add_argument("--json", help="把结果写成 JSON")
    ap.add_argument("--prefix", help="同目录有多份运行结果时，按文件名前缀指定目标（如 col1m.fa）")
    ap.add_argument("--strict", action="store_true", help="存在严重问题时以退出码 1 结束")
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

    recs = parse_intact(files["intact"])
    cls_rows = []
    for p in files["cls"]:
        cls_rows.extend(parse_cls(p))
    print(f"\n元素级记录: {len(recs)}   分类器条目: {len(cls_rows)}")

    report = {"dir": args.edta_dir, "n_elements": len(recs), "n_classified": len(cls_rows)}
    severe = 0

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

    # ---- 汇总
    print("\n" + "=" * 78)
    print(f"汇总：严重问题 {severe} 项（C2 长度异常 + C4 未分类元件）")
    print("=" * 78)
    report["severe_total"] = severe

    if args.json:
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print(f"[写出] {args.json}")

    if args.strict and severe:
        sys.exit(1)


if __name__ == "__main__":
    main()
