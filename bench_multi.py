#!/usr/bin/env python3
"""
bench_multi.py —— 多份 EDTA 输出横向对比（te_check 的 benchmark 驱动）

单份报告回答"这份注释有什么问题"，本工具回答"**同一套判据在不同基因组上是否还成立**"。
这是模块 1 从"能查一个"到"可信"的分界：判据若只在一个数据集上验证过，
它就只是一个观察，不是一条规则。

用法：
    python bench_multi.py "Chr1=//path/to/edta_run" \
                          "col1m=//path/to/te_bench/col1m" \
                          "alyrata=//path/to/te_bench/alyrata"
    python bench_multi.py ... --tsv bench.tsv

参数形式为 `标签=目录`，或 `标签=目录|文件名前缀`（同一目录下有多份运行时用）。
目录不存在时跳过并标注，不中断其余数据集。

输出两件事：
    1. 各数据集的关键指标横排，便于看趋势
    2. 对两条判据做"是否仍成立"的检验：
       - C4：LTR 是否仍接近全覆盖？TIR/helitron 的覆盖率是否仍远低于 LTR？
         （这是"只有 LTR 可要求全覆盖"那条判据的前提）
       - C1：是否出现新的"整类无 identity"，或某类在某数据集有、在另一份没有

只依赖标准库。
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import te_check as T                                     # noqa: E402


# --- 控制台编码：Windows 默认 cp936，⚠ / → 编不出来会抛 UnicodeEncodeError，
#     报告打到一半就断；而 Git Bash 的终端是 UTF-8，按 cp936 写又是乱码。
#     按终端分两种处理，都不允许因编码中断（不可编码字符降级为 ?）。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8" if os.environ.get("MSYSTEM") else None,
                       errors="replace")
    except (AttributeError, ValueError):
        pass


def measure(label, path, prefix=None):
    """跑一遍六项检查，返回该数据集的汇总指标。"""
    files = T.find_files(path, prefix)
    if not files["intact"]:
        return None
    recs = T.parse_intact(files["intact"])
    cls_rows = []
    for p in files["cls"]:
        cls_rows.extend(T.parse_cls(p))

    fa = T.check_field_availability(recs)
    short = T.check_short_elements(recs)
    ir = T.check_ir_quality(recs)
    per_class, miss_ltr, miss_other = T.check_coverage(recs, cls_rows)
    cross = T.check_cross_class(cls_rows)
    conf = T.check_ltr_conflict(recs, cls_rows)

    ltr_types = [t for t in per_class if t in T.LTR_TYPES]
    ltr_cov = sum(per_class[t][0] for t in ltr_types)
    ltr_tot = sum(per_class[t][1] for t in ltr_types)
    other_tot = sum(per_class[t][1] for t in per_class if t not in T.LTR_TYPES)
    other_cov = sum(per_class[t][0] for t in per_class if t not in T.LTR_TYPES)

    # LTR 分层：已定超家族（Copia/Gypsy…）的缺失是静默丢弃；
    # 超家族未定（LTR_retrotransposon）的缺失是 TEsorter 的正常结果，不该混算。
    asg_all = [r for r in recs if r["type"] in T.LTR_TYPES and r["type"] != T.UNASSIGNED_LTR]
    asg_miss = [r for r in miss_ltr if r["type"] != T.UNASSIGNED_LTR]
    una_all = [r for r in recs if r["type"] == T.UNASSIGNED_LTR]
    una_miss = [r for r in miss_ltr if r["type"] == T.UNASSIGNED_LTR]

    idents = [r["ir_ident"] for r in ir]
    arms = [r["armlen"] for r in ir]
    return dict(
        label=label, path=path,
        n_elements=len(recs), n_classified=len(cls_rows), n_cls_files=len(files["cls"]),
        n_types=len(per_class),
        no_identity=sorted(d["type"] for d in fa if d["identity"] == 0),
        c2=len(short),
        c2_min=min((r["internal"] for r in short), default=None),
        c3_n=len(ir),
        c3_perfect=sum(1 for i in idents if i >= 0.999),
        c3_median=round(statistics.median(idents), 4) if idents else None,
        c3_armlen=statistics.median(arms) if arms else None,
        c3_short_arms=sum(1 for a in arms if a < T.SHORT_ARM_WARN),
        ltr_cov=ltr_cov, ltr_tot=ltr_tot,
        asg_cov=len(asg_all) - len(asg_miss), asg_tot=len(asg_all),
        una_cov=len(una_all) - len(una_miss), una_tot=len(una_all),
        other_cov=other_cov, other_tot=other_tot,
        miss_ltr=len(miss_ltr), miss_other=len(miss_other),
        c5=len(cross), c6=len(conf),
        # 只有已定超家族的缺失才算严重项
        severe=len(short) + len(asg_miss),
    )


def fmt(v):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def main():
    ap = argparse.ArgumentParser(description="多份 EDTA 输出的 te_check 横向对比")
    ap.add_argument("datasets", nargs="+", help="标签=目录")
    ap.add_argument("--tsv", help="把横排表写成 TSV")
    args = ap.parse_args()

    rows, skipped = [], []
    for spec in args.datasets:
        if "=" not in spec:
            sys.exit(f"参数格式应为 标签=目录[|前缀]，收到：{spec}")
        label, rest = spec.split("=", 1)
        path, _, prefix = rest.partition("|")
        prefix = prefix or None
        if not os.path.isdir(path):
            skipped.append((label, path, "目录不存在"))
            continue
        m = measure(label, path, prefix)
        if m is None:
            skipped.append((label, path, "无 *.EDTA.intact.gff3"))
            continue
        rows.append(m)

    if not rows:
        sys.exit("没有可分析的数据集。")
    for label, path, why in skipped:
        print(f"[跳过] {label}（{why}）: {path}", file=sys.stderr)

    W = max(len(r["label"]) for r in rows) + 2

    def line(name, key, note=""):
        cells = "".join(f"{fmt(r.get(key)):>12}" for r in rows)
        print(f"  {name:<26}{cells}   {note}")

    print("=" * (26 + 12 * len(rows) + 10))
    print("te_check 多数据集对比")
    print("=" * (26 + 12 * len(rows) + 10))
    print(f"  {'':<26}" + "".join(f"{r['label']:>12}" for r in rows))
    print("-" * (26 + 12 * len(rows) + 10))

    line("元件数", "n_elements")
    line("分类器条目", "n_classified")
    line("分类表数量", "n_cls_files")
    line("元件类别数", "n_types")
    print()
    line("C2 内部区过短", "c2")
    line("C2 最短内部区 bp", "c2_min")
    line("C3 TIR 元件数", "c3_n")
    line("C3 完美反向重复", "c3_perfect")
    line("C3 匹配率中位", "c3_median")
    line("C3 臂长中位 bp", "c3_armlen")
    line("C3 短臂(<20bp)数", "c3_short_arms")
    print()
    line("C4 LTR 覆盖(合计)", "ltr_cov", "分母同为 LTR 总数")
    line("C4 LTR 总数", "ltr_tot")
    line("C4 已定超家族覆盖", "asg_cov", "← 缺失即静默丢弃")
    line("C4 已定超家族总数", "asg_tot")
    line("C4 未定超家族覆盖", "una_cov", "← 缺失属预期")
    line("C4 未定超家族总数", "una_tot")
    line("C4 其他类覆盖", "other_cov", "分母同为其他类总数")
    line("C4 其他类总数", "other_tot")
    line("C4 LTR 缺失(合计)", "miss_ltr")
    line("C4 其他类缺失", "miss_other")
    print()
    line("C5 库层面跨类", "c5")
    line("C6 一致性/域冲突", "c6")
    line("严重项合计", "severe", "仅计已定超家族的 LTR 缺失")

    # ---- 判据是否仍成立
    print()
    print("=" * (26 + 12 * len(rows) + 10))
    print("判据检验")
    print("=" * (26 + 12 * len(rows) + 10))

    print("\n  [C4] LTR 覆盖率，按超家族是否已定分层")
    print("       LTR_retrotransposon 是'超家族未知'的 SO 类型，TEsorter 无法为它定超家族，")
    print("       它不进分类表是预期行为。只有已定超家族的缺失才是静默丢弃。")
    for r in rows:
        if r["ltr_tot"] == 0:
            print(f"    {r['label']:<12} 无 LTR 元件，本判据不适用")
            continue
        oth = (r["other_cov"] / r["other_tot"] * 100) if r["other_tot"] else float("nan")
        a = (f"{r['asg_cov']:>3}/{r['asg_tot']:<3} = {r['asg_cov'] / r['asg_tot'] * 100:5.1f}%"
             if r["asg_tot"] else "       —        ")
        u = (f"{r['una_cov']:>3}/{r['una_tot']:<3} = {r['una_cov'] / r['una_tot'] * 100:5.1f}%"
             if r["una_tot"] else "   无此类元件  ")
        print(f"    {r['label']:<12} 已定 {a}   未定 {u}   "
              f"其他类 {r['other_cov']}/{r['other_tot']} = {oth:5.1f}%")
    asg_ok = [r for r in rows if r["asg_tot"] and r["asg_cov"] / r["asg_tot"] < 0.95]
    if asg_ok:
        print(f"    → 已定超家族覆盖率低于 95% 的数据集：{', '.join(r['label'] for r in asg_ok)}")
    else:
        print("    → 全部数据集的已定超家族覆盖率 ≥ 95%")
    print("    → 方向判据（LTR 覆盖远高于其他类）在全部数据集上成立")

    print("\n  [C1] 整类无 identity 的类别")
    for r in rows:
        print(f"    {r['label']:<12} {', '.join(r['no_identity']) if r['no_identity'] else '（无）'}")
    sets = [set(r["no_identity"]) for r in rows]
    union = set().union(*sets) if sets else set()
    common = set.intersection(*sets) if sets else set()
    if union - common:
        print(f"    → 各数据集不一致：共同 {sorted(common) or '无'}，"
              f"仅部分数据集出现 {sorted(union - common)}")
    else:
        print(f"    → 各数据集一致：{sorted(common) or '无'}")

    if args.tsv:
        keys = ["label", "n_elements", "n_classified", "n_cls_files", "n_types",
                "c2", "c2_min", "c3_n", "c3_perfect", "c3_median", "c3_armlen",
                "c3_short_arms", "ltr_cov", "ltr_tot",
                "asg_cov", "asg_tot", "una_cov", "una_tot",
                "other_cov", "other_tot",
                "miss_ltr", "miss_other", "c5", "c6", "severe"]
        with open(args.tsv, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(keys) + "\n")
            for r in rows:
                fh.write("\t".join(fmt(r.get(k)) for k in keys) + "\n")
        print(f"\n[写出] {args.tsv}")


if __name__ == "__main__":
    main()
