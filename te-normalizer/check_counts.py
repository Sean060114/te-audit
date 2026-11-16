"""
check_counts.py —— 模块 2：同一份注释，为什么会有三个不同的答案

它做什么：
    读一个 TE 注释 gff3（EDTA / RepeatMasker / 其他工具都行），回答两个问题：
        1. 这条染色体上"有几个 TE"？—— 分别用三种数法算，并解释差在哪
        2. TE "占多少"？—— 对比 RepeatMasker 报的 %masked 与 bed 合并后的真实覆盖

为什么需要它（2026-11-15 实测依据）：
    同一条拟南芥 Chr1、同一份 EDTA 注释：
        gff3 全部行   3666
        元素级条目    3461   ← 差 205，全是 sub_feature / region 行
        RM 汇总表     3453   ← 再差 8，RepeatMasker 合并了相邻碎片
    覆盖度看起来也有两个数：RM 报 9.24%，bed 逐行累加 9.38%。
    但把 bed 真正合并后是 9.24% —— 9.38% 是假数，因为 bed 没合并，嵌套元件被算了两次。
    所以结论是：计数有 3 个答案（真分歧），覆盖只有 1 个答案（假分歧）。
    用户看到任何一个数都无从判断它对不对 —— 这个脚本就是把所有数摆在一起，
    并指出每个差异是真是假、来源在哪。

怎么跑：
    python check_counts.py --gff3 chr1.fa.mod.EDTA.TEanno.gff3 \\
                           --sum  chr1.fa.mod.EDTA.TEanno.sum \\
                           --bed  chr1.fa.mod.EDTA.TEanno.bed
    python check_counts.py --gff3 x.gff3 --sum x.sum --out counts.csv
"""

import argparse
import csv
import os
import sys

# 复用模块 0 的角色判断：哪些名字是"元件"，哪些是"零件"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from normalize import describe_so  # noqa: E402


# ---------------------------------------------------------------------------
# 函数 1：读 gff3，只留四样东西 —— 位置、类型、分类、角色
# ---------------------------------------------------------------------------
def parse_gff3(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            try:
                start, end = int(c[3]), int(c[4])
            except ValueError:
                continue                      # 坐标不是数字，坏行
            cls = ""
            for kv in c[8].split(";"):
                if kv.startswith("classification="):
                    cls = kv.split("=", 1)[1]
                    break
            d = describe_so(c[2])             # 第 3 列是 SO 名，用它判断角色
            rows.append({"seqid": c[0], "type": c[2], "start": start, "end": end,
                         "cls": cls, "role": d["role"], "superfamily": d["superfamily"]})
    return rows


# ---------------------------------------------------------------------------
# 函数 2：读 RepeatMasker 汇总表（TEanno.sum）
# 格式：类名 + 计数 + 掩蔽碱基数 + 百分比，两级（目 / 科），末尾 total interspersed
# ---------------------------------------------------------------------------
def parse_sum(path):
    total_len, classes, total = None, [], None
    current_order = ""
    if not path or not os.path.exists(path):
        return total_len, classes, total
    in_classes = True        # 文件有两段：Repeat Classes（要）和 Repeat Stats（逐家族，不要）
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("Repeat Stats"):
                in_classes = False
                continue
            if not in_classes:
                continue
            if s.startswith("Total Length:"):
                total_len = int(s.split()[2])
                continue
            if s.startswith("=====") or s.startswith("Class ") or not s:
                continue
            parts = s.split()
            # total interspersed 行：total interspersed 3453 2797686 9.24%
            if len(parts) >= 5 and parts[0] == "total" and parts[1] == "interspersed":
                total = {"count": int(parts[2]), "bp": int(parts[3]), "pct": parts[4]}
                continue
            if len(parts) >= 2 and parts[1] == "--":
                current_order = parts[0]           # 目级表头，如 LINE / LTR / TIR
                continue
            if parts[0] == "Total":
                continue                           # 第一段末尾的小计行，不是分类
            if len(parts) == 4 and parts[1].isdigit():
                classes.append({"order": current_order, "label": parts[0],
                                "count": int(parts[1]), "bp": int(parts[2]), "pct": parts[3]})
    return total_len, classes, total


# ---------------------------------------------------------------------------
# 函数 3：bed 覆盖度
#   注意：EDTA 的 TEanno.bed 并没有合并 —— 3461 行里有 223 处区间重叠（嵌套元件）。
#   直接累加长度会把重叠区算两遍。所以这里同时给出两个数：
#       raw   = 逐行累加（会虚高）
#       union = 按染色体合并后的真实覆盖（唯一正确的覆盖度）
#   2026-11-15 实测：raw 2,837,354 bp(9.38%)  vs  union 2,796,034 bp(9.24%)
#   union 与 RepeatMasker 汇总只差 1,652 bp —— 覆盖度其实没有分歧，分歧在计数。
# ---------------------------------------------------------------------------
def bed_coverage(path):
    if not path or not os.path.exists(path):
        return None
    raw, ivs, overlap = 0, {}, 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.split("\t")
            if len(c) < 3:
                continue
            try:
                s, e = int(c[1]), int(c[2])
            except ValueError:
                continue
            raw += e - s
            ivs.setdefault(c[0], []).append((s, e))

    union = 0
    for chrom, seg in ivs.items():
        seg.sort()
        cs, ce = seg[0]
        for s, e in seg[1:]:
            if s <= ce:                      # 相交或首尾相接 → 合并
                overlap += 1
                ce = max(ce, e)
            else:
                union += ce - cs
                cs, ce = s, e
        union += ce - cs
    return {"raw": raw, "union": union, "rows": sum(len(v) for v in ivs.values()),
            "overlap": overlap}


# ---------------------------------------------------------------------------
# 函数 4：主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="TE 注释计数口径核查（模块 2）")
    ap.add_argument("--gff3", required=True, help="TE 注释 gff3")
    ap.add_argument("--sum", default="", help="RepeatMasker 汇总表（TEanno.sum）")
    ap.add_argument("--bed", default="", help="TE 注释 bed（可选，用来核对覆盖度；脚本会自行合并重叠）")
    ap.add_argument("--genome-size", type=int, default=0,
                    help="基因组长度；不给则从 sum 里读，再不行用 gff3 最大坐标兜底")
    ap.add_argument("--out", default="", help="逐类对照表 CSV 输出路径（可选）")
    args = ap.parse_args()

    rows = parse_gff3(args.gff3)
    if not rows:
        print("gff3 里没有可解析的注释行，检查文件格式。", file=sys.stderr)
        return 1
    total_len, sum_classes, sum_total = parse_sum(args.sum)
    bed = bed_coverage(args.bed)

    genome = args.genome_size or total_len or max(r["end"] for r in rows)

    # ---- 答案一：三个数 ----
    n_all = len(rows)
    n_element = sum(1 for r in rows if r["role"] in ("element", "fragment"))
    n_skip = n_all - n_element

    print("=" * 66)
    print("问题一：这条染色体上有几个 TE？")
    print("=" * 66)
    print(f"  数法 A  gff3 全部行        {n_all}")
    print(f"  数法 B  元素级条目         {n_element}"
          f"   （扣掉 {n_skip} 行非元件）")
    if sum_total:
        d = sum_total["count"] - n_element
        rel = "少" if d < 0 else ("多" if d > 0 else "持平 ")
        print(f"  数法 C  RepeatMasker 汇总  {sum_total['count']}"
              f"   （比 B {rel} {abs(d)}）")

    # 差异来源逐项列出来，这是整个脚本的重点
    detail = {}
    for r in rows:
        if r["role"] not in ("element", "fragment"):
            detail.setdefault(r["role"], {}).setdefault(r["type"], 0)
            detail[r["role"]][r["type"]] += 1
    if detail:
        print("\n  A 与 B 的差，具体是这些行：")
        for role, types in sorted(detail.items()):
            items = "  ".join(f"{t} {c}" for t, c in sorted(types.items(), key=lambda x: -x[1]))
            print(f"    {role:<12} {sum(types.values()):>4} 行：{items}")

    # ---- 答案二：覆盖度 ----
    print("\n" + "=" * 66)
    print("问题二：TE 占这条染色体的多少？")
    print("=" * 66)
    if sum_total:
        print(f"  数法 A  RepeatMasker 汇总  {sum_total['bp']:>9,} bp   "
              f"{sum_total['bp'] / genome * 100:.2f}%")
    if bed:
        print(f"  数法 B  bed 逐行累加     {bed['raw']:>9,} bp   "
              f"{bed['raw'] / genome * 100:.2f}%   ← 虚高，别用")
        print(f"  数法 C  bed 合并去重叠   {bed['union']:>9,} bp   "
              f"{bed['union'] / genome * 100:.2f}%   ← 真实覆盖")
        if bed["overlap"]:
            print(f"\n  B 为什么虚高：这 {bed['rows']} 行 bed 里有 {bed['overlap']} 处区间相交或首尾相接"
                  f"（嵌套元件、以及相邻的碎片）。")
            print(f"  直接累加长度会把重叠段算两遍，多出 {bed['raw'] - bed['union']:,} bp。")
    if sum_total and bed:
        d = bed["union"] - sum_total["bp"]
        pct = d / genome * 100
        print(f"\n  A 与 C 相差 {d:+,} bp（{pct:+.3f} 个百分点）", end="")
        if abs(pct) < 0.05:
            print(" —— 实质一致。")
            print("  结论：覆盖度这一项其实没有分歧，真正的分歧只在「计数」。")
            print("  差的那几百到一千多 bp 来自坐标约定（1-based 含端点 vs 0-based 半开）"
                  "与边界取整，不足 0.01%。")
        else:
            print("。需要逐位点核查。")

    # ---- 逐类对照 ----
    if sum_classes:
        print("\n" + "=" * 66)
        print("逐类对照：gff3 元素级 vs RepeatMasker 汇总")
        print("=" * 66)
        gff = {}
        for r in rows:
            if r["role"] in ("element", "fragment"):
                key = r["superfamily"]
                gff[key] = gff.get(key, 0) + 1
        # RM 把 LTR/unknown 写成 LTR 目下的 unknown，换成超家族名才能对齐
        table = []
        for c in sum_classes:
            key = c["label"]
            if key == "unknown" and c["order"] == "LTR":
                key = "LTR_retrotransposon"
            table.append((c["order"], c["label"], c["count"], c["bp"],
                          gff.get(key, 0), gff.get(key, 0) - c["count"]))
        print(f"  {'目':<8}{'科':<20}{'RM 计数':>8}{'gff3':>7}{'差':>6}   {'掩蔽 bp':>10}")
        for order, label, cnt, bp, g, diff in sorted(table, key=lambda x: -x[2]):
            flag = "" if diff == 0 else ("  ←" if abs(diff) > 2 else "  ·")
            print(f"  {order:<8}{label:<20}{cnt:>8}{g:>7}{diff:>+6}   {bp:>10,}{flag}")
        unmatched = set(gff) - {("LTR_retrotransposon" if c["label"] == "unknown"
                                 and c["order"] == "LTR" else c["label"])
                                for c in sum_classes}
        if unmatched:
            print(f"\n  只在 gff3 出现、汇总表没有的科：{', '.join(sorted(unmatched))}")

        if args.out:
            with open(args.out, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["目", "科", "RM计数", "gff3计数", "差", "掩蔽bp"])
                w.writerows([[o, l, c, g, d, b] for o, l, c, b, g, d in table])
            print(f"\n  对照表已写入 {args.out}")

    # ---- 一句话结论 ----
    print("\n" + "-" * 66)
    ans = [str(n_all), str(n_element)] + ([str(sum_total["count"])] if sum_total else [])
    print(f"计数问题：{len(ans)} 个答案（{' / '.join(ans)}）。报数时必须写清用的是哪一种。")
    if sum_total and bed:
        print(f"覆盖问题：其实是 1 个答案（{sum_total['bp'] / genome * 100:.2f}%），"
              f"但用错方法会得出 {bed['raw'] / genome * 100:.2f}% 的假数。")
    print("-" * 66)
    return 0


if __name__ == "__main__":
    sys.exit(main())
