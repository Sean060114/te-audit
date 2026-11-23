#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐元件台账：每个 EDTA 结构 LTR 元件，在 MSU 侧有没有对应的超家族判断。

与 compare_rice_te.py 的分工
----------------------------
那份脚本算的是**染色体级 bp 指标**（并集密度、召回、精确率），没有"元件↔元件"
的对应关系。这里要的是**逐条台账**：这一个元件，MSU 说它是什么、覆盖了多少。

复用（import 自 compare_rice_te）：
  bucket_msu    MSU 的 Note → 类别。顺序敏感的规则链，SINE 那个坑已修过。
  bucket_ed     EDTA 的 classification → 类别。
  merge         区间并集。
  covered_len   区间被覆盖的碱基数。

不复用：
  parse_ed      它的 SKIP_ED 含 repeat_region —— 那正是元件本体，一条都读不到；
                且不区分 method=structural|homology，会把同源碎片混进来。
  parse_msu     它把 Note 按染色体聚合成 bp，元件级拿不到"哪条 Note 压在它身上"。
  analyze       整体是 bp 级并集/交集，没有元件身份。

一条口径上的硬约束
------------------
覆盖率必须用**并集**：MSU 是碎片化注释（Chr1 就有 57k 行），按"单行覆盖 ≥50%"判
会有 207/420 元件被误判成"无对应"。实测并集口径下 10%/30%/50% 三档结果稳定。

用法
----
  python TE_case/link_msu.py \
    --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
    --msu ref/msu/osa1_r7.asm.repeat_masked.gff3.gz \
    --tsv TE_case/ltr_msu_ledger.tsv
"""
import argparse
import gzip
import os
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from urllib.parse import unquote

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import compare_rice_te as C          # noqa: E402  只取映射与区间工具

# 台账里写 ASCII 码、控制台里印中文：数据文件要能被下游直接读（本仓库既有 TSV
# 一律无 BOM、无中文表头），而人看的是屏幕。
VERDICT_CODE = {
    "超家族一致": "agree_superfamily",
    "超家族冲突": "conflict_superfamily",
    "MSU 只到 LTR 层": "msu_ltr_only",
    "EDTA 未定/MSU 定": "ed_unknown_msu_resolved",
    "双方均未定": "both_unresolved",
    "MSU 判为非 LTR": "msu_nonltr_label",
    "MSU 覆盖不足": "insufficient_cov",
    "空区间": "empty",
}

LTR_RESOLVED = (C.B_LTR_GYPSY, C.B_LTR_COPIA)
LTR_ALL = LTR_RESOLVED + (C.B_LTR_UNK,)


# --------------------------------------------------------------------- 读入
def read_elements(gff):
    """只取 repeat_region 行里的 LTR 元件。

    一个元件在 intact.gff3 里占 6 行（repeat_region / lTSD / lLTR / LTRRT /
    rLTR / rTSD），只有 repeat_region 代表元件本体，且带 tsd / motif /
    ltr_identity 三个字段。

    身份用 (chrom, start, end) 而不是 Name —— Name 是**家族名**，
    TE_00000794 在 Chr1+Chr12 上出现 9 次。
    """
    out = []
    with open(gff, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "repeat_region":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            cls = a.get("classification", "")
            if not cls.startswith("LTR"):
                continue
            out.append({
                # gff 的 1-based 闭区间原样保留（人工核对时看的就是这两个数）；
                # 所有区间算术在 classify 里转成 0-based 半开，与 MSU 侧一致。
                "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                "strand": f[6], "cls": cls,
                "ident": a.get("ltr_identity", "NA"),
                "tsd": a.get("tsd", "NA"),
                "motif": a.get("motif", "NA"),
                "name": a.get("Name", ""),      # 家族名，只作参考，不当主键
            })
    return out


def read_msu(path, cmap):
    """返回 EDTA seqid -> (starts, rows)，rows 为 (start0, end, bucket, note)。

    MSU 的 chr 名与 EDTA 的 seqid 不同（Chr1 vs OsChr1）。cmap 是
    MSU 名 → EDTA seqid；不在 cmap 里的染色体直接丢弃——**不要先归一再用
    原名去过滤**，那样会把所有行都过滤掉（本脚本初版就这样，MSU 行数 0）。
    """
    rows = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            c = cmap.get(p[0])
            if c is None:
                continue
            m = re.search(r"Note=([^;]+)", p[8])
            note = unquote(m.group(1)) if m else ""
            s, e = int(p[3]) - 1, int(p[4])
            if e <= s:
                continue
            rows[c].append((s, e, C.bucket_msu(note), note))
    out = {}
    for c, v in rows.items():
        v.sort(key=lambda x: x[0])
        out[c] = ([r[0] for r in v], v)
    return out


# --------------------------------------------------------------------- 判定
def classify(elem, msu, min_cov):
    """返回 (verdict, dominant, cov_by_bucket, n_hits, top_note)。

    判定顺序（**先认一致，再认冲突，最后才降级为"没解出来"**）：
      一致 > 冲突 > MSU 只到 LTR 层 > MSU 判为非 LTR > 覆盖不足
    顺序不能反：若先判 "MSU/unknown 达标" 再判冲突，一个 Gypsy 元件同时压着
    Gypsy 与 unknown 两条时会漏掉真实的冲突。
    """
    eth = C.bucket_ed(elem["cls"], "repeat_region")
    s0 = elem["start"] - 1                       # 1-based 闭 → 0-based 半开
    L = elem["end"] - s0
    if L <= 0:
        return "空区间", "", {}, 0, ""

    c = elem["chrom"]
    starts, rows = msu.get(c, ([], []))
    # 只需看 start 落在 [s - 最长行, e) 的行
    maxlen = max((e - s for s, e, _, _ in rows), default=0)
    i = bisect_left(starts, s0 - maxlen)
    cov = Counter()
    noted = defaultdict(Counter)
    n = 0
    while i < len(rows) and rows[i][0] < elem["end"]:
        s, e, b, note = rows[i]
        if e > s0:
            n += 1
            cov[b] += min(e, elem["end"]) - max(s, s0)
            noted[b][note] += 1
        i += 1

    # 上面按行累加过一遍，只能用来判断"有没有命中"；覆盖率必须用并集重算，
    # 否则同一桶里相邻/重叠的 MSU 行会被重复计入，超过 100%。
    real = _union_cov(elem, s0, msu, cov.keys())

    ok = {b: v for b, v in real.items() if v >= min_cov * L}
    if not ok:
        dom = max(real, key=real.get) if real else ""
        return "MSU 覆盖不足", dom, real, n, _top_note(noted, dom)

    ed_unres = eth == C.B_LTR_UNK or eth not in LTR_RESOLVED
    others = [b for b in ok if b in LTR_RESOLVED and b != eth]

    if eth in LTR_RESOLVED and eth in ok:
        return "超家族一致", eth, real, n, _top_note(noted, eth)
    if others and not ed_unres:
        # 两源都给出了超家族，且不同 —— 这才是真冲突
        dom = max(others, key=lambda b: ok[b])
        return "超家族冲突", dom, real, n, _top_note(noted, dom)
    if ed_unres and others:
        # EDTA 自己没定（LTR/unknown），MSU 定了 —— 不是冲突，是单边可判
        dom = max(others, key=lambda b: ok[b])
        return "EDTA 未定/MSU 定", dom, real, n, _top_note(noted, dom)
    if C.B_LTR_UNK in ok and not ed_unres:
        return "MSU 只到 LTR 层", C.B_LTR_UNK, real, n, _top_note(noted, C.B_LTR_UNK)
    if C.B_LTR_UNK in ok:
        return "双方均未定", C.B_LTR_UNK, real, n, _top_note(noted, C.B_LTR_UNK)

    dom = max(ok, key=ok.get)
    return "MSU 判为非 LTR", dom, real, n, _top_note(noted, dom)


def _union_cov(elem, s0, msu, buckets):
    """按桶算并集覆盖率（重叠的 MSU 行只算一次）。s0 是 0-based 起点。"""
    starts, rows = msu.get(elem["chrom"], ([], []))
    maxlen = max((e - s for s, e, _, _ in rows), default=0)
    i = bisect_left(starts, s0 - maxlen)
    per = defaultdict(list)
    while i < len(rows) and rows[i][0] < elem["end"]:
        s, e, b, _ = rows[i]
        if e > s0:
            per[b].append((max(s, s0), min(e, elem["end"])))
        i += 1
    out = {}
    for b, ivs in per.items():
        if b in buckets:
            out[b] = C.span(C.merge(ivs))
    return out


def _top_note(noted, bucket):
    if not bucket or bucket not in noted or not noted[bucket]:
        return ""
    return noted[bucket].most_common(1)[0][0]


# --------------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="EDTA *.intact.gff3")
    ap.add_argument("--msu", required=True, help="MSU repeat_masked.gff3(.gz)")
    ap.add_argument("--min-cov", type=float, default=0.3,
                    help="判定「有对应」所需的最小并集覆盖率（默认 0.3）")
    ap.add_argument("--tsv", default="", help="台账输出路径")
    ap.add_argument("--chrom", action="append", default=[],
                    help="EDTA seqid=MSU chr 的映射，如 --chrom OsChr1=Chr1；"
                         "缺省按去掉前缀 Os 处理")
    args = ap.parse_args()
    for f in (args.gff, args.msu):
        if not os.path.isfile(f):
            sys.exit(f"找不到输入文件：{f}\n"
                     f"仓库根目录下应能访问 rice/…intact.gff3 与 ref/msu/…gff3.gz；"
                     f"序列与参考数据不入库，换机器需先放回原位")

    elems = read_elements(args.gff)
    if not elems:
        sys.exit("没有读到 LTR 元件——检查 gff 是否为 *.intact.gff3，且含 repeat_region 行")

    cmap = {}
    for e in elems:
        ed = e["chrom"]
        pair = [x for x in args.chrom if x.split("=")[0] == ed]
        cmap[ed] = pair[0].split("=")[1] if pair else ed[2:] if ed.startswith("Os") else ed
    # cmap 是 EDTA seqid → MSU 名；read_msu 要的是反向表
    rev = {ms: ed for ed, ms in cmap.items()}
    msu = read_msu(args.msu, rev)

    ledger = []
    over = 0
    for e in elems:
        verdict, dom, cov, n, note = classify(e, msu, args.min_cov)
        L = e["end"] - e["start"] + 1
        if any(v > L for v in cov.values()):
            # 并集覆盖率不可能超过元件长度。超了就是同桶内区间被重复计入，
            # 是分区间的老毛病，宁可中止也不要写出一个 >1 的数。
            over += 1
        ledger.append({
            "chrom": e["chrom"], "start": e["start"], "end": e["end"],
            "strand": e["strand"], "length": L, "ed_class": e["cls"],
            "ltr_identity": e["ident"], "tsd": e["tsd"],
            "verdict": verdict, "verdict_code": VERDICT_CODE[verdict],
            "msu_dominant": dom,
            "cov_gypsy": round(cov.get(C.B_LTR_GYPSY, 0) / L, 3),
            "cov_copia": round(cov.get(C.B_LTR_COPIA, 0) / L, 3),
            "cov_unknown": round(cov.get(C.B_LTR_UNK, 0) / L, 3),
            "n_msu_hits": n, "msu_note": note, "family": e["name"],
        })

    if over:
        sys.exit(f"{over} 个元件的某桶覆盖率超过 100%——并集没生效，结果不可用")

    # ---- 汇总
    print(f"EDTA 结构 LTR 元件: {len(ledger)}；MSU 行: {sum(len(v[1]) for v in msu.values())}")
    print(f"判定阈值：并集覆盖率 ≥ {args.min_cov:.0%}\n")
    tab = Counter((r["ed_class"], r["verdict"]) for r in ledger)
    verdicts = ["超家族一致", "超家族冲突", "MSU 只到 LTR 层", "EDTA 未定/MSU 定",
                "双方均未定", "MSU 判为非 LTR", "MSU 覆盖不足"]
    classes = sorted({r["ed_class"] for r in ledger})
    w = max(len(v) for v in verdicts)
    print(" " * 14 + "".join(f"{v:<{w + 3}}" for v in verdicts) + "合计")
    for cl in classes:
        row = [tab.get((cl, v), 0) for v in verdicts]
        print(f"  {cl:<12}" + "".join(f"{x:<{w + 3}}" for x in row) + f"{sum(row)}")
    print(f"  {'合计':<12}" + "".join(f"{tab_total(tab, v):<{w + 3}}" for v in verdicts)
          + f"{len(ledger)}")

    # 两源都给出超家族时的严格一致率——这是能写进报告的分子分母
    both = [r for r in ledger if r["ed_class"] in LTR_RESOLVED
            and r["verdict"] in ("超家族一致", "超家族冲突")]
    agree = sum(1 for r in both if r["verdict"] == "超家族一致")
    print(f"\n两源都定超家族: {len(both)} 条，一致 {agree} = "
          f"{agree / len(both) * 100:.1f}%" if both else "\n两源都定超家族的样本为 0")

    if args.tsv:
        cols = ["chrom", "start", "end", "strand", "length", "ed_class",
                "ltr_identity", "tsd", "verdict_code", "msu_dominant",
                "cov_gypsy", "cov_copia", "cov_unknown", "n_msu_hits",
                "msu_note", "family"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in ledger:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n台账已写: {args.tsv}（{len(ledger)} 行）")


def tab_total(tab, v):
    return sum(n for (cl, vv), n in tab.items() if vv == v)


if __name__ == "__main__":
    main()
