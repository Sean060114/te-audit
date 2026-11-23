#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""这 28 条元件内部，是不是宿主基因的序列？

为什么做这个
------------
原计划拿 Repbase/DFam 的水稻 Mutator、CACTA consensus 比对，两条路都断了（见报告
§1）。剩下的判据里，只有一条既**有独立数据源**又**有零模型**：

  真 TIR 转座子占据的是基因贫乏的重复区；如果一条 4 kb 的"元件"整段盖住一个完整
  的宿主基因（含全部外显子），那它更可能是把两段反向重复接起来时顺手裹进去的
  宿主序列。

零模型是关键：水稻基因密度约每 10 kb 一个，随便取一个 4 kb 窗口也有相当概率压到
基因。所以不能只看"有没有基因"，要看**比同长随机窗口多还是少**。这也是它比
"MSU 覆盖率"那类指标强的地方——后者用同一个量定义分组又拿来比较，是循环的。

判据
----
  contained   基因模型的两端都落在元件内（完整被裹进去）
  exon_frac   元件长度中被外显子覆盖的比例
随机窗口：长度按元件长度分布重采样，每染色体 4,000 次。

用法
----
  python TE_case/gene_occupancy.py \
    --gff /tmp/msu/all_models.gff3.gz \
    --ledger TE_case/tir_msu_ledger.tsv
"""
import argparse
import gzip
import random
import statistics as st
import sys
from bisect import bisect_left
from collections import defaultdict

import os

# --- 控制台编码：Windows 默认 cp936，⚠ / → 编不出来会抛 UnicodeEncodeError，
#     报告打到一半就断；而 Git Bash 的终端是 UTF-8，按 cp936 写又是乱码。
#     按终端分两种处理，都不允许因编码中断（不可编码字符降级为 ?）。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8" if os.environ.get("MSYSTEM") else None,
                       errors="replace")
    except (AttributeError, ValueError):
        pass


def parse_genes(path, chroms):
    """返回 chrom -> (gene_starts, genes, exon_starts, exons)，均 0-based 半开。"""
    genes = defaultdict(list)
    exons = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    truncated = False
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
      try:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in chroms:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if e <= s:
                continue
            if f[2] == "gene":
                genes[f[0]].append((s, e))
            elif f[2] == "exon":
                exons[f[0]].append((s, e))
      except EOFError:
        # 下载被截断时仍然可用：gff 按染色体顺序写，先读到的染色体是完整的。
        # 这里只警告，不静默——调用方必须自己确认缺了哪条。
        truncated = True
    if truncated:
        print(f"  警告：{path} 不完整（EOFError），只统计已读到的染色体", file=sys.stderr)
    out = {}
    for c in chroms:
        g = sorted(set(genes.get(c, [])))
        x = sorted(set(exons.get(c, [])))
        out[c] = ([a for a, _ in g], g, [a for a, _ in x], x)
    return out


def union_len(iv):
    if not iv:
        return 0
    iv = sorted(iv)
    tot = 0
    cs, ce = iv[0]
    for a, b in iv[1:]:
        if a <= ce:
            ce = max(ce, b)
        else:
            tot += ce - cs
            cs, ce = a, b
    return tot + ce - cs


def measure(g, c, s, e):
    """返回 (完整被包含的基因数, 溢出边界的基因数, 外显子覆盖比例)。"""
    gs, gv, xs, xv = g[c]
    i = bisect_left(gs, s - 200000)          # 最长基因模型不会超过 200 kb
    inside = over = 0
    while i < len(gv) and gv[i][0] < e:
        a, b = gv[i]
        i += 1
        if b <= s:
            continue
        if a >= s and b <= e:
            inside += 1
        else:
            over += 1
    j = bisect_left(xs, s - 200000)
    iv = []
    while j < len(xv) and xv[j][0] < e:
        a, b = xv[j]
        j += 1
        aa, bb = max(a, s), min(b, e)
        if bb > aa:
            iv.append((aa, bb))
    return inside, over, union_len(iv) / (e - s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--draws", type=int, default=4000)
    ap.add_argument("--dedup", action="store_true",
                    help="合并起止差 ≤2 bp 的近重复记录（1,824 条里有 13 对，"
                         "62 条锚里有 3 对）——不合并会把同一位点算两次")
    ap.add_argument("--only-chrom", nargs="+", default=[],
                    help="只做这几条 EDTA 染色体，如 --only-chrom OsChr1")
    args = ap.parse_args()

    import csv
    led = list(csv.DictReader(open(args.ledger, encoding="utf-8"), delimiter="\t"))
    if args.only_chrom:
        want = set(args.only_chrom)
        led = [r for r in led if r["chrom"] in want]
        if not led:
            sys.exit(f"--only-chrom {args.only_chrom} 过滤后为空")
    if args.dedup:
        keep = []
        for r in sorted(led, key=lambda x: (x["chrom"], int(x["start"]), int(x["end"]))):
            s0, e0 = int(r["start"]), int(r["end"])
            if any(k["chrom"] == r["chrom"] and abs(int(k["start"]) - s0) <= 2
                   and abs(int(k["end"]) - e0) <= 2 for k in keep):
                continue
            keep.append(r)
        print(f"去重：{len(led)} → {len(keep)} 条（合并起止差 ≤2 bp 的近重复记录）")
        led = keep

    # EDTA 的 OsChr1 ↔ MSU 的 Chr1；只做这一次归一，别在两处各写一遍
    for r in led:
        r["msu"] = r["chrom"][2:] if r["chrom"].startswith("Os") else r["chrom"]
    chroms = sorted({r["msu"] for r in led})
    print(f"读 {args.gff}（{', '.join(chroms)}）…")
    g = parse_genes(args.gff, set(chroms))
    for c in chroms:
        print(f"  {c}: 基因 {len(g[c][1])}，外显子 {len(g[c][3])}")

    # ---- 元件
    print(f"\n{'分组':<10}{'n':>3}  完整基因 中位/最大   外显子占比 中位   压实基因的条数")
    by = defaultdict(list)
    for r in led:
        c, s, e = r["msu"], int(r["start"]) - 1, int(r["end"])
        ins, over, ex = measure(g, c, s, e)
        by[r["verdict_code"]].append((ins, over, ex, e - s))
    for v, rows in sorted(by.items(), key=lambda x: -len(x[1])):
        ins = sorted(x[0] for x in rows)
        ex = sorted(x[2] for x in rows)
        print(f"{v:<10}{len(rows):>3}  {st.median(ins):^10.1f}/{max(ins):<6} "
              f"{st.median(ex):^14.3f}   {sum(1 for x in rows if x[1])}")

    # ---- 零模型
    random.seed(0)
    null_ins, null_ex = [], []
    lens = [int(r["length"]) for r in led]
    for c in chroms:
        span = max(b for _, b in g[c][1]) + 1000
        for _ in range(args.draws):
            n = random.choice(lens)
            s = random.randint(0, max(0, span - n - 1))
            ins, over, ex = measure(g, c, s, s + n)
            null_ins.append(ins)
            null_ex.append(ex)
    null_ins.sort()
    null_ex.sort()
    print(f"\n零模型（{len(null_ex)} 个同长随机窗口）："
          f"完整基因中位 {st.median(null_ins):.0f}，"
          f"含 ≥1 个完整基因的窗口占 "
          f"{sum(1 for x in null_ins if x) / len(null_ins) * 100:.0f}%")
    print(f"  外显子占比 中位 {st.median(null_ex):.3f}，"
          f"P90 {null_ex[int(len(null_ex) * 0.9)]:.3f}")

    # ---- 逐组给出"比随机多/少"的方向
    print(f"\n{'分组':<10}{'n':>3}  完整基因 ≥1 的占比   vs 随机   "
          f"外显子占比 P(随机更高)")
    base = sum(1 for x in null_ins if x) / len(null_ins)
    for v, rows in sorted(by.items(), key=lambda x: -len(x[1])):
        withg = sum(1 for x in rows if x[0]) / len(rows)
        pex = [sum(1 for y in null_ex if y < x[2]) / len(null_ex) for x in rows]
        print(f"{v:<10}{len(rows):>3}  {withg * 100:>10.0f}%          "
              f"{(withg - base) * 100:+.0f} pp     {st.median(pex):.2f}")


if __name__ == "__main__":
    main()
