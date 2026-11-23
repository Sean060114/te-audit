#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TIR 元件的"家族性"检验：这 62 条锚，在 EDTA 自己检出的 1,824 条结构 TIR 里
有没有同族同伴？

为什么需要这个而不是直接查库
--------------------------
原本要拿 Repbase/DFam 的水稻 Mutator、CACTA consensus 比对。两条路都不通：
  Repbase  商业授权，不能下载、不能分发；
  Dfam     3.9 版**没有一个 Oryza 家族**（可用 API 复现，见报告 §1）。
所以退一步问一个不用外部库也能回答的问题：

  **如果这 28 条是真的 TIR 转座子家族成员，同一家族的其它拷贝应该也在同一份
  基因组里，并且大概率也被 EDTA 检出。**

这个检验的效力与局限
--------------------
效力：真家族的臂序列（TIR）在拷贝间保守，80% 相似度的臂在 1,824 条里应当找得到
伴。孤立结构（两个偶然的反向重复被接成一个元件）不会有伴。
局限：它只能证明"EDTA 内部有没有同伴"，不能证明"这个家族叫什么名字"。同伴为 0
不等于假阳性——一个单拷贝的古老元件同样可能没有同伴。

做法
----
1. 两两比较臂序列：先按共享 8-mer 数预筛，再在 −8..+8 的位移窗口内做无缺口比对
   （臂是末端重复，5' 端锚定，只考虑小 indel 造成的位移）。
2. 相似度 ≥0.70 且长度比 ≥0.75 记为同族，并查集聚类。
3. 对每条锚报出它所在簇的大小（含自身）。

用法
----
  python TE_case/tir_family_cluster.py \
    --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
    --ledger TE_case/tir_msu_ledger.tsv \
    --tsv TE_case/tir_family.tsv
"""
import argparse
import sys
from collections import Counter, defaultdict

MIN_IDENT = 0.70      # 同族阈值：与"强档"下界一致，不另发明一套
MIN_LENRATIO = 0.75   # 臂长差超过 1/4 就不是同一个末端重复
SEED = 8              # 预筛 k-mer
MIN_SHARED = 4        # 无关联的两条 150 bp 臂期望共享 0.34 个 8-mer，4 个已经很硬
MAX_SHIFT = 8


def parse_tir(gff):
    """返回 [(chrom, start, end, cls, lt, rt)]，一行一个元件。"""
    out = []
    with open(gff, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[2].endswith("_TIR_transposon"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tir = a.get("TIR", "")
            if "_" not in tir:
                continue
            lt, rt = tir.split("_", 1)
            out.append((f[0], int(f[3]), int(f[4]), a.get("classification", ""), lt, rt))
    return out


def kmers(s, k=SEED):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def ident_ungapped(a, b):
    """无缺口比对的最高相似度（滑动窗口 −8..+8），只比重叠部分。

    为什么不做带缺口比对：这里只用来聚类，臂长 100–300 bp、错配率 ~20%，
    小 indel 造成的位移已由窗口覆盖；引入动态规划换不来更可靠的分簇，只会让
    运行时间从分钟级变成小时级。
    """
    best = 0.0
    for sh in range(-MAX_SHIFT, MAX_SHIFT + 1):
        if sh >= 0:
            x, y = a[sh:], b
        else:
            x, y = a, b[-sh:]
        n = min(len(x), len(y))
        if n < 20:
            continue
        same = sum(1 for i in range(n) if x[i] == y[i])
        s = same / n
        if s > best:
            best = s
            if best == 1.0:
                break
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--ledger", default="", help="含 verdict_code 的台账，用来分组统计")
    ap.add_argument("--tsv", default="")
    args = ap.parse_args()

    els = parse_tir(args.gff)
    if not els:
        sys.exit("没有读到 TIR 元件")
    print(f"结构 TIR 元件 {len(els)} 条；臂长中位 "
          f"{sorted(min(len(e[4]), len(e[5])) for e in els)[len(els) // 2]} bp")

    # ---- 预筛：8-mer 倒排
    idx = defaultdict(list)
    ks = []
    for i, e in enumerate(els):
        s = kmers(e[4])
        ks.append(s)
        for k in s:
            idx[k].append(i)

    # ---- 聚类
    parent = list(range(len(els)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    cand = 0
    hit = 0
    for i, e in enumerate(els):
        shared = Counter()
        for k in ks[i]:
            for j in idx[k]:
                if j > i:
                    shared[j] += 1
        for j, n in shared.items():
            if n < MIN_SHARED:
                continue
            cand += 1
            la, ra = len(e[4]), len(els[j][4])
            if min(la, ra) / max(la, ra) < MIN_LENRATIO:
                continue
            if ident_ungapped(e[4], els[j][4]) >= MIN_IDENT:
                hit += 1
                union(i, j)
    print(f"预筛候选对 {cand}；判为同族的边 {hit}")

    cluster = defaultdict(list)
    for i in range(len(els)):
        cluster[find(i)].append(i)
    sizes = {i: len(cluster[find(i)]) for i in range(len(els))}
    print(f"簇数 {len(cluster)}；最大簇 {max(len(v) for v in cluster.values())} 条")

    # ---- 按台账分组
    groups = {}
    if args.ledger:
        import csv
        with open(args.ledger, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                groups[(r["chrom"], int(r["start"]), int(r["end"]))] = r["verdict_code"]

    rows = []
    for i, e in enumerate(els):
        key = (e[0], e[1], e[2])
        rows.append({"chrom": e[0], "start": e[1], "end": e[2], "cls": e[3],
                     "arm_l": len(e[4]), "arm_r": len(e[5]),
                     "family_size": sizes[i],
                     "in_ledger": key in groups,
                     "verdict_code": groups.get(key, "")})

    if args.ledger:
        print("\n按 MSU 判定分组，看「同族同伴数」（家族大小）的分布：")
        by = defaultdict(list)
        for r in rows:
            if r["in_ledger"]:
                by[r["verdict_code"]].append(r["family_size"])
        print(f"  {'判定':<26}{'n':>3}  家族大小 众数 / 中位 / 最大 / 孤立(=1) 占比")
        for v, s in sorted(by.items(), key=lambda x: -len(x[1])):
            c = Counter(s)
            med = sorted(s)[len(s) // 2]
            iso = sum(1 for x in s if x == 1) / len(s) * 100
            print(f"  {v:<26}{len(s):>3}  {c.most_common(1)[0][0]:^6} {med:^6} {max(s):^6}   {iso:4.0f}%")

        print("\n臂长 ≥100 的子集（与长臂质疑直接相关）：")
        sub = [r for r in rows if r["in_ledger"] and min(r["arm_l"], r["arm_r"]) >= 100]
        by2 = defaultdict(list)
        for r in sub:
            by2[r["verdict_code"]].append(r["family_size"])
        for v, s in sorted(by2.items(), key=lambda x: -len(x[1])):
            iso = sum(1 for x in s if x == 1) / len(s) * 100
            print(f"  {v:<26}{len(s):>3}  中位 {sorted(s)[len(s)//2]:^4} 最大 {max(s):^4} "
                  f"孤立 {iso:4.0f}%")

    if args.tsv:
        cols = ["chrom", "start", "end", "cls", "arm_l", "arm_r", "family_size",
                "verdict_code"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                if args.ledger and not r["in_ledger"]:
                    continue
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n结果已写: {args.tsv}")


if __name__ == "__main__":
    main()
