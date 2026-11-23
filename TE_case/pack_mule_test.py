#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pack-MULE 判别：元件里那个完整基因，是"被捕获的基因副本"还是"元件边界画错了"？

为什么这条能分开
----------------
水稻 Mutator 家族以捕获宿主基因著称（Pack-MULE）。一个 4 kb 的 MULE 内部装一个
完整基因，是该家族的**正常形态**；但如果那个基因在基因组里只有一份，就不可能是
"捕获副本"——只能是元件边界把唯一的原基因裹了进去，也就是结构假阳性。

判别式：
  该基因在元件之外**还有 ≥1 个高相似拷贝** → 捕获副本存在 → Pack-MULE 解释成立
  该基因在元件之外**没有任何拷贝**     → 捕获无从谈起 → 边界画错

它不依赖任何外部库，用的是本机基因组序列 + MSU 基因模型，两者都已在手。

零模型
------
水稻基因组本身有大量片段重复，随便挑一个基因也常有 paralog。所以必须给基线：
随机抽 20 个同规模基因，用同一套流程数它们的旁系拷贝数。**只有当目标基因的
拷贝数明显低于基线时，"唯一基因"才说明问题。**

做法
----
对每个外显子取 40 bp 种子在全基因组 `str.find` 枚举候选位点（40-mer 在 71 Mb 里
期望命中 0 次，实际命中都可疑），再在 ±6 位移内做无缺口比对，相似度 ≥0.85 记为
同源；位点 5 kb 内聚类。

用法
----
  python TE_case/pack_mule_test.py \
    --models osa1_r7.all_models.gff3.gz \
    --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa
"""
import argparse
import random
import re
import statistics as st
import sys
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

MSU_GFF = None          # 由 --models 传入
# 种子长度定 20 而不是 40。40-mer 在 71 Mb 里期望命中 0 次，看着很干净，但对
# 90% 相似度的旁系拷贝几乎必然漏检（需要连续 40 bp 全同，概率 0.9^40 = 1.5%）；
# 20-mer 的随机命中率仍是 71e6/4^20 ≈ 7e-5，可以忽略，而 90% 相似度下每 20 bp
# 就有 12% 的机会全同 —— 750 bp 的基因期望命中 90 个种子左右。
SEED = 20
MAX_SHIFT = 6
MIN_IDENT = 0.85
CLUSTER = 5000

# 上一轮点名的那几个元件：(EDTA seqid, start, end, 基因, 类型)
TARGETS = [
    ("OsChr1", 4586706, 4590758, "LOC_Os01g09130", "宿主"),
    ("OsChr1", 29178746, 29182856, "LOC_Os01g50800", "宿主"),
    ("OsChr1", 37694047, 37696711, "LOC_Os01g64949", "宿主"),
    ("OsChr1", 39320381, 39324964, "LOC_Os01g67660", "宿主"),
    ("OsChr1", 39321160, 39324198, "LOC_Os01g67660", "宿主"),
    ("OsChr1", 36609926, 36614941, "LOC_Os01g63170", "转座酶"),
    ("OsChr12", 1202488, 1206553, "LOC_Os12g03160", "转座酶"),
    ("OsChr12", 2628957, 2633977, "LOC_Os12g05700", "转座酶"),
]


def read_fasta(paths):
    seq, name, buf = {}, None, []
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        seq[name] = "".join(buf).upper()
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if name:
            seq[name] = "".join(buf).upper()
            name, buf = None, []
    return seq


def gene_models(path):
    """返回 gene_id -> {chrom, start, end, exons:[(s,e)], note}。"""
    out = {}
    mrna2gene = {}
    op = __import__("gzip").open(path, "rt", encoding="utf-8", errors="replace") \
        if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")
    with op as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in ("Chr1", "Chr12"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            s, e = int(f[3]) - 1, int(f[4])
            if f[2] == "gene":
                out[a.get("ID", "")] = {"chrom": f[0], "start": s, "end": e,
                                        "exons": [], "note": a.get("Note", "")}
            elif f[2] == "mRNA":
                mrna2gene[a.get("ID", "")] = a.get("Parent", "")
            elif f[2] == "exon":
                g = mrna2gene.get(a.get("Parent", ""))
                if g in out:
                    out[g]["exons"].append((s, e))
    for v in out.values():
        v["exons"].sort()
    return out


def ident(x, y):
    best = 0.0
    for sh in range(-MAX_SHIFT, MAX_SHIFT + 1):
        if sh >= 0:
            a, b = x[sh:], y
        else:
            a, b = x, y[-sh:]
        n = min(len(a), len(b))
        if n < 30:
            continue
        s = sum(1 for i in range(n) if a[i] == b[i]) / n
        if s > best:
            best = s
    return best


def loci_of(exons, genome, self_locus):
    """返回该基因在元件之外的旁系位点：[(chrom, pos, 命中外显子数, 最高相似度)]。

    每个外显子取首/中/尾三个 40 bp 种子做全基因组精确枚举——40-mer 在 71 Mb 里
    期望命中 0 次，凡是命中的都值得算一遍相似度。命中按 5 kb 分块聚合，同一块内
    命中多个外显子才算一个位点。
    """
    hits = defaultdict(dict)          # (chrom, blk) -> {exon_idx: best_ident}
    for ei, (s, e) in enumerate(exons):
        sq = genome[self_locus[0]][s:e]
        if len(sq) < SEED + 10:
            continue
        for k in range(0, len(sq) - SEED + 1, SEED):
            seed = sq[k:k + SEED]
            if len(seed) < SEED:
                continue
            for c, g in genome.items():
                pos = g.find(seed)
                n_hit = 0
                while pos >= 0 and n_hit < 200:      # 低复杂度种子直接放弃
                    n_hit += 1
                    base = pos - k
                    if base >= 0:
                        left = max(0, base - MAX_SHIFT)
                        seg = g[left: base + len(sq) + MAX_SHIFT]
                        best = 0.0
                        for d in range(0, min(2 * MAX_SHIFT + 1, len(seg) - len(sq) + 1)):
                            v = ident(seg[d:d + len(sq)], sq)
                            if v > best:
                                best = v
                        if best >= MIN_IDENT:
                            key = (c, base // CLUSTER)
                            d = hits[key]
                            if best > d.get(ei, 0.0):
                                d[ei] = best
                    pos = g.find(seed, pos + 1)
    out = []
    n_ex = len(exons)
    sc, ss, se = self_locus
    for (c, blk), d in hits.items():
        p0 = blk * CLUSTER
        # 自身必须按**整个基因跨度**排除：基因往往横跨多个 5 kb 块，只排起点所在的
        # 那一块，其余外显子会被当成"元件外的旁系位点"，相似度当然是 1.000。
        # 初版就栽在这里，基线一度显示 100% 的基因都有"近同源副本"。
        if c == sc and not (p0 + CLUSTER < ss - CLUSTER or p0 > se + CLUSTER):
            continue
        mean_id = sum(d.values()) / len(d)
        out.append({"chrom": c, "pos": p0, "n_hit": len(d), "n_ex": n_ex,
                    "mean_ident": round(mean_id, 3),
                    "frac": round(len(d) / n_ex, 2)})
    # 排序：先按"命中外显子比例"，再按平均相似度
    return sorted(out, key=lambda x: (-x["frac"], -x["mean_ident"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="C:/Users/Sean/AppData/Local/Temp/msu/all_models.gff3.gz")
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--controls", type=int, default=10)
    args = ap.parse_args()

    print("读基因组 …")
    genome = read_fasta(args.fasta)
    # 基因模型用 Chr1，基因组 FASTA 用 OsChr1。**必须统一成一套名字，不能加别名**：
    # 同一序列挂在两个键下会被遍历两次，而自身排除只比对染色体名，于是自身的命中
    # 换了个键就变成"旁系位点"，相似度 1.000——初版就这样让基线虚高到 100%。
    genome = {(k[2:] if k.startswith("Os") else k): v for k, v in genome.items()}
    for k, v in genome.items():
        print(f"  {k}: {len(v):,} bp")
    print("读 MSU 基因模型 …")
    gm = gene_models(args.models)
    print(f"  基因 {len(gm)} 个")

    print(f"\n{'元件':<20}{'基因':<18}{'类型':<7}{'外显子':>4}{'旁系位点':>7}  "
          f"{'最佳位点 外显子比例/平均相似度':<20}  近同源副本")
    rows = []
    for c, s, e, gid, kind in TARGETS:
        g = gm.get(gid)
        if not g:
            print(f"{c}:{s:<12} {gid} 未在基因模型中找到")
            continue
        ex = g["exons"]
        loci = loci_of(ex, genome, (g["chrom"], g["start"], g["end"]))
        best = loci[0] if loci else None
        near = [x for x in loci if x["frac"] >= 0.5 and x["mean_ident"] >= 0.95]
        rows.append((c, s, gid, kind, len(ex), len(loci), near))
        b = loci[0] if loci else None
        print(f"{c}:{s:<12} {gid:<18}{kind:<7}{len(ex):>4}{len(loci):>7}  "
              f"{(str(b['n_hit']) + '/' + str(b['n_ex']) + ' @ ' + format(b['mean_ident'], '.3f')) if b else '-':<20}  "
              f"{'是' if near else '否'}")

    # ---- 基线：随机同规模基因
    random.seed(7)
    pool = [k for k, v in gm.items() if 5 <= len(v["exons"]) <= 9]
    pick = random.sample(pool, min(args.controls, len(pool)))
    print(f"\n基线：随机 {len(pick)} 个 5–9 外显子的基因，同样流程数旁系位点")
    base, base_near = [], []
    for gid in pick:
        g = gm[gid]
        loci = loci_of(g["exons"], genome, (g["chrom"], g["start"], g["end"]))
        base.append(len(loci))
        base_near.append(any(x["frac"] >= 0.5 and x["mean_ident"] >= 0.95 for x in loci))
    base.sort()
    print(f"  旁系位点数 中位 {st.median(base):.0f}，均值 {sum(base) / len(base):.1f}，"
          f"为 0 的占 {sum(1 for x in base if x == 0) / len(base) * 100:.0f}%")
    print(f"  含「近同源副本」（过半外显子、平均相似度 ≥0.95）的占 "
          f"{sum(base_near) / len(base_near) * 100:.0f}%")


if __name__ == "__main__":
    main()
