#!/usr/bin/env python3
"""用基因组序列独立复算 EDTA 结构法 LTR 元件的 TSD 与末端基序。

为什么需要它：EDTA 的 intact.gff3 自带 `tsd=`、`motif=`、`ltr_identity=` 三个
字段，但那是 EDTA 自己算的。要判断「这批注释可信」，必须能从序列重新算一遍，
而不是引用它自己的字段。本脚本只做可独立验证的两项（TSD、末端基序）；
`ltr_identity` 需要 LTR 子坐标，而 LTR_retriever 的 pass.list 未随产物回传，
暂时无法复算——这一点在报告里必须写明。

字段语义（据 EDTA v2.2.2 bin/make_gff3.pl 与 LTR_retriever bin/make_gff3.pl
源码确认，2026-11-23）：
  1. 坐标 1-based 闭区间。
  2. `repeat_region` 的起止 = 5′TSD 第一个碱基 … 3′TSD 最后一个碱基，即
     **TSD 含在元件坐标内**；`tsd=NA` 时才退化为左右 LTR 的起止。
  3. `tsd=` 写的是**元件方向**上的序列：负链元件需反向互补后再比对。
  4. `motif=` 是 LTR 的末端基序（如 TGCA），拆成「前 2 nt + 后 2 nt」分别落在
     5′LTR 起点与 3′LTR 终点，**不是**边界处连续出现的 4 个碱基。
  5. 一个 LTR 元件在 intact.gff3 里占 6 行（repeat_region / lTSD / lLTR /
     LTRRT / rLTR / rTSD），**只有 repeat_region 行代表元件**。

用法：
    python verify_ltr_structure.py \
        --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
        --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa
"""

import argparse
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(paths):
    seqs = {}
    for p in paths:
        name, buf = None, []
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(buf).upper()
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf).upper()
    return seqs


def parse_gff(path):
    """只取 repeat_region 行——一个 LTR 元件占 6 行，其余 5 行不是元件。"""
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "repeat_region":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if not a.get("classification", "").startswith("LTR"):
                continue
            out.append({
                "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                "strand": f[6], "tsd": a.get("tsd", "NA"),
                "motif": a.get("motif", ""), "ident": a.get("ltr_identity", ""),
            })
    return out


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--tol", type=int, default=1, help="3′TSD 允许的错配数（默认 1）")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    elems = parse_gff(args.gff)
    if not elems:
        sys.exit("没有找到 LTR 的 repeat_region 行")

    n = ok5 = ok3 = ok3tol = okboth = okmot = 0
    idents = []
    for e in elems:
        chrom = e["chrom"]
        if chrom not in seqs:
            # seqid 可能带前缀（OsChr1 vs Chr1），退回后缀匹配
            cand = [k for k in seqs if k.endswith(chrom) or chrom.endswith(k)]
            if not cand:
                continue
            chrom = cand[0]
        sq = seqs[chrom]
        s, en = e["start"], e["end"]
        t = e["tsd"]
        if not t or t == "NA":
            continue
        n += 1
        L = len(t)
        left = sq[s - 1:s - 1 + L]
        right = sq[en - L:en]
        if e["strand"] == "-":            # tsd 按元件方向书写
            left, right = rc(left), rc(right)
        if left == t:
            ok5 += 1
        if right == t:
            ok3 += 1
        if hamming(right, t) <= args.tol:
            ok3tol += 1
        if left == t and hamming(right, t) <= args.tol:
            okboth += 1
        m = e["motif"]
        if len(m) == 4:
            # 前 2 nt 在 5′LTR 起点（紧接 5′TSD），后 2 nt 在 3′LTR 终点（紧邻 3′TSD）
            if (sq[s - 1 + L:s + 1 + L] == m[:2]
                    and sq[en - L - 2:en - L] == m[2:]):
                okmot += 1
        if e["ident"]:
            try:
                idents.append(float(e["ident"]))
            except ValueError:
                pass

    pct = lambda k: f"{k}/{n} ({k / n * 100:.1f}%)" if n else "n/a"
    print(f"LTR 元件（repeat_region 行）: {len(elems)}；带 TSD 可复算: {n}")
    print(f"  5′ TSD 与序列完全一致          : {pct(ok5)}")
    print(f"  3′ TSD 与序列完全一致          : {pct(ok3)}")
    print(f"  3′ TSD 允许 {args.tol} 处错配      : {pct(ok3tol)}")
    print(f"  两侧同时成立（5′ 严格 + 3′ 容错）: {pct(okboth)}")
    print(f"  末端基序（motif 前 2 + 后 2）    : {pct(okmot)}")
    if idents:
        idents.sort()
        print(f"  ltr_identity 中位/最小（字段值，未复算）: "
              f"{idents[len(idents) // 2]:.4f} / {idents[0]:.4f}")


if __name__ == "__main__":
    main()
