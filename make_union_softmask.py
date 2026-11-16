#!/usr/bin/env python3
"""
生成用于基因预测的软屏蔽基因组。

为什么不用 EDTA 单独的软屏蔽：
  EDTA 只覆盖 9.20%，而 Araport11 参照为 16.70%（见 EDTA_vs_Araport11_评估报告.md）。
  直接用 EDTA 的软屏蔽跑 MAKER/AUGUSTUS，会有 3.2 Mb 的真实 TE 以大写形式暴露，
  被当成基因外显子候选，产生大量假阳性基因模型。

本脚本取 EDTA ∪ Araport11 的并集做软屏蔽（重复区转小写），
  并在表头标注掩码来源与比例，便于追溯。

输入： chr1.fa.mod
       araport11_TE.gff3（Chr1）
       chr1.fa.mod.EDTA.TEanno.gff3
输出： chr1.fa.mod.union.softmasked.fa
"""
import sys
import numpy as np

FASTA_IN = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod"
FASTA_OUT = "chr1.fa.mod.union.softmasked.fa"
AP_GFF = "araport11_TE.gff3"
ED_GFF = ("//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/"
          "chr1.fa.mod.EDTA.TEanno.gff3")
SKIP = {"target_site_duplication", "long_terminal_repeat", "primer_binding_site"}
WRAP = 60


def read_fasta(path):
    name, buf = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip()
            else:
                buf.append(line.strip())
    return name, "".join(buf)


def main():
    name, seq = read_fasta(FASTA_IN)
    n = len(seq)
    print(f"输入序列: {name[:60]}")
    print(f"长度: {n:,} bp")

    ap = np.zeros(n, dtype=bool)
    for line in open(AP_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != "Chr1" or f[2] != "transposable_element":
            continue
        ap[int(f[3]) - 1:int(f[4])] = True

    ed = np.zeros(n, dtype=bool)
    for line in open(ED_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] in SKIP:
            continue
        ed[int(f[3]) - 1:int(f[4])] = True

    union = ap | ed
    print(f"\nAraport11 掩码 : {ap.sum():>10,} bp ({ap.sum()/n*100:5.2f}%)")
    print(f"EDTA 掩码      : {ed.sum():>10,} bp ({ed.sum()/n*100:5.2f}%)")
    print(f"并集           : {union.sum():>10,} bp ({union.sum()/n*100:5.2f}%)")
    print(f"EDTA 新增      : {(union & ~ap).sum():>10,} bp "
          f"({(union & ~ap).sum()/n*100:5.2f}%)")
    print(f"Araport11 补充 : {(union & ~ed).sum():>10,} bp "
          f"({(union & ~ed).sum()/n*100:5.2f}%)")

    lower_lut = str.maketrans("ACGTNacgtn", "acgtnacgtn")
    chars = list(seq)
    idx = np.flatnonzero(union)
    for i in idx:
        chars[i] = chars[i].translate(lower_lut)
    out_seq = "".join(chars)

    with open(FASTA_OUT, "w", encoding="utf-8", newline="\n") as fh:
        # chr1.fa.mod 的头被 EDTA 预处理截成了 ">1"，这里补回原始 TAIR10 描述
        desc = ("dna:chromosome chromosome:TAIR10:1:1:"
                f"{n}:1 REF")
        fh.write(f">1 {desc} softmasked=Araport11_TE_union_EDTA_"
                 f"{union.sum()/n*100:.2f}pct\n")
        for i in range(0, n, WRAP):
            fh.write(out_seq[i:i + WRAP] + "\n")

    # 自检
    low = sum(1 for c in out_seq if c.islower())
    print(f"\n输出: {FASTA_OUT}")
    print(f"小写碱基数: {low:,}  与并集一致: {'是' if low == union.sum() else '否'}")


if __name__ == "__main__":
    sys.exit(main())
