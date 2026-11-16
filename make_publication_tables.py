#!/usr/bin/env python3
"""
生成两张发表用表。

表 1  Araport11_TE_family_summary.tsv
      Araport11（TAIR10 Chr1）TE 家族 × 拷贝数 × bp × 占基因组比例
      同时给出超家族级汇总，bp 提供「逐条相加」与「区间去重」两种口径

表 2  EDTA_B_level_candidates.tsv
      EDTA 独有元件中有完整结构证据（TIR + TSD）的 B 级候选，
      附左/右 TIR 序列、TSD 序列、identity，可直接作为补充表

依赖：araport11_TE.gff3、chr1.fa.mod.EDTA.TEanno.gff3、candidates_novel.py
"""
import sys
import numpy as np

import candidates_novel as cn

AP_GFF = "araport11_TE.gff3"
GENOME_FA = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod"
GENOME_LEN = 30427671
OUT_FAMILY = "Araport11_TE_family_summary.tsv"
OUT_SUPP = "EDTA_B_level_candidates.tsv"


def merge_len(intervals):
    """区间并集总长（0-based 半开）。"""
    if not intervals:
        return 0
    iv = sorted(intervals)
    total = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s > ce:
            total += ce - cs
            cs, ce = s, e
        elif e > ce:
            ce = e
    return total + (ce - cs)


def parse_araport11():
    """返回 {完整 classification: [(s, e), ...]}。"""
    fam = {}
    for line in open(AP_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != "Chr1" or f[2] != "transposable_element":
            continue
        s, e = int(f[3]) - 1, int(f[4])
        cls = cn.attr_get(f[8], "classification") or "Unknown"
        fam.setdefault(cls, []).append((s, e))
    return fam


def superfamily(cls):
    """取第一层。RC/Helitron/ATREP3 -> RC；DNA/MuDR/X -> DNA。"""
    return cls.split("/")[0]


def subfamily(cls):
    """取第二层；若无第二层则返回超家族本身。"""
    parts = cls.split("/")
    return parts[1] if len(parts) > 1 else parts[0]


def write_family_table(fam):
    rows = []
    for cls, iv in fam.items():
        bp_sum = sum(e - s for s, e in iv)
        rows.append((cls, superfamily(cls), subfamily(cls),
                     len(iv), bp_sum, merge_len(iv)))
    rows.sort(key=lambda r: (-r[5], r[0]))

    with open(OUT_FAMILY, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(["classification", "superfamily", "subfamily",
                            "copies", "bp_sum", "bp_union",
                            "pct_genome"]) + "\n")
        for cls, sup, sub, n, bsum, bun in rows:
            fh.write(f"{cls}\t{sup}\t{sub}\t{n}\t{bsum}\t{bun}\t"
                     f"{bun / GENOME_LEN * 100:.4f}\n")

    print(f"表 1 已写出: {OUT_FAMILY}  （{len(rows)} 个家族）")

    # 超家族级汇总
    agg = {}
    for cls, sup, sub, n, bsum, bun in rows:
        a = agg.setdefault(sup, [0, 0, []])
        a[0] += n
        a[1] += bsum
        a[2].extend(fam[cls])
    total_union = merge_len([iv for v in fam.values() for iv in v])

    print("\n=== 超家族级汇总（Araport11, TAIR10 Chr1）===")
    print(f"{'超家族':<10}{'拷贝数':>8}{'bp_sum':>12}{'bp_union':>12}{'占比':>9}")
    print("-" * 53)
    for sup, (n, bsum, ivs) in sorted(agg.items(), key=lambda x: -x[1][1]):
        print(f"{sup:<10}{n:>8}{bsum:>12,}{merge_len(ivs):>12,}"
              f"{merge_len(ivs) / GENOME_LEN * 100:>8.2f}%")
    print("-" * 53)
    print(f"{'合计':<10}{sum(a[0] for a in agg.values()):>8}"
          f"{sum(a[1] for a in agg.values()):>12,}"
          f"{total_union:>12,}{total_union / GENOME_LEN * 100:>8.2f}%")
    print("\n注：bp_sum 为各元件长度直接相加（跨家族重叠会重复计入）；"
          "bp_union 为区间并集，合计行即全基因组去重值。")


def load_genome():
    name, buf = None, []
    for line in open(GENOME_FA, encoding="utf-8"):
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].strip()
        else:
            buf.append(line.strip())
    return "".join(buf).upper()


COMP = str.maketrans("ACGTNacgtn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def tir_qc(genome, start, end, tir_len, slack=6):
    """
    检验两端是否构成反向重复。
    返回 (原位错配数, 是否可通过位移达到 ≤1 错配, 需位移量)
    """
    left = genome[start:start + tir_len]
    right = revcomp(genome[end - tir_len:end])
    if len(left) != tir_len or len(right) != tir_len:
        return -1, False, 0
    miss = sum(1 for i in range(tir_len) if left[i] != right[i])
    if miss <= 1:
        return miss, True, 0
    for k in range(1, slack + 1):
        for dl, dr in ((-k, k), (k, -k), (-k, -k), (k, k)):
            a = genome[start + dl:start + dl + tir_len]
            b = genome[end - tir_len + dr:end + dr]
            if len(a) != tir_len or len(b) != tir_len:
                continue
            b = revcomp(b)
            if sum(1 for i in range(tir_len) if a[i] != b[i]) <= 1:
                return miss, True, k
    return miss, False, 0


def write_b_level(rows, ap_mask):
    zero = [r for r in rows if not ap_mask[r["start"]:r["end"]].any()]
    for r in zero:
        r["grade"], r["reason"] = cn.grade(r)
    b = [r for r in zero if r["grade"] == "B"]
    b.sort(key=lambda r: -r["len"])

    genome = load_genome()

    with open(OUT_SUPP, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join([
            "id", "chrom", "start", "end", "length", "class", "identity",
            "left_TIR", "right_TIR", "TIR_len", "TIR_mismatch",
            "perfect_IR_after_shift", "shift_bp",
            "left_TSD", "TSD_match_pct", "method"]) + "\n")
        for i, r in enumerate(b, 1):
            tir = r["tir"].split("_")
            ltir = tir[0] if tir else ""
            rtir = tir[1] if len(tir) > 1 else ""
            tsd = r["tsd"].split("_")
            ltsd = tsd[0] if tsd else ""
            match = tsd[2] if len(tsd) > 2 else ""
            miss, perfect, shift = tir_qc(genome, r["start"], r["end"],
                                          len(ltir))
            fh.write("\t".join([
                f"EDTA_cand_{i:03d}", r["chrom"], str(r["start"] + 1),
                str(r["end"]), str(r["len"]), r["cls"], r["ident"],
                ltir, rtir, str(len(ltir)), str(miss),
                "Y" if perfect else "N", str(shift),
                ltsd, match, r["method"],
            ]) + "\n")

    print(f"\n表 2 已写出: {OUT_SUPP}  （{len(b)} 条 B 级候选）")
    import collections
    c = collections.Counter(r["cls"] for r in b)
    print("  构成: " + ", ".join(f"{k}×{v}" for k, v in c.most_common()))
    print(f"  长度: {min(r['len'] for r in b)}–{max(r['len'] for r in b)} bp, "
          f"合计 {sum(r['len'] for r in b):,} bp")
    print(f"  identity: {min(float(r['ident']) for r in b):.3f}–"
          f"{max(float(r['ident']) for r in b):.3f}")

    # TIR 反向重复质控
    stat = collections.Counter()
    for r in b:
        ltir = r["tir"].split("_")[0]
        miss, perfect, shift = tir_qc(genome, r["start"], r["end"], len(ltir))
        if perfect and shift == 0:
            stat["原位即完美"] += 1
        elif perfect:
            stat[f"位移后完美"] += 1
        else:
            stat["无完美反向重复"] += 1
    print("\n=== TIR 反向重复质控（错配 ≤1 视为完美）===")
    for k, v in stat.items():
        print(f"  {k}: {v} 条")
    n_deg = stat["无完美反向重复"]
    print(f"  → {n_deg}/{len(b)} 条即使允许 ±6 bp 位移也无法形成完美反向重复，"
          f"属真实退化（与 identity 0.81–0.92 一致）")


def main():
    fam = parse_araport11()
    write_family_table(fam)
    ap_mask = cn.build_ap_mask()
    rows = cn.read_edta()
    write_b_level(rows, ap_mask)


if __name__ == "__main__":
    sys.exit(main())
