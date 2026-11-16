#!/usr/bin/env python3
"""
交付包自检。全部通过则输出 PASS，任一失败输出 FAIL 并以非零退出。

检查项：
  1  GFF3 结构：9 列、坐标 1..GENOME_LEN、start<=end、ID 唯一
  2  final.gff3 记录数与来源一致（Araport11 TE + 子片段 + EDTA 补充）
  3  补充元件与 Araport11 零重叠
  4  mask.bed 覆盖全部注释区间，且并集 bp 与软屏蔽文件小写碱基数一致
  5  软屏蔽基因组：序列长度、seqid、大小写构成
  6  家族汇总表 TOTAL 行与主注释一致
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# 脚本既可放在仓库根（数据在 ./delivery/），也可随包分发放在包内（数据在同级）。
# 自动判定，保证同一份文件在两种位置都能直接运行。
_sub = os.path.join(HERE, "delivery")
D = _sub if os.path.isdir(_sub) else HERE
GENOME_LEN = 30427671
EXPECT = {
    "TE_Araport11": 7135,
    "frag_Araport11": 7890,
    "supplement": 43,
    "mask_bp": 6020523,
    "ap_union_bp": 5081906,
}

fails, checks = [], 0


def ok(cond, msg):
    global checks
    checks += 1
    if not cond:
        fails.append(msg)
        print(f"  ✗ {msg}")
    else:
        print(f"  ✓ {msg}")


def parse_gff(path, seqid="Chr1"):
    recs, bad = [], []
    for ln, line in enumerate(open(path, encoding="utf-8"), 1):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            bad.append((ln, f"列数={len(f)}"))
            continue
        if f[0] != seqid:
            bad.append((ln, f"seqid={f[0]}"))
            continue
        try:
            s, e = int(f[3]), int(f[4])
        except ValueError:
            bad.append((ln, "坐标非整数"))
            continue
        if not (1 <= s <= e <= GENOME_LEN):
            bad.append((ln, f"坐标越界 {s}-{e}"))
        recs.append((f[2], s, e, f[6], f[8]))
    return recs, bad


def merge_iv(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def main():
    print("=== 1. GFF3 结构 ===")
    fin, bad = parse_gff(os.path.join(D, "TAIR10.Chr1.TE.final.gff3"))
    ok(not bad, f"final.gff3 无结构错误（错 {len(bad)} 行）")
    if bad[:3]:
        print("     前 3 个问题:", bad[:3])
    ids = []
    for _, _, _, _, a in fin:
        for kv in a.split(";"):
            if kv.startswith("ID="):
                ids.append(kv[3:])
    ok(len(ids) == len(set(ids)), f"ID 唯一性（{len(ids)} 个 ID）")

    for name in ("TE.Araport11", "TE.EDTA", "TE.supplement"):
        p = os.path.join(D, f"TAIR10.Chr1.{name}.gff3")
        r, b = parse_gff(p)
        ok(not b, f"{name}.gff3 无结构错误（{len(r)} 行）")

    print("\n=== 2. 记录数 ===")
    ap_recs, _ = parse_gff(os.path.join(D, "TAIR10.Chr1.TE.Araport11.gff3"))
    n_te = sum(1 for t, *_ in ap_recs if t == "transposable_element")
    n_fr = sum(1 for t, *_ in ap_recs if t == "transposon_fragment")
    sp_recs, _ = parse_gff(os.path.join(D, "TAIR10.Chr1.TE.supplement.gff3"))
    ok(n_te == EXPECT["TE_Araport11"], f"Araport11 TE = {n_te}")
    ok(n_fr == EXPECT["frag_Araport11"], f"Araport11 子片段 = {n_fr}")
    ok(len(sp_recs) == EXPECT["supplement"], f"EDTA 补充 = {len(sp_recs)}")
    ok(len(fin) == n_te + n_fr + len(sp_recs),
       f"final = {len(fin)} = {n_te}+{n_fr}+{len(sp_recs)}")

    print("\n=== 3. 零重叠 ===")
    ap_iv = [(s - 1, e) for t, s, e, _, _ in ap_recs
             if t == "transposable_element"]
    sp_iv = [(s - 1, e) for _, s, e, _, _ in sp_recs]
    # 用区间求交
    inter = 0
    a = merge_iv(ap_iv)
    for s, e in sp_iv:
        for as_, ae in a:
            lo, hi = max(s, as_), min(e, ae)
            if hi > lo:
                inter += hi - lo
    ok(inter == 0, f"补充与 Araport11 重叠 = {inter} bp（应为 0）")

    print("\n=== 4. 掩码 BED ===")
    bed = []
    for line in open(os.path.join(D, "TAIR10.Chr1.TE.mask.bed"),
                     encoding="utf-8"):
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        bed.append((int(f[1]), int(f[2])))
    mb = sum(e - s for s, e in bed)
    ok(mb == EXPECT["mask_bp"], f"BED 区间 {len(bed)} 条, {mb:,} bp")
    ok(all(bed[i][1] <= bed[i + 1][0] for i in range(len(bed) - 1)),
       "BED 已排序且互不重叠")
    # 覆盖性：所有元件区间应落在 BED 内
    covered = 0
    allv = [(s - 1, e) for t, s, e, _, _ in ap_recs if t != "transposon_fragment"]
    allv += sp_iv
    import bisect
    starts = [s for s, _ in bed]
    for s, e in allv:
        i = bisect.bisect_right(starts, s) - 1
        if i >= 0 and bed[i][1] >= e:
            covered += 1
    ok(covered == len(allv), f"注释区间被 BED 完整覆盖 {covered}/{len(allv)}")

    print("\n=== 5. 软屏蔽基因组 ===")
    fa = os.path.join(D, "TAIR10.Chr1.genome.softmasked.fa")
    hdr, n, low, chars = None, 0, 0, {}
    for line in open(fa, encoding="utf-8"):
        if line.startswith(">"):
            hdr = line[1:].strip()
            continue
        s = line.strip()
        n += len(s)
        for c in s:
            chars[c] = chars.get(c, 0) + 1
        low += sum(1 for c in s if c.islower())
    ok(n == GENOME_LEN, f"序列长度 {n:,}")
    ok(hdr.split()[0] == "Chr1", f"seqid = {hdr.split()[0]}")
    ok(low == EXPECT["mask_bp"], f"小写碱基 {low:,} = BED 并集")
    non_acgt = sum(v for k, v in chars.items() if k.upper() not in "ACGT")
    ok(non_acgt < n * 0.01, f"非 ACGT 碱基 {non_acgt:,}（{non_acgt/n*100:.2f}%）")

    print("\n=== 6. 家族汇总 ===")
    tot = None
    nfam = 0
    for line in open(os.path.join(D, "TAIR10.Chr1.TE.family_summary.tsv"),
                     encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "classification":
            continue
        if f[0] == "TOTAL":
            tot = f
        else:
            nfam += 1
    ok(tot is not None, f"存在 TOTAL 行（家族 {nfam} 个）")
    if tot:
        ok(int(tot[5]) == EXPECT["ap_union_bp"],
           f"TOTAL bp_union = {int(tot[5]):,}")
        ok(int(tot[3]) == EXPECT["TE_Araport11"],
           f"TOTAL copies = {int(tot[3])}")

    print("\n" + "=" * 46)
    if fails:
        print(f"FAIL：{len(fails)}/{checks} 项未通过")
        for m in fails:
            print("   -", m)
        return 1
    print(f"PASS：{checks} 项检查全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
