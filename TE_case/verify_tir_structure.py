#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TIR/MITE 结构证据的独立复算与分级。

和 verify_ltr_structure.py 同一套办法：不信任 gff 里的字段，用基因组序列把三件事
逐条重算——两翼 TIR 是否真在两端、是否构成反向重复、TSD 是否真在两侧。

为什么要单独一份而不是改 LTR 那份
----------------------------------
两者差三处，且都是"照抄会算错"的那种：

1. **坐标含不含 TSD 正好相反。** LTR 的 repeat_region 起止**含**两侧 TSD；
   TIR 的元件行起止**不含** TSD——TSD 在元件外侧（实测 1,824/1,824 两侧都能在
   元件外原样找到）。照抄 LTR 的取值逻辑，TSD 会整体错位到元件内部。
2. **右翼 TIR 是基因组正向写的，不是反向互补。** 字段里
   `TIR=GAGAAGAGAAG_ATTCTCTTCGC` 的两段都是 5'→3'（正链）序列；
   判断"是不是反向重复"要拿左翼和 **rc(右翼)** 比。
3. **`identity=` 不是"两翼逐位相似度"的精确值**，也**不是 TIR-Learner 的 TIR identity**
   （后者在官方样例里只有 90.0 / 100.0 两档）。它是 EDTA 从 bed 第 6 列写进来的值，
   与逐位 hamming 平均差 1.8 pp，个别差到 0.6——见报告 5.2 的 6 条例外。

分级不拍阈值，用随机命中期望
----------------------------
TIR 只有 11–20 bp，而 2 bp 的 TSD 随机出现概率是 1/16。所以"有 TIR、有 TSD"
本身不是证据，必须问：**这么长、这么像的一对反向重复，在一个 70.8 Mb 的基因组里
随机能撞上几次？**

    期望次数 E = N_pos × Σ_{k=0..m} C(L,k)·3^k / 4^L

L 为臂长、m 为错配数、N_pos 取两染色体总长。E 越小越不可能是巧合：

分级最后并没有用 E 当唯一判据（见 tier_of 的注释）：E 只做参考列。理由是
"1,556 bp、84% 相似"这类长臂按 E 排序靠后，却显然是硬证据。

用法
----
  python TE_case/verify_tir_structure.py \
    --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
    --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa \
    --tsv TE_case/tir_structure.tsv
"""
import argparse
import math
import sys
from collections import Counter

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

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


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


def parse_tir_gff(path):
    """只取结构 TIR 元件行（type 以 _TIR_transposon 结尾）。

    一个 TIR 元件在 gff 里就是一行（不像 LTR 占 6 行），行内同时带
    TIR=、TSD=、identity= 三个字段。
    """
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[2].endswith("_TIR_transposon"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tir = a.get("TIR", "")
            tsd = a.get("TSD", "")
            if "_" not in tir:
                continue
            lt, rt = tir.split("_", 1)
            parts = tsd.split("_")
            ts, te = (parts[0], parts[1]) if len(parts) >= 2 else ("", "")
            out.append({
                "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                "strand": f[6], "type": f[2],
                "cls": a.get("classification", ""),
                "lt": lt, "rt": rt, "tsd_l": ts, "tsd_r": te,
                "tsd_score": parts[2] if len(parts) >= 3 else "",
                "ident": a.get("identity", "NA"),
            })
    return out


def log10_chance(L, m):
    """一对长度 L、错配 ≤m 的反向重复在一个位置上随机出现的概率（log10）。

    Σ_{k≤m} C(L,k)·3^k/4^L。用 lgamma 在对数域累加，避免 L 大时溢出。
    """
    if L <= 0:
        return 0.0
    terms = []
    for k in range(0, min(m, L) + 1):
        # lgamma 返回自然对数，后面几项都是 log10 —— 必须换算，否则整条式子
        # 差一个 ln(10) 的因子（初版就这样，把 P 算大了约 2.3 个数量级）
        log_c = (math.lgamma(L + 1) - math.lgamma(k + 1)
                 - math.lgamma(L - k + 1)) / math.log(10)
        terms.append(log_c + k * math.log10(3) - L * math.log10(4))
    mx = max(terms)
    return mx + math.log10(sum(10 ** (t - mx) for t in terms))


def tier_of(arm, sim):
    """按"臂长 + 相似度"两级判，不用 E 当唯一判据。

    为什么不能只用 E：E 是"某个位置与某个对位"的期望次数，按它排序会把
    1,556 bp / 84% 相似的长臂判成"弱"（错配多，E 不低），而 156 个错配散在
    1.5 kb 上、两端还能对上，这在基因组里近乎不可能。长度本身必须是判据。

    为什么也不能只看相似度：11 bp 的臂即使 100% 一致，在 71 Mb 里预期也撞得上
    十几次（71e6 / 4^11 ≈ 17）。短臂必须靠"错配少 + TSD 长"补。

    规则（可对评审逐条解释）：
      强  臂长 ≥100 且相似度 ≥0.70（长反向重复本身就极难偶然出现）
          或 臂长 ≥20 且相似度 ≥0.90（20 bp 内错配 ≤2）
      中  臂长 ≥16 且相似度 ≥0.85
      弱  其余
    """
    # 长臂还要求实测相似度 ≥0.70：57 条长臂里有 6 条的实测只有 0.35–0.51
    # （接近随机的 0.25），而 identity 字段写着 0.86–0.96 —— 只看长度会把这 6 条
    # 放进"强"。实测分布在此处是断开的（0.51 → 0.71），0.70 落在空洞里。
    if (arm >= 100 and sim >= 0.70) or (arm >= 20 and sim >= 0.90):
        return "强"
    if arm >= 16 and sim >= 0.85:
        return "中"
    return "弱"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="EDTA *.intact.gff3")
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--tsv", default="", help="逐元件输出")
    args = ap.parse_args()

    seq = read_fasta(args.fasta)
    elems = parse_tir_gff(args.gff)
    if not elems:
        sys.exit("没有读到结构 TIR 元件（type 以 _TIR_transposon 结尾的行）")
    npos = sum(len(v) for v in seq.values())

    rows = []
    for e in elems:
        sq = seq.get(e["chrom"])
        if sq is None:
            print(f"跳过 {e['chrom']}：FASTA 里没有这条序列", file=sys.stderr)
            continue
        s, en = e["start"], e["end"]
        lt, rt = e["lt"], e["rt"]
        L = en - s + 1

        left_ok = sq[s - 1:s - 1 + len(lt)] == lt
        right_ok = sq[en - len(rt):en] == rt

        # 反向重复：左翼 vs rc(右翼)。取两臂较短的作比对长度。
        n = min(len(lt), len(rt))
        miss = sum(1 for i in range(n) if lt[i] != rc(rt)[i]) if n else n
        sim = 1 - miss / n if n else 0.0
        log_e = log10_chance(n, miss) + math.log10(npos)

        # TSD：元件外侧两端
        ts_ok = bool(e["tsd_l"]) and sq[max(0, s - 1 - len(e["tsd_l"])):s - 1] == e["tsd_l"]
        te_ok = bool(e["tsd_r"]) and sq[en:en + len(e["tsd_r"])] == e["tsd_r"]

        rows.append({
            "chrom": e["chrom"], "start": s, "end": en, "length": L,
            "cls": e["cls"], "tir_l": len(lt), "tir_r": len(rt),
            "ir_ident_measured": round(sim, 3),
            "ident_field": e["ident"],
            "ir_mismatch": miss,
            "field_gap": round(float(e["ident"]) - sim, 3) + 0.0
                         if e["ident"] not in ("NA", "") else "",
            "chance_log10": round(log_e, 2),
            "tier": tier_of(n, sim),
            "left_tir_ok": int(left_ok), "right_tir_ok": int(right_ok),
            "tsd_len": len(e["tsd_l"]),
            "tsd_left_ok": int(ts_ok), "tsd_right_ok": int(te_ok),
            "tsd_score": e["tsd_score"],
        })

    n = len(rows)
    print(f"结构 TIR 元件: {n}；基因组 {npos:,} bp")
    print(f"两翼 TIR 在基因组上原样复现: 左 {sum(r['left_tir_ok'] for r in rows)}/{n}"
          f"，右 {sum(r['right_tir_ok'] for r in rows)}/{n}")
    print(f"TSD 在元件外侧复现: 左 {sum(r['tsd_left_ok'] for r in rows)}/{n}"
          f"，右 {sum(r['tsd_right_ok'] for r in rows)}/{n}")

    diffs = [abs(r["ir_ident_measured"] - float(r["ident_field"]))
             for r in rows if r["ident_field"] not in ("NA", "")]
    if diffs:
        print(f"逐位相似度 vs identity 字段: 平均 |差| {sum(diffs) / len(diffs):.3f}"
              f"，>0.05 的 {sum(1 for d in diffs if d > 0.05)} 条")
        gap = [r for r in rows if r["field_gap"] != "" and r["field_gap"] > 0.20]
        print(f"  其中字段比实测高 >0.20 的 {len(gap)} 条 —— 这些行的 identity 不可信，"
              f"分级一律以实测为准")
        for r in gap:
            print(f"    {r['chrom']}:{r['start']}-{r['end']} 臂长 {r['tir_l']} "
                  f"实测 {r['ir_ident_measured']:.3f} 字段 {r['ident_field']}")

    print("\n按臂长 + 相似度分级（强：臂≥100 且相似度≥0.70，或 臂≥20 且≥0.90）：")
    tier = Counter(r["tier"] for r in rows)
    for t in ("强", "中", "弱"):
        print(f"  {t}: {tier.get(t, 0):5}  {tier.get(t, 0) / n * 100:5.1f}%")
    print("  注：分级只看反向重复本身。TSD 不作证据位——1,824/1,824 的两臂完全相同，")
    print("      它是 TIR-Learner 的接受条件（v3 只收无 indel 的 TIR/TSD），用它当证据是循环。")
    lg = sorted(r["chance_log10"] for r in rows)
    print(f"  参考：随机命中期望的中位 log10 E = {lg[len(lg) // 2]:.1f}"
          f"（E 越负越难偶然出现；该模型只做单对位比较，未校正搜索全部间距）")

    print("\n强/中档的 TSD 长度（描述性，不作证据）：")
    for t in ("强", "中"):
        sub = [r for r in rows if r["tier"] == t]
        if not sub:
            continue
        c = Counter("≥5 bp" if r["tsd_len"] >= 5 else f"{r['tsd_len']} bp" for r in sub)
        print(f"  {t}（{len(sub)} 条）: {dict(sorted(c.items()))}")

    print("\n按超家族的强度分布：")
    tab = Counter((r["cls"], r["tier"]) for r in rows)
    classes = sorted({r["cls"] for r in rows})
    print("  " + " " * 12 + "".join(f"{t:<7}" for t in ("强", "中", "弱")) + "合计")
    for c in classes:
        v = [tab.get((c, t), 0) for t in ("强", "中", "弱")]
        print(f"  {c:<12}" + "".join(f"{x:<7}" for x in v) + str(sum(v)))

    if args.tsv:
        cols = ["chrom", "start", "end", "length", "cls", "tir_l", "tir_r",
                "ir_ident_measured", "ident_field", "field_gap", "ir_mismatch",
                "chance_log10", "tier", "left_tir_ok", "right_tir_ok",
                "tsd_len", "tsd_left_ok", "tsd_right_ok", "tsd_score"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n逐元件结果: {args.tsv}（{n} 行）")


if __name__ == "__main__":
    main()
