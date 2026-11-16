#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EDTA 注释 vs Araport11 官方 TE 注释（TAIR10 Chr1，30,427,671 bp）

两套命名差异：
  Araport11 : classification=LTR/Copia/ATCOPIA24   (class/subclass/family)
  EDTA      : classification=DNA/DTC ...           (TEsorter 三字母码)
              且元件类别写在 feature type 列（helitron / Gypsy_LTR_retrotransposon ...）

输出：
  1. 总体覆盖度
  2. 类级召回（Araport 各类碱基被 EDTA 覆盖的比例 + 被 EDTA 判成什么）
  3. 类级精确率（EDTA 各类碱基被 Araport 支持的比例）
  4. 亚类混淆矩阵（bp 级）
  5. 元件级召回（Araport 元件被 EDTA 覆盖 ≥50% 的比例）
  6. EDTA 中与 Araport 零重叠的元件
"""
import os
import re
import sys
from collections import defaultdict

import numpy as np

GENOME_LEN = 30427671
AP_GFF = "araport11_TE.gff3"
ED_GFF = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod.EDTA.TEanno.gff3"
SKIP_ED = {"target_site_duplication", "long_terminal_repeat", "primer_binding_site"}

# TEsorter 三字母码 -> 可读名
AAA = {"DTA": "hAT", "DTC": "CACTA", "DTH": "PIF_Harbinger",
       "DTM": "Mutator", "DTT": "Tc1_Mariner", "DTE": "hAT",
       "DTR": "hAT", "DTX": "Crypton", "RIT": "RIT", "Ngaro": "Ngaro"}


def split_class(c):
    """返回 (class, subclass)"""
    if not c:
        return "Unknown", ""
    parts = [p.strip() for p in c.split("/") if p.strip()]
    if not parts:
        return "Unknown", ""
    top = parts[0]
    low = c.lower()
    if "helitron" in low or top in ("RC", "Helitron"):
        return "RC", "Helitron"
    if top == "MITE":
        return "DNA", AAA.get(parts[1], parts[1]) if len(parts) > 1 else "MITE"
    if top in ("LTR", "DNA", "LINE", "SINE"):
        sub = parts[1] if len(parts) > 1 else ""
        return top, AAA.get(sub, sub)
    if low in ("unknown", "other"):
        return "Unknown", ""
    return "Other", ""


def merge(iv):
    """合并左闭右开区间"""
    if not iv:
        return []
    iv.sort()
    out = [(iv[0][0], iv[0][1])]
    for s, e in iv[1:]:
        if s > out[-1][1]:
            out.append((s, e))
        elif e > out[-1][1]:
            out[-1] = (out[-1][0], e)
    return out


def parse_ap(path, chrom="Chr1"):
    """返回 {class: [(s,e)]}, {subclass: [(s,e)]}"""
    bc, bs = defaultdict(list), defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            if p[2] != "transposable_element":
                # transposon_fragment 经核验 100% 嵌套在 transposable_element 内，
                # 且多数无 classification，纳入会污染类别统计
                continue
            m = re.search(r"classification=([^;]+)", p[8])
            c, s = split_class(m.group(1) if m else "")
            a, b = int(p[3]) - 1, int(p[4])
            a, b = max(0, a), min(GENOME_LEN, b)
            if b <= a:
                continue
            bc[c].append((a, b))
            if s:
                bs[f"{c}/{s}"].append((a, b))
    return bc, bs


def parse_ed(path, chrom="1"):
    bc, bs = defaultdict(list), defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom or p[2] in SKIP_ED:
                continue
            m = re.search(r"classification=([^;]+)", p[8])
            c, s = split_class(m.group(1) if m else "")
            a, b = int(p[3]) - 1, int(p[4])
            a, b = max(0, a), min(GENOME_LEN, b)
            if b <= a:
                continue
            bc[c].append((a, b))
            if s:
                bs[f"{c}/{s}"].append((a, b))
    return bc, bs


def label_array(groups):
    """把 {label: intervals} 画成每碱基标签数组（长区间优先），返回 (arr, code->label)"""
    arr = np.zeros(GENOME_LEN, dtype=np.int16)
    mapping = {0: "none"}
    items = [(k, v) for k, v in groups.items() if v]
    # 按区间长度升序写入，长区间后写 → 覆盖短区间
    flat = []
    for k, v in items:
        for s, e in v:
            flat.append((e - s, k, s, e))
    flat.sort()
    labels = [k for k, _ in items]
    lab2code = {k: i + 1 for i, k in enumerate(labels)}
    for k, code in lab2code.items():
        mapping[code] = k
    for ln, k, s, e in flat:
        arr[s:e] = lab2code[k]
    return arr, mapping


def bincount2(a, b, K):
    """统计 (a,b) 组合计数，返回二维矩阵"""
    m = (a > 0) & (b > 0)
    idx = a[m].astype(np.int64) * K + b[m].astype(np.int64)
    cnt = np.bincount(idx, minlength=K * K)
    return cnt.reshape(K, K)


def main():
    for p in (AP_GFF, ED_GFF):
        if not os.path.exists(p):
            sys.exit("找不到: " + p)

    ap_c, ap_s = parse_ap(AP_GFF)
    ed_c, ed_s = parse_ed(ED_GFF)

    ap_all = np.zeros(GENOME_LEN, dtype=bool)
    for v in ap_c.values():
        for s, e in merge(v):
            ap_all[s:e] = True
    ed_all = np.zeros(GENOME_LEN, dtype=bool)
    for v in ed_c.values():
        for s, e in merge(v):
            ed_all[s:e] = True

    # 预先构建 EDTA 各类掩码（类数很少）
    ed_masks = {}
    for c, iv in ed_c.items():
        m = np.zeros(GENOME_LEN, dtype=bool)
        for s, e in merge(iv):
            m[s:e] = True
        ed_masks[c] = m

    L = GENOME_LEN
    o = []
    add = o.append
    add("=" * 84)
    add(f"EDTA vs Araport11 - TAIR10 Chr1, {L:,} bp")
    add("=" * 84)

    add("\n[一] 总体覆盖度")
    add(f"  Araport11 TE        {ap_all.sum():>12,} bp  {ap_all.sum()/L*100:>6.2f}%")
    add(f"  EDTA                {ed_all.sum():>12,} bp  {ed_all.sum()/L*100:>6.2f}%")
    inter = int((ap_all & ed_all).sum())
    add(f"  两者交集            {inter:>12,} bp")
    add(f"  -> EDTA 对 Araport 的碱基召回 : {inter/ap_all.sum()*100:5.1f}%")
    add(f"  -> EDTA 的碱基精确率          : {inter/ed_all.sum()*100:5.1f}%")
    add(f"  -> EDTA 独有(Araport 未注)    : {int((ed_all&~ap_all).sum()):>12,} bp "
        f"({(ed_all&~ap_all).sum()/ed_all.sum()*100:.1f}% of EDTA)")

    # ---------- 类级 ----------
    add("\n[二] 类级召回：Araport11 每一类的碱基有多少被 EDTA 覆盖")
    add(f"  {'class':<12}{'元件':>7}{'bp':>12}{'被覆盖bp':>12}{'召回':>8}   EDTA 给的类别")
    add("  " + "-" * 80)
    rows = []
    for c, iv in ap_c.items():
        mg = merge(iv)
        bp = sum(e - s for s, e in mg)
        if bp == 0:
            continue
        tmp = np.zeros(L, dtype=bool)
        for s, e in mg:
            tmp[s:e] = True
        cov = int((tmp & ed_all).sum())
        mix = {k: int((tmp & v).sum()) for k, v in ed_masks.items()}
        mix = {k: v for k, v in mix.items() if v}
        tot = sum(mix.values())
        top = sorted(mix.items(), key=lambda x: -x[1])[:3]
        desc = ", ".join(f"{k} {v/tot*100:.0f}%" for k, v in top)
        rows.append((c, len(iv), bp, cov, cov / bp, desc))
    rows.sort(key=lambda r: -r[2])
    for c, n, bp, cov, r, desc in rows:
        add(f"  {c:<12}{n:>7}{bp:>12,}{cov:>12,}{r*100:>7.1f}%   {desc or '-'}")

    add("\n[三] 类级精确率：EDTA 每一类的碱基有多少落在 Araport11 TE 内")
    add(f"  {'class':<12}{'元件':>7}{'bp':>12}{'被支持bp':>12}{'精确率':>8}   同类别支持")
    add("  " + "-" * 80)
    er = []
    for c, iv in ed_c.items():
        mg = merge(iv)
        bp = sum(e - s for s, e in mg)
        if bp == 0:
            continue
        tmp = np.zeros(L, dtype=bool)
        for s, e in mg:
            tmp[s:e] = True
        hit = int((tmp & ap_all).sum())
        same = 0
        if c in ap_c:
            apc = np.zeros(L, dtype=bool)
            for s, e in merge(ap_c[c]):
                apc[s:e] = True
            same = int((tmp & apc).sum())
        er.append((c, len(iv), bp, hit, hit / bp, same / bp))
    er.sort(key=lambda r: -r[2])
    for c, n, bp, hit, r, same in er:
        add(f"  {c:<12}{n:>7}{bp:>12,}{hit:>12,}{r*100:>7.1f}%   {same*100:>6.1f}%")

    # ---------- 亚类混淆 ----------
    add("\n[四] 亚类混淆矩阵（bp；行为 Araport11，列为 EDTA，仅列主要亚类）")
    ap_arr, ap_map = label_array(ap_s)
    ed_arr, ed_map = label_array(ed_s)
    K = max(len(ap_map), len(ed_map)) + 1
    M = bincount2(ap_arr, ed_arr, K)
    ap_names = dict(ap_map)
    ed_names = dict(ed_map)
    ap_tot = M.sum(axis=1)
    ed_tot = M.sum(axis=0)
    keep_r = [i for i in range(1, K) if ap_tot[i] > 20000 and i in ap_names]
    keep_r.sort(key=lambda i: -ap_tot[i])
    keep_c = [j for j in range(1, K) if ed_tot[j] > 20000 and j in ed_names]
    keep_c.sort(key=lambda j: -ed_tot[j])
    hdr = "".ljust(26) + "".join(f"{ed_names[j][:13]:>14}" for j in keep_c) + f"{'未注':>12}"
    add(hdr)
    for i in keep_r:
        line = f"{ap_names[i][:24]:<26}"
        covered = 0
        for j in keep_c:
            v = M[i, j]
            covered += v
            line += f"{v:>14,}"
        line += f"{ap_tot[i]-covered:>12,}"
        add(line)

    # ---------- 元件级召回 ----------
    add("\n[五] 元件级召回（Araport11 元件被 EDTA 覆盖 ≥50% 视为“找到”）")
    cum = np.concatenate(([0], np.cumsum(ed_all, dtype=np.int64)))
    add(f"  {'class':<12}{'元件数':>8}{'≥50%覆盖':>10}{'召回':>8}{'≥10%覆盖':>10}")
    add("  " + "-" * 56)
    for c, iv in sorted(ap_c.items(), key=lambda x: -len(x[1])):
        n = g50 = g10 = 0
        for s, e in iv:
            ln = e - s
            if ln <= 0:
                continue
            cov = int(cum[e] - cum[s])
            n += 1
            if cov / ln >= 0.5:
                g50 += 1
            if cov / ln >= 0.1:
                g10 += 1
        if n:
            add(f"  {c:<12}{n:>8}{g50:>10}{g50/n*100:>7.1f}%{g10:>10}")

    # ---------- 零重叠元件 ----------
    add("\n[六] EDTA 中与 Araport11 零重叠的元件")
    n0 = bp0 = 0
    per = defaultdict(lambda: [0, 0])
    for c, iv in ed_c.items():
        for s, e in iv:
            if not ap_all[s:e].any():
                n0 += 1
                bp0 += e - s
                per[c][0] += 1
                per[c][1] += e - s
    tot_n = sum(len(v) for v in ed_c.values())
    add(f"  合计 {n0} / {tot_n} 个元件，{bp0:,} bp ({bp0/L*100:.2f}% 基因组)")
    for c, (n, bp) in sorted(per.items(), key=lambda x: -x[1][1]):
        add(f"    {c:<12}{n:>6} 个   {bp:>10,} bp")

    txt = "\n".join(o)
    print(txt)
    with open("araport11_vs_edta.txt", "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    print("\n[写出] araport11_vs_edta.txt")

    # 类级明细表（机器可读，供交付包与投稿附表使用）
    rows = []
    for src, d, other in (("Araport11", ap_c, ed_all), ("EDTA", ed_c, ap_all)):
        for cls, iv in d.items():
            m = np.zeros(L, dtype=bool)
            for s, e in iv:
                m[s:e] = True
            bp = int(m.sum())
            cov = int((m & other).sum())
            rows.append((src, cls, len(iv), bp, cov, cov / bp if bp else 0.0))
    order = {"Araport11": 0, "EDTA": 1}
    rows.sort(key=lambda r: (order[r[0]], -r[3]))
    with open("araport11_vs_edta.per_class.tsv", "w",
              encoding="utf-8", newline="\n") as f:
        f.write("source\tclass\tn_elements\tbp\tcovered_bp\trate\t"
                "rate_meaning\n")
        for src, cls, n, bp, cov, rate in rows:
            f.write(f"{src}\t{cls}\t{n}\t{bp}\t{cov}\t{rate:.4f}\t"
                    f"{'recall' if src == 'Araport11' else 'precision'}\n")
    print("[写出] araport11_vs_edta.per_class.tsv")


if __name__ == "__main__":
    main()
