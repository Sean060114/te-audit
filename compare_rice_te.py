#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EDTA 注释 vs MSU RGAP7 官方重复序列注释（水稻 IRGSP-1.0 Chr1，43,270,923 bp）

为什么用 MSU RGAP7 作真值集
    它与 EDTA 无谱系关系，且基于**同一套 IRGSP-1.0 组装**，坐标可直接比对。
    地位相当于拟南芥那边的 Araport11（见 compare_araport11.py）。

两侧的分类体系完全不同，这是本脚本要处理的核心问题
    EDTA : feature type = TEsorter 名（Gypsy_LTR_retrotransposon ...）
           classification = TEsorter 三字母码（LTR/Gypsy、DNA/DTT、MITE/DTT ...）
    MSU  : feature type 一律 repeat_region
           类别藏在 Name 类型码（ORSiTEMT = MITE）与 Note 自由文本
           （"Ty3-gypsy-like retrotransposon, putative"，URL 编码）

    若不做映射直接比类别，得到的"差异"大部分是命名差异，不是生物学差异。
    因此本脚本：
      1) 把所有条目归入一套双桶（细桶 + 三合桶），规则写在 MSU_RULES / ED_MAP，
         逐条可审；
      2) 显式统计**无法映射**的碱基量，作为一切解释的上限；
      3) 结论优先用三合桶（LTR / DNA / 其他），细桶只描述不改判。

指标口径
    所有指标都是区间精确计算，不做"每碱基贴标签、长区间覆盖短区间"的近似：
      - 每类先合并重叠区间（EDTA 允许嵌套，MSU 也允许）
      - 类间重叠用双向游标求交
      - 因此**各类 bp 之和可能大于总并集**，这是注释本身的重叠，不是算错；
        每个指标各自精确，总并集单独给出。

用法
    python compare_rice_te.py [--msu ref/msu/osa1_r7.asm.repeat_masked.gff3.gz]
                              [--edta rice/edta_out/Oryza_sativa.Chr1.fa.mod.EDTA.TEanno.gff3]
                              [--tsv rice/msu_vs_edta.per_class.tsv]
                              [--json rice/msu_vs_edta.json]
"""
import argparse
import gzip
import json
import os
import re
import sys
from bisect import bisect_right
from collections import Counter, defaultdict
from urllib.parse import unquote

GENOME_LEN = 43270923          # IRGSP-1.0 Chr1
MSU_CHROM = "Chr1"
ED_CHROM = "OsChr1"

# EDTA 输出里这些 feature 不是元件本体，参与统计会虚增覆盖
SKIP_ED = {"target_site_duplication", "long_terminal_repeat",
           "primer_binding_site", "repeat_region"}

# ---------------------------------------------------------------- 类别映射
# 细桶。两侧都映射到同一套名字，才能对齐。
B_LTR_GYPSY, B_LTR_COPIA, B_LTR_UNK = "LTR/Gypsy", "LTR/Copia", "LTR/unknown"
B_MITE, B_TIR = "DNA/MITE", "DNA/TIR"
B_HELI, B_LINE, B_SINE = "Helitron", "LINE", "SINE"
B_NONTE, B_UNCLS = "Other/non-TE", "Unclassified"

FINE = [B_LTR_GYPSY, B_LTR_COPIA, B_LTR_UNK, B_MITE, B_TIR,
        B_HELI, B_LINE, B_SINE, B_NONTE, B_UNCLS]

# 三合桶：结论以它为准。MITE 与 TIR 在两侧的边界不一致（见报告），合并后再比更稳。
COARSE = {B_LTR_GYPSY: "LTR", B_LTR_COPIA: "LTR", B_LTR_UNK: "LTR",
          B_MITE: "DNA", B_TIR: "DNA",
          B_HELI: "其他", B_LINE: "其他", B_SINE: "其他",
          B_NONTE: "其他", B_UNCLS: "其他"}

# MSU 的 Note 规则：**顺序敏感，先匹配先归属**。
# 关键教训：先判"是不是 TE"，再判"是不是端粒/着丝粒"——因为
# "Rice chromosome 8 centromere-like LTR retrotransposon" 是 TE，
# 若先匹配 centromere 就会被错划成非 TE。
MSU_RULES = [
    (r"helitron",                                        B_HELI),
    (r"\bmite\b",                                        B_MITE),
    (r"gypsy|ty3",                                       B_LTR_GYPSY),
    (r"copia|ty1",                                       B_LTR_COPIA),
    # SINE 的写法带编号（p-SINE2-r2010、SINE1-4），词尾边界不成立，
    # 只能用词首匹配。早先写 \bsine\b 导致 1832 条被误归进 LTR/unknown 与 DNA/TIR。
    (r"p-sine|\bsine",                                   B_SINE),
    (r"\bline|\bl1\b",                                   B_LINE),
    (r"retrotransposon|\bltr\b",                         B_LTR_UNK),
    (r"transposon|en/spm|anacb|\bmule\b|mutator|"
     r"\bpif\b|harbinger|tc1|mariner|\bhat\b|stowaway",   B_TIR),
    (r"telomere|centromere|satellite|rrna|\bdna\b",       B_NONTE),
]

# EDTA 的 TEsorter 分类 → 细桶。三字母码含义见 TE_Sequence_Ontology。
ED_MAP = {
    "LTR/Gypsy": B_LTR_GYPSY, "LTR/Copia": B_LTR_COPIA, "LTR/unknown": B_LTR_UNK,
    "LTR/TR": B_LTR_UNK, "LTR/ERV1": B_LTR_UNK,
    "MITE/DTT": B_MITE, "MITE/DTH": B_MITE, "MITE/DTM": B_MITE,
    "MITE/DTC": B_MITE, "MITE/DTE": B_MITE, "MITE/DTX": B_MITE,
    "DNA/DTA": B_TIR, "DNA/DTC": B_TIR, "DNA/DTH": B_TIR, "DNA/DTM": B_TIR,
    "DNA/DTT": B_TIR, "DNA/DTE": B_TIR, "DNA/DTR": B_TIR, "DNA/DTX": B_TIR,
    "DNA/Helitron": B_HELI, "RC/Helitron": B_HELI, "Helitron": B_HELI,
    "LINE/L1": B_LINE, "LINE/RTE": B_LINE, "LINE/I": B_LINE,
    "SINE/tRNA": B_SINE, "SINE/5S": B_SINE,
}
# feature type 兜底（classification 缺失时用）
FT_MAP = {
    "helitron": B_HELI,
    "Gypsy_LTR_retrotransposon": B_LTR_GYPSY,
    "Copia_LTR_retrotransposon": B_LTR_COPIA,
    "LTR_retrotransposon": B_LTR_UNK,
    "L1_LINE_retrotransposon": B_LINE,
    "RTE_LINE_retrotransposon": B_LINE,
    "SINE_retrotransposon": B_SINE,
    "repeat_fragment": B_UNCLS,
}


def bucket_msu(note):
    t = note.lower()
    for pat, b in MSU_RULES:
        if re.search(pat, t):
            return b
    return B_UNCLS


def bucket_ed(cls_attr, ftype):
    if cls_attr in ED_MAP:
        return ED_MAP[cls_attr]
    if ftype in FT_MAP:
        return FT_MAP[ftype]
    if ftype.endswith("_TIR_transposon"):
        return B_TIR
    if ftype.endswith("_LTR_retrotransposon"):
        return B_LTR_UNK
    if ftype.endswith("_LINE_retrotransposon"):
        return B_LINE
    return B_UNCLS


# ---------------------------------------------------------------- 区间工具
def merge(iv):
    """合并重叠与相邻区间，返回升序不重叠列表"""
    if not iv:
        return []
    iv.sort()
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def span(iv):
    return sum(e - s for s, e in iv)


def inter_len(a, b):
    """两个已合并升序区间表的交集总长（双向游标）"""
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def covered_len(s, e, mg, starts):
    """[s,e) 被合并区间表 mg 覆盖的碱基数"""
    tot = 0
    k = bisect_right(starts, s) - 1
    if k < 0:
        k = 0
    while k < len(mg) and mg[k][0] < e:
        a, b = mg[k]
        ss, ee = max(a, s), min(b, e)
        if ee > ss:
            tot += ee - ss
        k += 1
    return tot


def clamp(a, b):
    a, b = max(0, a), min(GENOME_LEN, b)
    return (a, b) if b > a else None


# ---------------------------------------------------------------- 解析
def parse_msu(path, chrom=MSU_CHROM):
    """返回 fine_bucket -> [(s,e)], 以及 Note 原始分布的 bp 统计"""
    groups = defaultdict(list)
    note_bp = Counter()
    note_bucket = Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            m = re.search(r"Note=([^;]+)", p[8])
            note = unquote(m.group(1)) if m else ""
            iv = clamp(int(p[3]) - 1, int(p[4]))
            if not iv:
                continue
            b = bucket_msu(note)
            groups[b].append(iv)
            ln = iv[1] - iv[0]
            note_bp[note] += ln
            note_bucket[(b, note)] += ln
    return groups, note_bp, note_bucket


def parse_ed(path, chrom=ED_CHROM):
    groups = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom or p[2] in SKIP_ED:
                continue
            m = re.search(r"classification=([^;]+)", p[8])
            cls = m.group(1) if m else ""
            iv = clamp(int(p[3]) - 1, int(p[4]))
            if not iv:
                continue
            groups[bucket_ed(cls, p[2])].append(iv)
    return groups


# ---------------------------------------------------------------- 输出
def sec(t):
    return "\n" + "-" * 86 + "\n" + t + "\n" + "-" * 86


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msu", default="ref/msu/osa1_r7.asm.repeat_masked.gff3.gz")
    ap.add_argument("--edta",
                    default="rice/edta_out/Oryza_sativa.Chr1.fa.mod.EDTA.TEanno.gff3")
    ap.add_argument("--tsv", default="rice/msu_vs_edta.per_class.tsv")
    ap.add_argument("--json", default="rice/msu_vs_edta.json")
    a = ap.parse_args()

    for p in (a.msu, a.edta):
        if not os.path.exists(p):
            sys.exit("找不到: " + p)

    msu_gr, note_bp, note_bucket = parse_msu(a.msu)
    ed_gr = parse_ed(a.edta)

    msu_mg = {k: merge(v) for k, v in msu_gr.items()}
    ed_mg = {k: merge(v) for k, v in ed_gr.items()}
    msu_all = merge([iv for v in msu_gr.values() for iv in v])
    ed_all = merge([iv for v in ed_gr.values() for iv in v])
    msu_bp, ed_bp = span(msu_all), span(ed_all)
    inter = inter_len(msu_all, ed_all)

    out = []
    add = out.append
    L = GENOME_LEN
    add("=" * 86)
    add(f"EDTA vs MSU RGAP7 — 水稻 IRGSP-1.0 Chr1, {L:,} bp")
    add("=" * 86)
    add("真值集：osa1_r7.asm.repeat_masked.gff3（MSU/UGA，与 EDTA 无谱系关系，同一组装）")

    # ---------- 一、总体
    add(sec("[一] 总体覆盖度（区间精确，非按碱基贴标签）"))
    add(f"  全染色体                          {L:>14,} bp")
    add(f"  MSU 重复注释（合并后）            {msu_bp:>14,} bp  {msu_bp / L * 100:>6.2f}%")
    add(f"  EDTA 注释（合并后）               {ed_bp:>14,} bp  {ed_bp / L * 100:>6.2f}%")
    add(f"  两者交集                          {inter:>14,} bp")
    add(f"  -> EDTA 对 MSU 的碱基召回         {inter / msu_bp * 100:>13.1f}%")
    add(f"  -> EDTA 的碱基精确率              {inter / ed_bp * 100:>13.1f}%")
    add(f"  -> EDTA 独有（MSU 未注释）        {ed_bp - inter:>14,} bp"
        f"  ({(ed_bp - inter) / ed_bp * 100:.1f}% of EDTA)")
    add(f"  -> MSU 独有（EDTA 未注释）        {msu_bp - inter:>14,} bp"
        f"  ({(msu_bp - inter) / msu_bp * 100:.1f}% of MSU)")

    # ---------- 二、映射体检
    add(sec("[二] 类别映射体检（这一节决定后面能解释到什么程度）"))
    add(f"  Chr1 上 MSU 共 {len(note_bp)} 种不同 Note 文本，按规则归入以下桶：")
    add(f"  {'桶':<16}{'条目bp':>14}{'不同Note数':>11}")
    per_bucket_notes = Counter()
    for (b, _n), bp in note_bucket.items():
        per_bucket_notes[b] += 1
    for b in FINE:
        bp = sum(v for (bb, _n), v in note_bucket.items() if bb == b)
        if bp:
            add(f"  {b:<16}{bp:>14,}{per_bucket_notes[b]:>11}")
    unm = [n for n in note_bp if bucket_msu(n) == B_UNCLS]
    add("")
    add(f"  落入 Unclassified 的不同 Note 数：{len(unm)}"
        f"（占 MSU bp 的 "
        f"{sum(v for (b, _n), v in note_bucket.items() if b == B_UNCLS) / msu_bp * 100:.1f}%）")
    add("  说明：Unclassified 包含 Note 为空、'Unclassified repetitive sequence'、")
    add("        以及规则未覆盖的少数文本。它是解释上限，不是 EDTA 的漏检。")

    # ---------- 三、细桶：召回
    add(sec("[三] 细桶召回：MSU 每类碱基有多少被 EDTA 覆盖，以及 EDTA 把它判成什么"))
    add(f"  {'MSU桶':<16}{'bp':>13}{'被EDTA覆盖':>13}{'召回':>8}   EDTA 主判（占比）")
    add("  " + "-" * 82)
    rows = []
    for b in FINE:
        if b not in msu_mg:
            continue
        bp = span(msu_mg[b])
        cov = inter_len(msu_mg[b], ed_all)
        # EDTA 对该桶的贡献排序
        contrib = []
        for eb, emg in ed_mg.items():
            x = inter_len(msu_mg[b], emg)
            if x:
                contrib.append((x, eb))
        contrib.sort(reverse=True)
        top = "；".join(f"{eb} {x / bp * 100:.0f}%" for x, eb in contrib[:2]) or "—"
        add(f"  {b:<16}{bp:>13,}{cov:>13,}{cov / bp * 100:>7.1f}%   {top}")
        rows.append(dict(bucket=b, side="MSU", bp=bp, covered=cov,
                         rate=cov / bp * 100, top=top))

    # ---------- 四、细桶：精确率
    add(sec("[四] 细桶精确率：EDTA 每类碱基有多少被 MSU 支持"))
    add(f"  {'EDTA桶':<16}{'bp':>13}{'被MSU支持':>13}{'精确率':>8}   同桶支持")
    add("  " + "-" * 82)
    for b in FINE:
        if b not in ed_mg:
            continue
        bp = span(ed_mg[b])
        sup = inter_len(ed_mg[b], msu_all)
        same = inter_len(ed_mg[b], msu_mg.get(b, []))
        add(f"  {b:<16}{bp:>13,}{sup:>13,}{sup / bp * 100:>7.1f}%   "
            f"{same / bp * 100:>6.1f}%")
        rows.append(dict(bucket=b, side="EDTA", bp=bp, covered=sup,
                         rate=sup / bp * 100, same=same / bp * 100))

    # ---------- 五、三合桶
    add(sec("[五] 三合桶（结论以此为准：MITE 与 TIR 在两侧边界不一致，合并后再比）"))
    msu_coarse, ed_coarse = defaultdict(list), defaultdict(list)
    for b, mg in msu_mg.items():
        msu_coarse[COARSE[b]].extend(mg)
    for b, mg in ed_mg.items():
        ed_coarse[COARSE[b]].extend(mg)
    msu_coarse = {k: merge(v) for k, v in msu_coarse.items()}
    ed_coarse = {k: merge(v) for k, v in ed_coarse.items()}
    add(f"  {'组':<8}{'MSU bp':>13}{'EDTA bp':>13}{'交集':>13}{'召回':>8}{'精确率':>9}")
    add("  " + "-" * 78)
    for g in ("LTR", "DNA", "其他"):
        m, e = msu_coarse.get(g, []), ed_coarse.get(g, [])
        mb, eb = span(m), span(e)
        x = inter_len(m, e)
        r = x / mb * 100 if mb else float("nan")
        p = x / eb * 100 if eb else float("nan")
        add(f"  {g:<8}{mb:>13,}{eb:>13,}{x:>13,}{r:>7.1f}%{p:>8.1f}%")
        rows.append(dict(bucket=g, side="coarse", msu_bp=mb, edta_bp=eb,
                         inter=x, recall=r, precision=p))

    # ---------- 六、元件级
    add(sec("[六] 元件级召回：MSU 每条重复序列被 EDTA 覆盖 ≥50% 的比例"))
    add(f"  {'MSU桶':<16}{'元件数':>9}{'≥50%':>9}{'召回':>8}{'≥10%':>9}")
    add("  " + "-" * 76)
    ed_starts = [s for s, _ in ed_all]
    for b in FINE:
        ivs = msu_gr.get(b, [])
        if not ivs:
            continue
        n50 = n10 = 0
        for s, e in ivs:
            c = covered_len(s, e, ed_all, ed_starts)
            frac = c / (e - s)
            n50 += frac >= 0.5
            n10 += frac >= 0.1
        add(f"  {b:<16}{len(ivs):>9,}{n50:>9,}{n50 / len(ivs) * 100:>7.1f}%"
            f"{n10:>9,}")
        rows.append(dict(bucket=b, side="MSU_elements", n=len(ivs),
                         ge50=n50, recall_el=n50 / len(ivs) * 100, ge10=n10))

    # ---------- 七、EDTA 独有
    add(sec("[七] EDTA 独有：与 MSU 零重叠的 EDTA 元件（候选新发现，也是精确率的主要扣分项）"))
    add(f"  {'EDTA桶':<16}{'元件数':>9}{'零重叠':>9}{'占比':>8}{'零重叠bp':>13}")
    add("  " + "-" * 76)
    for b, mg in ed_mg.items():
        ivs = ed_gr[b]
        zero = [iv for iv in ivs if inter_len([iv], msu_all) == 0]
        if not ivs:
            continue
        add(f"  {b:<16}{len(ivs):>9,}{len(zero):>9,}"
            f"{len(zero) / len(ivs) * 100:>7.1f}%{sum(e - s for s, e in zero):>13,}")
        rows.append(dict(bucket=b, side="EDTA_only", n=len(ivs),
                         nov=len(zero), nov_pct=len(zero) / len(ivs) * 100))

    # ---------- 八、独有部分的构成
    add(sec("[八] 双方独有碱基的构成（解释精确率/召回率的关键，不能用单一数字概括）"))
    add("  EDTA 独有 %s bp 按 EDTA 类别拆开：" % f"{ed_bp - inter:,}")
    add(f"  {'EDTA桶':<16}{'未被MSU覆盖bp':>15}{'占EDTA独有':>11}")
    add("  " + "-" * 74)
    for b in sorted(ed_mg, key=lambda x: -(span(ed_mg[x]) - inter_len(ed_mg[x], msu_all))):
        unc = span(ed_mg[b]) - inter_len(ed_mg[b], msu_all)
        add(f"  {b:<16}{unc:>15,}{unc / (ed_bp - inter) * 100:>10.1f}%")
        rows.append(dict(bucket=b, side="EDTA_uncovered", bp=unc,
                         pct=unc / (ed_bp - inter) * 100))
    add("")
    add("  MSU 独有 %s bp 按 MSU 类别拆开：" % f"{msu_bp - inter:,}")
    add(f"  {'MSU桶':<16}{'未被EDTA覆盖bp':>15}{'占MSU独有':>11}")
    add("  " + "-" * 74)
    for b in sorted(msu_mg, key=lambda x: -(span(msu_mg[x]) - inter_len(msu_mg[x], ed_all))):
        unc = span(msu_mg[b]) - inter_len(msu_mg[b], ed_all)
        if not unc:
            continue
        add(f"  {b:<16}{unc:>15,}{unc / (msu_bp - inter) * 100:>10.1f}%")
        rows.append(dict(bucket=b, side="MSU_uncovered", bp=unc,
                         pct=unc / (msu_bp - inter) * 100))

    # ---------- 写出
    print("\n".join(out))
    with open(a.tsv, "w", encoding="utf-8", newline="\n") as f:
        keys = ["side", "bucket", "bp", "covered", "rate"]
        f.write("\t".join(keys + ["extra"]) + "\n")
        for r in rows:
            extra = json.dumps({k: v for k, v in r.items() if k not in keys},
                               ensure_ascii=False)
            f.write("\t".join(str(r.get(k, "")) for k in keys) + "\t" + extra + "\n")
    print(f"\n[写出] {a.tsv}")

    with open(a.json, "w", encoding="utf-8", newline="\n") as f:
        json.dump(dict(genome_len=L, msu_bp=msu_bp, ed_bp=ed_bp, inter=inter,
                       rows=rows,
                       note_bp_top=note_bp.most_common(30)), f,
                  ensure_ascii=False, indent=1)
    print(f"[写出] {a.json}")


if __name__ == "__main__":
    main()
