#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EDTA 注释 vs MSU RGAP7 官方重复序列注释（水稻 IRGSP-1.0，多染色体）

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

多染色体（v2）
    不同染色体的坐标空间互相重叠，**绝不能把区间跨染色体合并**。
    因此本版改为「逐染色体独立计算，再按碱基量汇总」：
      - 每个染色体各自 merge / 求交，互不干扰
      - 汇总表只累加 bp 与元件数，不碰坐标
      - --fasta 提供每条染色体的真实长度（缺失时退化为 MSU 最大坐标）

指标口径
    所有指标都是区间精确计算，不做"每碱基贴标签、长区间覆盖短区间"的近似：
      - 每类先合并重叠区间（EDTA 允许嵌套，MSU 也允许）
      - 类间重叠用双向游标求交
      - 因此**各类 bp 之和可能大于总并集**，这是注释本身的重叠，不是算错；
        每个指标各自精确，总并集单独给出。

用法
    # 单染色体（与旧版等价）
    python compare_rice_te.py --edta rice/edta_out/Oryza_sativa.Chr1.fa.mod.EDTA.TEanno.gff3

    # 多染色体（Chr1+Chr12 同库运行）
    python compare_rice_te.py --edta rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.TEanno.gff3 \
        --fasta rice/Oryza_sativa.Chr1_Chr12.fa \
        --chrom OsChr1=Chr1 --chrom OsChr12=Chr12 \
        --tsv rice/msu_vs_edta.per_class.tsv --json rice/msu_vs_edta.json
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
    iv = sorted(iv)
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


# ---------------------------------------------------------------- 解析
def parse_msu(path, wanted):
    """wanted: MSU 染色体名集合。返回 chrom -> (groups, note_bp, note_bucket)"""
    groups = {c: defaultdict(list) for c in wanted}
    note_bp = {c: Counter() for c in wanted}
    note_bucket = {c: Counter() for c in wanted}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in wanted:
                continue
            c = p[0]
            m = re.search(r"Note=([^;]+)", p[8])
            note = unquote(m.group(1)) if m else ""
            s, e = int(p[3]) - 1, int(p[4])
            if e <= s:
                continue
            b = bucket_msu(note)
            groups[c][b].append((s, e))
            note_bp[c][note] += e - s
            note_bucket[c][(b, note)] += e - s
    return groups, note_bp, note_bucket


def parse_ed(path, wanted):
    """wanted: EDTA seqid 集合"""
    groups = {c: defaultdict(list) for c in wanted}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in wanted or p[2] in SKIP_ED:
                continue
            m = re.search(r"classification=([^;]+)", p[8])
            cls = m.group(1) if m else ""
            s, e = int(p[3]) - 1, int(p[4])
            if e <= s:
                continue
            groups[p[0]][bucket_ed(cls, p[2])].append((s, e))
    return groups


def fasta_lengths(path):
    """返回 seqid -> 长度"""
    out, name, n = {}, None, 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    out[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name:
        out[name] = n
    return out


# ---------------------------------------------------------------- 单染色体分析
def analyze(chrom_ed, chrom_msu, L, msu_gr, ed_gr, note_bp, note_bucket):
    """一条染色体上的全部指标。返回 dict（含 rows 与可直接打印的段落文本）"""
    msu_mg = {k: merge(v) for k, v in msu_gr.items()}
    ed_mg = {k: merge(v) for k, v in ed_gr.items()}
    msu_all = merge([iv for v in msu_gr.values() for iv in v])
    ed_all = merge([iv for v in ed_gr.values() for iv in v])
    msu_bp, ed_bp = span(msu_all), span(ed_all)
    inter = inter_len(msu_all, ed_all)

    lines, rows = [], []

    def add(t=""):
        lines.append(t)

    add("=" * 86)
    add(f"{chrom_ed} / {chrom_msu}   EDTA vs MSU RGAP7 — {L:,} bp")
    add("=" * 86)

    # ---- 总体
    add("\n" + "-" * 86 + "\n[总体] 覆盖度（区间精确，非按碱基贴标签）\n" + "-" * 86)
    add(f"  染色体长度                        {L:>14,} bp")
    add(f"  MSU 重复注释（合并后）            {msu_bp:>14,} bp  {msu_bp / L * 100:>6.2f}%")
    add(f"  EDTA 注释（合并后）               {ed_bp:>14,} bp  {ed_bp / L * 100:>6.2f}%")
    add(f"  两者交集                          {inter:>14,} bp")
    add(f"  -> EDTA 对 MSU 的碱基召回         {inter / msu_bp * 100:>13.1f}%")
    add(f"  -> EDTA 的碱基精确率              {inter / ed_bp * 100:>13.1f}%")
    add(f"  -> EDTA 独有（MSU 未注释）        {ed_bp - inter:>14,} bp"
        f"  ({(ed_bp - inter) / ed_bp * 100:.1f}% of EDTA)")
    add(f"  -> MSU 独有（EDTA 未注释）        {msu_bp - inter:>14,} bp"
        f"  ({(msu_bp - inter) / msu_bp * 100:.1f}% of MSU)")
    rows.append(dict(chrom=chrom_ed, side="overall", bucket="ALL",
                     bp=None, covered=None, rate=None,
                     L=L, msu_bp=msu_bp, ed_bp=ed_bp, inter=inter,
                     recall=inter / msu_bp * 100,
                     precision=inter / ed_bp * 100))

    # ---- 映射体检
    add("\n" + "-" * 86 + "\n[映射体检] 这一节决定后面能解释到什么程度\n" + "-" * 86)
    unm_bp = sum(v for (b, _n), v in note_bucket.items() if b == B_UNCLS)
    unm_notes = len({n for n in note_bp if bucket_msu(n) == B_UNCLS})
    add(f"  该染色体 MSU 共 {len(note_bp)} 种不同 Note 文本，{unm_notes} 种落入 Unclassified")
    add(f"  Unclassified 占 MSU bp 的 {unm_bp / msu_bp * 100:.1f}%"
        f"（{unm_bp:,} bp）")
    add("  说明：Unclassified 含 Note 为空、'Unclassified repetitive sequence'、")
    add("        以及规则未覆盖的少数文本。它是解释上限，不是 EDTA 的漏检。")

    # ---- 细桶召回
    add("\n" + "-" * 86 + "\n[细桶召回] MSU 每类碱基有多少被 EDTA 覆盖，以及 EDTA 把它判成什么\n" + "-" * 86)
    add(f"  {'MSU桶':<16}{'bp':>13}{'被EDTA覆盖':>13}{'召回':>8}   EDTA 主判（占比）")
    add("  " + "-" * 82)
    for b in FINE:
        if b not in msu_mg:
            continue
        bp = span(msu_mg[b])
        cov = inter_len(msu_mg[b], ed_all)
        contrib = []
        for eb, emg in ed_mg.items():
            x = inter_len(msu_mg[b], emg)
            if x:
                contrib.append((x, eb))
        contrib.sort(reverse=True)
        top = "；".join(f"{eb} {x / bp * 100:.0f}%" for x, eb in contrib[:2]) or "—"
        add(f"  {b:<16}{bp:>13,}{cov:>13,}{cov / bp * 100:>7.1f}%   {top}")
        rows.append(dict(chrom=chrom_ed, side="MSU", bucket=b, bp=bp, covered=cov,
                         rate=cov / bp * 100, top=top))

    # ---- 细桶精确率
    add("\n" + "-" * 86 + "\n[细桶精确率] EDTA 每类碱基有多少被 MSU 支持\n" + "-" * 86)
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
        rows.append(dict(chrom=chrom_ed, side="EDTA", bucket=b, bp=bp, covered=sup,
                         rate=sup / bp * 100, same=same / bp * 100))

    # ---- 三合桶
    add("\n" + "-" * 86 + "\n[三合桶] 结论以此为准（MITE 与 TIR 在两侧边界不一致，合并后再比）\n" + "-" * 86)
    msu_coarse, ed_coarse = defaultdict(list), defaultdict(list)
    for b, mg in msu_mg.items():
        msu_coarse[COARSE[b]].extend(mg)
    for b, mg in ed_mg.items():
        ed_coarse[COARSE[b]].extend(mg)
    msu_coarse = {k: merge(v) for k, v in msu_coarse.items()}
    ed_coarse = {k: merge(v) for k, v in ed_coarse.items()}
    add(f"  {'组':<8}{'MSU bp':>13}{'EDTA bp':>13}{'交集':>13}{'召回':>8}{'精确率':>9}")
    add("  " + "-" * 78)
    coarse_out = {}
    for g in ("LTR", "DNA", "其他"):
        m, e = msu_coarse.get(g, []), ed_coarse.get(g, [])
        mb, eb = span(m), span(e)
        x = inter_len(m, e)
        r = x / mb * 100 if mb else float("nan")
        p = x / eb * 100 if eb else float("nan")
        add(f"  {g:<8}{mb:>13,}{eb:>13,}{x:>13,}{r:>7.1f}%{p:>8.1f}%")
        coarse_out[g] = dict(msu_bp=mb, edta_bp=eb, inter=x, recall=r, precision=p)
        rows.append(dict(chrom=chrom_ed, side="coarse", bucket=g, bp=eb,
                         covered=x, rate=p, msu_bp=mb, edta_bp=eb,
                         inter=x, recall=r, precision=p))

    # ---- 元件级
    add("\n" + "-" * 86 + "\n[元件级召回] MSU 每条重复序列被 EDTA 覆盖 ≥50% 的比例\n" + "-" * 86)
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
        rows.append(dict(chrom=chrom_ed, side="MSU_elements", bucket=b,
                         bp=None, covered=None, rate=n50 / len(ivs) * 100,
                         n=len(ivs), ge50=n50, ge10=n10))

    # ---- 独有构成
    add("\n" + "-" * 86 + "\n[独有构成] 解释精确率/召回率的关键，不能用单一数字概括\n" + "-" * 86)
    ed_only = {}
    add(f"  EDTA 独有 {ed_bp - inter:,} bp 按 EDTA 类别拆开：")
    add(f"  {'EDTA桶':<16}{'未被MSU覆盖bp':>15}{'占EDTA独有':>11}")
    add("  " + "-" * 74)
    for b in sorted(ed_mg, key=lambda x: -(span(ed_mg[x]) - inter_len(ed_mg[x], msu_all))):
        unc = span(ed_mg[b]) - inter_len(ed_mg[b], msu_all)
        pct = unc / (ed_bp - inter) * 100 if ed_bp > inter else 0.0
        add(f"  {b:<16}{unc:>15,}{pct:>10.1f}%")
        ed_only[b] = dict(bp=unc, pct=pct)
        rows.append(dict(chrom=chrom_ed, side="EDTA_uncovered", bucket=b,
                         bp=unc, covered=None, rate=pct))
    msu_only = {}
    add("")
    add(f"  MSU 独有 {msu_bp - inter:,} bp 按 MSU 类别拆开：")
    add(f"  {'MSU桶':<16}{'未被EDTA覆盖bp':>15}{'占MSU独有':>11}")
    add("  " + "-" * 74)
    for b in sorted(msu_mg, key=lambda x: -(span(msu_mg[x]) - inter_len(msu_mg[x], ed_all))):
        unc = span(msu_mg[b]) - inter_len(msu_mg[b], ed_all)
        msu_only[b] = dict(bp=unc, pct=unc / (msu_bp - inter) * 100 if msu_bp > inter else 0.0)
        if not unc:
            continue
        add(f"  {b:<16}{unc:>15,}{unc / (msu_bp - inter) * 100:>10.1f}%")
        rows.append(dict(chrom=chrom_ed, side="MSU_uncovered", bucket=b,
                         bp=unc, covered=None, rate=unc / (msu_bp - inter) * 100))

    return dict(chrom=chrom_ed, msu_chrom=chrom_msu, L=L,
                msu_bp=msu_bp, ed_bp=ed_bp, inter=inter,
                recall=inter / msu_bp * 100, precision=inter / ed_bp * 100,
                ed_only_bp=ed_bp - inter, msu_only_bp=msu_bp - inter,
                unclassified_pct=unm_bp / msu_bp * 100,
                coarse=coarse_out, ed_only=ed_only, msu_only=msu_only,
                lines=lines, rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msu", default="ref/msu/osa1_r7.asm.repeat_masked.gff3.gz")
    ap.add_argument("--edta",
                    default="rice/edta_out/Oryza_sativa.Chr1.fa.mod.EDTA.TEanno.gff3")
    ap.add_argument("--fasta", default=None,
                    help="基因组 FASTA，用于取每条序列的真实长度（推荐给出）")
    ap.add_argument("--chrom", action="append", default=None,
                    metavar="ED=MSU",
                    help="EDTA seqid=MSU 染色体名，可重复。缺省为 OsChr1=Chr1")
    ap.add_argument("--tsv", default="rice/msu_vs_edta.per_class.tsv")
    ap.add_argument("--json", default="rice/msu_vs_edta.json")
    a = ap.parse_args()

    for p in (a.msu, a.edta):
        if not os.path.exists(p):
            sys.exit("找不到: " + p)

    pairs = [("OsChr1", "Chr1")]
    if a.chrom:
        pairs = []
        for spec in a.chrom:
            if "=" not in spec:
                sys.exit(f"--chrom 需要 ED=MSU 形式，收到 {spec!r}")
            ed, msu = spec.split("=", 1)
            pairs.append((ed.strip(), msu.strip()))

    # 长度：优先 FASTA，其次退化为 MSU 最大坐标
    lengths = {}
    if a.fasta:
        if not os.path.exists(a.fasta):
            sys.exit("找不到 FASTA: " + a.fasta)
        lengths = fasta_lengths(a.fasta)
        miss = [ed for ed, _ in pairs if ed not in lengths]
        if miss:
            sys.exit(f"FASTA 里没有这些序列：{miss}；实际有 {sorted(lengths)}")

    wanted_msu = {m for _, m in pairs}
    wanted_ed = {e for e, _ in pairs}
    msu_g, note_bp, note_bucket = parse_msu(a.msu, wanted_msu)
    ed_g = parse_ed(a.edta, wanted_ed)

    results, all_rows = [], []
    for ed_c, msu_c in pairs:
        L = lengths.get(ed_c)
        if L is None:
            # 退化：用 MSU 该染色体最大坐标（会略短于真实长度，故告警）
            mx = max((e for v in msu_g[msu_c].values() for _, e in v), default=0)
            L = mx
            print(f"[警告] {ed_c} 未提供 FASTA 长度，用 MSU 最大坐标 {L:,} 代替"
                  f"（真实长度略大，百分比会偏高）", file=sys.stderr)
        if not msu_g[msu_c]:
            print(f"[警告] MSU 中没有 {msu_c} 的记录，跳过", file=sys.stderr)
            continue
        if not ed_g[ed_c]:
            print(f"[警告] EDTA 中没有 {ed_c} 的记录，跳过", file=sys.stderr)
            continue
        r = analyze(ed_c, msu_c, L, msu_g[msu_c], ed_g[ed_c],
                    note_bp[msu_c], note_bucket[msu_c])
        results.append(r)
        all_rows.extend(r["rows"])
        print("\n".join(r["lines"]))

    if len(results) > 1:
        out = ["", "=" * 86,
               "跨染色体汇总（不同染色体坐标互不重叠，此处只累加 bp，不做区间合并）",
               "=" * 86]
        tot_L = sum(r["L"] for r in results)
        tot_m = sum(r["msu_bp"] for r in results)
        tot_e = sum(r["ed_bp"] for r in results)
        tot_x = sum(r["inter"] for r in results)
        out.append(f"  {'染色体':<10}{'长度':>12}{'MSU%':>8}{'EDTA%':>8}"
                   f"{'召回':>8}{'精确率':>8}{'Uncl%':>8}")
        out.append("  " + "-" * 70)
        for r in results:
            out.append(f"  {r['chrom']:<10}{r['L']:>12,}"
                       f"{r['msu_bp'] / r['L'] * 100:>7.2f}%"
                       f"{r['ed_bp'] / r['L'] * 100:>7.2f}%"
                       f"{r['recall']:>7.1f}%{r['precision']:>7.1f}%"
                       f"{r['unclassified_pct']:>7.1f}%")
        out.append("  " + "-" * 70)
        out.append(f"  {'合计':<10}{tot_L:>12,}"
                   f"{tot_m / tot_L * 100:>7.2f}%{tot_e / tot_L * 100:>7.2f}%"
                   f"{tot_x / tot_m * 100:>7.1f}%{tot_x / tot_e * 100:>7.1f}%")

        out.append("")
        out.append("  三合桶逐染色体（LTR 组是唯一两侧口径可比的粒度）：")
        out.append(f"  {'染色体':<10}{'组':<6}{'MSU bp':>13}{'EDTA bp':>13}"
                   f"{'交集':>13}{'召回':>8}{'精确率':>9}")
        out.append("  " + "-" * 78)
        for r in results:
            for g in ("LTR", "DNA", "其他"):
                c = r["coarse"][g]
                out.append(f"  {r['chrom']:<10}{g:<6}{c['msu_bp']:>13,}"
                           f"{c['edta_bp']:>13,}{c['inter']:>13,}"
                           f"{c['recall']:>7.1f}%{c['precision']:>8.1f}%")
        out.append("")
        out.append("  EDTA 独有部分的类别构成（判断精确率扣分来自哪里）：")
        hdr = "  " + f"{'染色体':<10}" + "".join(f"{b.split('/')[-1][:11]:>12}" for b in FINE)
        out.append(hdr)
        out.append("  " + "-" * (10 + 12 * len(FINE)))
        for r in results:
            line = f"  {r['chrom']:<10}"
            for b in FINE:
                p = r["ed_only"].get(b, {}).get("pct", 0.0)
                line += f"{p:>11.1f}%" if p >= 0.05 else f"{'—':>12}"
            out.append(line)
        print("\n".join(out))
        print()

    with open(a.tsv, "w", encoding="utf-8", newline="\n") as f:
        keys = ["chrom", "side", "bucket", "bp", "covered", "rate"]
        f.write("\t".join(keys + ["extra"]) + "\n")
        for r in all_rows:
            extra = json.dumps({k: v for k, v in r.items() if k not in keys},
                               ensure_ascii=False)
            f.write("\t".join("" if r.get(k) is None else str(r.get(k))
                              for k in keys) + "\t" + extra + "\n")
    print(f"[写出] {a.tsv}")

    payload = dict(chroms=[{k: v for k, v in r.items()
                            if k not in ("lines", "rows")} for r in results],
                   rows=all_rows)
    with open(a.json, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    print(f"[写出] {a.json}")


if __name__ == "__main__":
    main()
