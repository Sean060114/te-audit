#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐元件台账：每个 EDTA 结构元件，在 MSU 侧有没有对应的类别/超家族判断。

两类元件，一套判定（--kind）
---------------------------
  ltr  LTR 反转录转座子   repeat_region 行，细到 Gypsy / Copia
  tir  TIR/MITE 转座子    *_TIR_transposon 行，细到 DTA/DTC/DTH/DTM/DTT

与 compare_rice_te.py 的分工
----------------------------
那份脚本算的是**染色体级 bp 指标**（并集密度、召回、精确率），没有"元件↔元件"
的对应关系。这里要的是**逐条台账**：这一个元件，MSU 说它是什么、覆盖了多少。

复用（import 自 compare_rice_te）：
  bucket_msu    MSU 的 Note → 类别。顺序敏感的规则链，SINE 那个坑已修过。
  bucket_ed     EDTA 的 classification → 类别。
  merge         区间并集。
  span          区间表总长。

不复用：
  parse_ed      它的 SKIP_ED 含 repeat_region —— 那正是元件本体，一条都读不到；
                且不区分 method=structural|homology，会把同源碎片混进来。
  parse_msu     它把 Note 按染色体聚合成 bp，元件级拿不到"哪条 Note 压在它身上"。
  analyze       整体是 bp 级并集/交集，没有元件身份。

三条口径上的硬约束（都是踩过的坑）
----------------------------------
1. **覆盖率必须用并集**：MSU 是碎片化注释（Chr1 就有 57k 行），按"单行覆盖 ≥50%"
   判会有 207/420 元件被误判成"无对应"。实测并集口径下 10%/30%/50% 三档稳定。

2. **没有 Note 的 MSU 行不是 TE**。MSU 的 repeat_masked.gff3 里 22% 的行（Chr1
   12,761/57,134）只有 Name 没有 Note，Name 全部是 `(CGG)n` / `AT_rich` /
   `SUBTEL_sat` 这类简单重复与卫星——162 种形态无一例外。`bucket_msu("")` 把它们
   归进 Unclassified，于是任何一个跨在微卫星上的元件都可能在 "MSU 判为非 TE"
   里记一笔。**实测：30% 档下 62 条 TIR 锚的判定一条都不变**（这些行都不到 30%
   覆盖），10% 档下会有 1 条从"非 TE"翻成类别判定——所以这条修的是语义，不是
   数字，但语义错了迟早会在别的样本上变成数字错。

3. **"Unclassified transposon, putative" 不能并进 DNA/TIR**。它占 DNA/TIR 桶的
   4,009 行（该桶最大单项），可它恰恰是 MSU 也说不清超家族的那批——把它算作
   "MSU 说是 TIR"，等于拿"不知道"当"知道"。单列成 DNA/unknown，语义等于 LTR 的
   LTR/unknown（只到类层）。

为什么 TIR 只能给到"类级"结论
------------------------------
LTR 侧的比较停在超家族（Gypsy vs Copia），因为 MSU 的 Note 里写着 gypsy/ty3、
copia/ty1。TIR 侧 MSU 写着 AnacB12、MITE-adh、Tourist-like 这些**家族名**，其中
能可靠折回 DT 超家族的只有 En/Spm、CACTA、MULE/Mutator 三类（Tourist、Castaway、
Ditto、Gaijin、Crackle、Wanderer 都是家族名，换算过去是猜）。所以：
  - 主结论用**类级**（Class II DNA transposon）：MSU 判 DNA 转座子就算一致；
  - 超家族一致率只在"两源都给出可映射超家族"的子集上算，n 单独报，不作结论。

用法
----
  python TE_case/link_msu.py --kind ltr \
    --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
    --msu ref/msu/osa1_r7.asm.repeat_masked.gff3.gz \
    --tsv TE_case/ltr_msu_ledger.tsv

  python TE_case/link_msu.py --kind tir \
    --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
    --msu ref/msu/osa1_r7.asm.repeat_masked.gff3.gz \
    --anchor-tsv TE_case/tir_structure.tsv --tier 强,中 \
    --tsv TE_case/tir_msu_ledger.tsv
"""
import argparse
import csv
import gzip
import os
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from urllib.parse import unquote

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import compare_rice_te as C          # noqa: E402  只取映射与区间工具


# --- 控制台编码：Windows 默认 cp936，⚠ / → 编不出来会抛 UnicodeEncodeError，
#     报告打到一半就断；而 Git Bash 的终端是 UTF-8，按 cp936 写又是乱码。
#     按终端分两种处理，都不允许因编码中断（不可编码字符降级为 ?）。
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8" if os.environ.get("MSYSTEM") else None,
                       errors="replace")
    except (AttributeError, ValueError):
        pass

# 台账里写 ASCII 码、控制台里印中文：数据文件要能被下游直接读（本仓库既有 TSV
# 一律无 BOM、无中文表头），而人看的是屏幕。
VERDICT_LTR = {
    "agree": "agree_superfamily",
    "conflict": "conflict_superfamily",
    "coarse": "msu_ltr_only",
    "ed_unk": "ed_unknown_msu_resolved",
    "both_unk": "both_unresolved",
    "nonte": "msu_nonltr_label",
    # 细分：MSU 在这块只有简单重复/低复杂度，根本没给 TE 注释。与"判成了别的
    # 转座子"是两件事，混在一个码里会让人以为 MSU 反着标了。
    "nonte_only": "msu_nonte_label",
    "insuff": "insufficient_cov",
    "empty": "empty",
}
VERDICT_TIR = {
    "agree": "agree_superfamily",
    "agree_class": "agree_class",
    "conflict": "conflict_class",
    "coarse": "msu_class_only",
    "ed_unk": "ed_unknown_msu_resolved",
    "both_unk": "both_unresolved",
    "nonte": "msu_nonte_label",
    "insuff": "insufficient_cov",
    "empty": "empty",
}

LTR_RESOLVED = (C.B_LTR_GYPSY, C.B_LTR_COPIA)
LTR_ALL = LTR_RESOLVED + (C.B_LTR_UNK,)

# MSU 只到"转座子"这一层的通用 Note。match 从头锚定，不会误伤 AnacB12 这类带
# 家族名的写法。占该桶 4,009/约 13k 行。
B_DNA_UNK = "DNA/unknown"
GENERIC_TIR = re.compile(r"^(unclassified )?transposon\b")

# 只收**可靠**的超家族换算。en/spm、cacta 同属 DTC（CACTA）；MULE/Mutator 属 DTM。
# 其余（Tourist、Castaway、Ditto、Gaijin、Crackle、Wanderer、MITE-adh、AnacA/B/C）
# 是家族名，不是超家族，换算过去是猜——不放进这张表，宁可让 n 变小。
MSU_SUPERFAMILY = [
    (r"en/spm", "DTC"),
    (r"cacta|cac1|hipa", "DTC"),
    (r"\bmule\b|mutator", "DTM"),
]


# --------------------------------------------------------------------- 读入
def read_elements(gff, kind):
    """读结构元件。一个元件 = 一行。

    LTR 一个元件在 intact.gff3 里占 6 行（repeat_region / lTSD / lLTR / LTRRT /
    rLTR / rTSD），只有 repeat_region 代表元件本体，且带 tsd / motif /
    ltr_identity 三个字段。TIR 元件只有一行，字段写在行内（TIR= / TSD= /
    identity=）。

    身份用 (chrom, start, end) 而不是 Name —— Name 是**家族名**，
    TE_00000794 在 Chr1+Chr12 上出现 9 次。
    """
    out = []
    with open(gff, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if kind == "ltr" and ftype != "repeat_region":
                continue
            if kind == "tir" and not ftype.endswith("_TIR_transposon"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            cls = a.get("classification", "")
            if kind == "ltr" and not cls.startswith("LTR"):
                continue
            # TIR 的字段名是大写（TSD= / TIR=），LTR 是小写（tsd=）；两处都取
            tsd = a.get("tsd") or a.get("TSD") or "NA"
            if kind == "tir" and "_" in tsd:
                # TSD=左_右_打分 —— 第三段是打分不是长度，见 基准_TIR字段语义校准.md
                tsd = tsd.split("_")[0]
            out.append({
                # gff 的 1-based 闭区间原样保留（人工核对时看的就是这两个数）；
                # 所有区间算术在 _coverage 里转成 0-based 半开，与 MSU 侧一致。
                "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                "strand": f[6], "ftype": ftype, "cls": cls,
                "ident": a.get("ltr_identity") or a.get("identity") or "NA",
                "tsd": tsd, "motif": a.get("motif", "NA"),
                "name": a.get("Name", ""),      # 家族名，只作参考，不当主键
            })
    return out


def msu_bucket(note):
    """MSU Note → 桶。在 compare_rice_te.bucket_msu 之上补两条口径（见文件头）。"""
    if not note:
        return C.B_NONTE            # 简单重复/低复杂度/卫星，不是 TE
    b = C.bucket_msu(note)
    if b == C.B_TIR and GENERIC_TIR.match(note.lower()):
        return B_DNA_UNK            # MSU 也只到"转座子"，没说超家族
    return b


def msu_superfamily(note):
    t = note.lower()
    for pat, sf in MSU_SUPERFAMILY:
        if re.search(pat, t):
            return sf
    return ""


def ed_superfamily(cls):
    """"DNA/DTM" / "MITE/DTC" → DTM / DTC。"""
    return cls.split("/")[-1] if "/" in cls else ""


def read_msu(path, cmap):
    """返回 EDTA seqid -> (starts, rows)，rows 为 (start0, end, bucket, note)。

    MSU 的 chr 名与 EDTA 的 seqid 不同（Chr1 vs OsChr1）。cmap 是
    MSU 名 → EDTA seqid；不在 cmap 里的染色体直接丢弃——**不要先归一再用
    原名去过滤**，那样会把所有行都过滤掉（本脚本初版就这样，MSU 行数 0）。
    """
    rows = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            c = cmap.get(p[0])
            if c is None:
                continue
            m = re.search(r"Note=([^;]+)", p[8])
            note = unquote(m.group(1)) if m else ""
            s, e = int(p[3]) - 1, int(p[4])
            if e <= s:
                continue
            rows[c].append((s, e, msu_bucket(note), note))
    out = {}
    for c, v in rows.items():
        v.sort(key=lambda x: x[0])
        out[c] = ([r[0] for r in v], v)
    return out


# --------------------------------------------------------------------- 区间
def _coverage(elem, msu):
    """返回 (L, bucket_cov, note_cov, noted, n_hits)。

    bucket_cov  桶级并集覆盖的碱基数（同一桶内相邻/重叠行只算一次）
    note_cov    (桶, Note) 级并集 —— 判超家族时必需：一个桶里可能压着十条不同
                Note，取"桶覆盖率"再配"最常见的那条 Note"会把少数派意见放大。
    """
    s0 = elem["start"] - 1                       # 1-based 闭 → 0-based 半开
    L = elem["end"] - s0
    starts, rows = msu.get(elem["chrom"], ([], []))
    if not rows:
        return L, {}, {}, defaultdict(Counter), 0
    maxlen = max(e - s for s, e, _, _ in rows)
    i = bisect_left(starts, s0 - maxlen)
    per_bucket = defaultdict(list)
    per_note = defaultdict(list)
    noted = defaultdict(Counter)
    n = 0
    while i < len(rows) and rows[i][0] < elem["end"]:
        s, e, b, note = rows[i]
        if e > s0:
            n += 1
            iv = (max(s, s0), min(e, elem["end"]))
            per_bucket[b].append(iv)
            per_note[(b, note)].append(iv)
            noted[b][note] += 1
        i += 1
    bucket_cov = {b: C.span(C.merge(iv)) for b, iv in per_bucket.items()}
    note_cov = {k: C.span(C.merge(iv)) for k, iv in per_note.items()}
    return L, bucket_cov, note_cov, noted, n


def _boundary(elem, msu, min_share):
    """最贴合的一条 MSU 行（与元件交集最大）的边界差，1-based 口径。

    返回 (ds, de, note)：ds = MSU起点 − EDTA起点，de = MSU终点 − EDTA终点。
    这是检验"两侧坐标约定是否一致"最直接的一列——LTR 侧因为 EDTA 把 TSD 算进
    元件、MSU 不算，众数恰好是 (+5,−5)；TIR 侧两侧都不含 TSD，应当更贴。
    """
    s0 = elem["start"] - 1
    starts, rows = msu.get(elem["chrom"], ([], []))
    if not rows:
        return "", "", ""
    maxlen = max(e - s for s, e, _, _ in rows)
    i = bisect_left(starts, s0 - maxlen)
    best = None
    while i < len(rows) and rows[i][0] < elem["end"]:
        s, e, b, note = rows[i]
        i += 1
        if b == C.B_NONTE:
            continue
        ov = min(e, elem["end"]) - max(s, s0)
        if ov >= min_share * (elem["end"] - s0) and (best is None or ov > best[0]):
            best = (ov, s, e, note)
    if best is None:
        return "", "", ""
    _, s, e, note = best
    return (s + 1) - elem["start"], e - elem["end"], note


def _top_note(noted, bucket):
    if not bucket or bucket not in noted or not noted[bucket]:
        return ""
    return noted[bucket].most_common(1)[0][0]


# --------------------------------------------------------------------- 判定
def classify_ltr(elem, msu, min_cov):
    """返回 (verdict_key, dominant_bucket, bucket_cov, n_hits, top_note)。

    判定顺序（**先认一致，再认冲突，最后才降级为"没解出来"**）：
      一致 > 冲突 > MSU 只到 LTR 层 > MSU 判为非 LTR > 覆盖不足
    顺序不能反：若先判 "MSU/unknown 达标" 再判冲突，一个 Gypsy 元件同时压着
    Gypsy 与 unknown 两条时会漏掉真实的冲突。
    """
    eth = C.bucket_ed(elem["cls"], elem["ftype"])
    L, real, _nc, noted, n = _coverage(elem, msu)
    if L <= 0:
        return "empty", "", {}, 0, ""

    ok = {b: v for b, v in real.items() if v >= min_cov * L}
    if not ok:
        dom = max(real, key=real.get) if real else ""
        return "insuff", dom, real, n, _top_note(noted, dom)

    ed_unres = eth == C.B_LTR_UNK or eth not in LTR_RESOLVED
    others = [b for b in ok if b in LTR_RESOLVED and b != eth]

    if eth in LTR_RESOLVED and eth in ok:
        return "agree", eth, real, n, _top_note(noted, eth)
    if others and not ed_unres:
        # 两源都给出了超家族，且不同 —— 这才是真冲突
        dom = max(others, key=lambda b: ok[b])
        return "conflict", dom, real, n, _top_note(noted, dom)
    if ed_unres and others:
        # EDTA 自己没定（LTR/unknown），MSU 定了 —— 不是冲突，是单边可判
        dom = max(others, key=lambda b: ok[b])
        return "ed_unk", dom, real, n, _top_note(noted, dom)
    if C.B_LTR_UNK in ok and not ed_unres:
        return "coarse", C.B_LTR_UNK, real, n, _top_note(noted, C.B_LTR_UNK)
    if C.B_LTR_UNK in ok:
        return "both_unk", C.B_LTR_UNK, real, n, _top_note(noted, C.B_LTR_UNK)
    if C.B_NONTE in ok and not any(b in LTR_ALL for b in ok):
        # MSU 在这里只认到简单重复/低复杂度 —— 不是"判成了别的转座子"，
        # 是"这块位置 MSU 没给 TE 注释"。两者不能混为一谈。
        return "nonte_only", C.B_NONTE, real, n, _top_note(noted, C.B_NONTE)

    dom = max(ok, key=ok.get)
    return "nonte", dom, real, n, _top_note(noted, dom)


CLASSII = (C.B_TIR, C.B_MITE, B_DNA_UNK)
CLASSI = (C.B_LTR_GYPSY, C.B_LTR_COPIA, C.B_LTR_UNK, C.B_LINE, C.B_SINE)


def classify_tir(elem, msu, min_cov):
    """TIR 版判定，层级比 LTR 低一档（类级为主，超家族为辅）。

    判定顺序同样是 一致 > 冲突 > 降级：
      agree            两源都给出可映射的超家族且相同（DTC/DTC、DTM/DTM…）
      agree_class      两源都说 Class II DNA 转座子，但 MSU 只给家族名/没给超家族
      conflict         两源都定了，MSU 说它是 Class I（LTR/LINE/SINE）或 Helitron
      coarse           MSU 只到"Unclassified transposon"这一层
      nonte            MSU 在该位置只有简单重复/低复杂度
      insuff           覆盖不足

    agree 与 agree_class 分开报，是因为前者能当结论、后者只能当"不矛盾"。
    """
    eth = C.bucket_ed(elem["cls"], elem["ftype"])
    L, real, note_cov, noted, n = _coverage(elem, msu)
    if L <= 0:
        return "empty", "", "", real, 0, ""

    ok = {b: v for b, v in real.items() if v >= min_cov * L}
    if not ok:
        dom = max(real, key=real.get) if real else ""
        return "insuff", dom, "", real, n, _top_note(noted, dom)

    ed_sf = ed_superfamily(elem["cls"])
    dna = {b: v for b, v in ok.items() if b in (C.B_TIR, C.B_MITE)}
    ci = {b: v for b, v in ok.items() if b in CLASSI or b == C.B_HELI}

    # --- 超家族层：找"映射到同一 DT 码"的 Note，看它自己覆盖够不够
    if ed_sf and dna:
        hit_cov = 0
        hit_note = ""
        hit_b = ""
        for (b, note), v in note_cov.items():
            if b in (C.B_TIR, C.B_MITE) and msu_superfamily(note) == ed_sf and v > hit_cov:
                hit_cov, hit_note, hit_b = v, note, b
        if hit_cov >= min_cov * L:
            return "agree", hit_b, ed_sf, real, n, hit_note
    # MSU 给出超家族、但与 EDTA 不同 —— 真冲突（例：EDTA 说 DTM，MSU 说 CACTA）
    if ed_sf and dna:
        rival = [(v, note, b) for (b, note), v in note_cov.items()
                 if b in (C.B_TIR, C.B_MITE) and v >= min_cov * L
                 and msu_superfamily(note) and msu_superfamily(note) != ed_sf]
        if rival:
            v, note, b = max(rival)
            return "conflict", b, msu_superfamily(note), real, n, note

    if dna:
        dom = max(dna, key=dna.get)
        return "agree_class", dom, "", real, n, _top_note(noted, dom)
    if ci:
        dom = max(ci, key=ci.get)
        return "conflict", dom, "", real, n, _top_note(noted, dom)
    if B_DNA_UNK in ok:
        return "coarse", B_DNA_UNK, "", real, n, _top_note(noted, B_DNA_UNK)
    if C.B_NONTE in ok:
        return "nonte", C.B_NONTE, "", real, n, _top_note(noted, C.B_NONTE)
    if C.B_LTR_UNK in ok:
        return "both_unk", C.B_LTR_UNK, "", real, n, _top_note(noted, C.B_LTR_UNK)

    dom = max(ok, key=ok.get)
    return "conflict", dom, "", real, n, _top_note(noted, dom)


# --------------------------------------------------------------------- 主流程
COV_COLS = [("cov_gypsy", C.B_LTR_GYPSY), ("cov_copia", C.B_LTR_COPIA),
            ("cov_ltr_unknown", C.B_LTR_UNK), ("cov_tir", C.B_TIR),
            ("cov_mite", C.B_MITE), ("cov_dna_unknown", B_DNA_UNK),
            ("cov_sine", C.B_SINE), ("cov_line", C.B_LINE),
            ("cov_helitron", C.B_HELI), ("cov_nonte", C.B_NONTE)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="EDTA *.intact.gff3")
    ap.add_argument("--msu", required=True, help="MSU repeat_masked.gff3(.gz)")
    ap.add_argument("--kind", choices=("ltr", "tir"), default="ltr")
    ap.add_argument("--min-cov", type=float, default=0.3,
                    help="判定「有对应」所需的最小并集覆盖率（默认 0.3）")
    ap.add_argument("--anchor-tsv", default="",
                    help="含 tier 列的分级表（如 TE_case/tir_structure.tsv），按 "
                         "(chrom,start,end) 关联；只影响 tier 列与 --tier 过滤")
    ap.add_argument("--tier", default="",
                    help="只保留这些档，逗号分隔，如 强,中")
    ap.add_argument("--tsv", default="", help="台账输出路径")
    ap.add_argument("--chrom", action="append", default=[],
                    help="EDTA seqid=MSU chr 的映射，如 --chrom OsChr1=Chr1；"
                         "缺省按去掉前缀 Os 处理")
    args = ap.parse_args()
    for f in (args.gff, args.msu):
        if not os.path.isfile(f):
            sys.exit(f"找不到输入文件：{f}\n"
                     f"仓库根目录下应能访问 rice/…intact.gff3 与 ref/msu/…gff3.gz；"
                     f"序列与参考数据不入库，换机器需先放回原位")

    elems = read_elements(args.gff, args.kind)
    if not elems:
        what = "LTR 元件（repeat_region 行，classification 以 LTR 开头）" \
            if args.kind == "ltr" else "TIR 元件（*_TIR_transposon 行）"
        sys.exit(f"没有读到 {what}——检查 gff 是否为 *.intact.gff3")

    tiers = {}
    if args.anchor_tsv:
        with open(args.anchor_tsv, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tiers[(r["chrom"], int(r["start"]), int(r["end"]))] = r.get("tier", "")
    want = {t.strip() for t in args.tier.split(",") if t.strip()}
    if want:
        before = len(elems)
        elems = [e for e in elems
                 if tiers.get((e["chrom"], e["start"], e["end"]), "") in want]
        print(f"按档过滤 {sorted(want)}：{before} → {len(elems)} 个元件")
        if not elems:
            sys.exit("过滤后为空——检查 anchor-tsv 的坐标口径是否为 1-based 闭区间")

    cmap = {}
    for e in elems:
        ed = e["chrom"]
        pair = [x for x in args.chrom if x.split("=")[0] == ed]
        cmap[ed] = pair[0].split("=")[1] if pair else ed[2:] if ed.startswith("Os") else ed
    # cmap 是 EDTA seqid → MSU 名；read_msu 要的是反向表
    rev = {ms: ed for ed, ms in cmap.items()}
    msu = read_msu(args.msu, rev)

    classify = classify_ltr if args.kind == "ltr" else classify_tir
    codes = VERDICT_LTR if args.kind == "ltr" else VERDICT_TIR

    ledger = []
    over = 0
    for e in elems:
        res = classify(e, msu, args.min_cov)
        if args.kind == "ltr":
            key, dom, cov, n, note = res
            msu_sf = ""
        else:
            key, dom, msu_sf, cov, n, note = res
        L = e["end"] - e["start"] + 1
        if any(v > L for v in cov.values()):
            # 并集覆盖率不可能超过元件长度。超了就是同桶内区间被重复计入，
            # 是分区间的老毛病，宁可中止也不要写出一个 >1 的数。
            over += 1
        ds, de, _bn = _boundary(e, msu, 0.5)
        ds9, de9, _bn9 = _boundary(e, msu, 0.9)
        # "混合"标记：同一元件上既有 Class II 桶又有 Class I/Helitron 桶过阈值。
        # 判定顺序是"一致 > 冲突"，所以这种情况会被判成一致，冲突被盖住——必须
        # 单独暴露出来，否则一致率里混着双标位点。
        okb = {b for b, v in cov.items() if v >= args.min_cov * L}
        mixed = int(any(b in CLASSII for b in okb)
                    and any(b in CLASSI or b == C.B_HELI for b in okb))
        ledger.append({
            "chrom": e["chrom"], "start": e["start"], "end": e["end"],
            "strand": e["strand"], "length": L, "ed_class": e["cls"],
            "ed_superfamily": ed_superfamily(e["cls"]),
            "tier": tiers.get((e["chrom"], e["start"], e["end"]), ""),
            "identity": e["ident"], "tsd": e["tsd"],
            "verdict_code": codes[key], "msu_dominant": dom,
            "msu_superfamily": msu_sf, "mixed": mixed,
            "n_msu_hits": n, "msu_note": note,
            "bd_ds": ds, "bd_de": de, "bd90_ds": ds9, "bd90_de": de9,
            "family": e["name"],
            **{c: round(cov.get(b, 0) / L, 3) for c, b in COV_COLS},
        })

    if over:
        sys.exit(f"{over} 个元件的某桶覆盖率超过 100%——并集没生效，结果不可用")

    # ---- 汇总
    kind_cn = "LTR 元件" if args.kind == "ltr" else "TIR/MITE 元件"
    print(f"\nEDTA 结构 {kind_cn}: {len(ledger)}；"
          f"MSU 行: {sum(len(v[1]) for v in msu.values())}")
    print(f"判定阈值：并集覆盖率 ≥ {args.min_cov:.0%}\n")
    order = list(dict.fromkeys(codes.values()))
    tab = Counter((r["ed_class"], r["verdict_code"]) for r in ledger)
    classes = sorted({r["ed_class"] for r in ledger})
    w = max(len(v) for v in order)
    print(" " * 14 + "".join(f"{v:<{w + 3}}" for v in order) + "合计")
    for cl in classes:
        row = [tab.get((cl, v), 0) for v in order]
        print(f"  {cl:<12}" + "".join(f"{x:<{w + 3}}" for x in row) + f"{sum(row)}")
    print(f"  {'合计':<12}"
          + "".join(f"{sum(1 for r in ledger if r['verdict_code'] == v):<{w + 3}}"
                    for v in order) + f"{len(ledger)}")

    # ---- 率
    n = len(ledger)
    tot = Counter(r["verdict_code"] for r in ledger)
    nmix = sum(1 for r in ledger if r["mixed"])
    print()
    if nmix:
        print(f"双标位点：{nmix} 条同时被 MSU 标成 Class II 桶与 Class I/Helitron 桶"
              f"（判定按「一致优先」取了前者，真实冲突数应理解为 ≥"
              f"{tot.get('conflict_superfamily', 0) + tot.get('conflict_class', 0)}）")
    if args.kind == "ltr":
        both = tot["agree_superfamily"] + tot["conflict_superfamily"]
        if both:
            print(f"两源都定超家族: {both} 条，一致 {tot['agree_superfamily']} = "
                  f"{tot['agree_superfamily'] / both * 100:.1f}%")
    else:
        # 类级分母：MSU 确实在这块位置给了"这是 Class II 转座子"或"这是别的东西"
        # 的判断。coarse（只说 Unclassified transposon）与 nonte（只有简单重复）
        # 都不算 MSU 表了态 —— 前者是"不知道"，后者是"没注释"。
        judged = (tot["agree_superfamily"] + tot["agree_class"]
                  + tot["conflict_class"])
        if judged:
            agree = tot["agree_superfamily"] + tot["agree_class"]
            print(f"类级（Class II DNA 转座子）可判 {judged} 条，一致 {agree} = "
                  f"{agree / judged * 100:.1f}%  "
                  f"（超家族层 {tot['agree_superfamily']} + 仅类层 {tot['agree_class']}）")
            print(f"  其中冲突 {tot['conflict_class']} 条")
        print(f"类级一致率（分母取全部 {n} 条）: "
              f"{(tot['agree_superfamily'] + tot['agree_class']) / n * 100:.1f}%"
              f"  —— 把 MSU 未表态的 {tot['msu_class_only'] + tot['msu_nonte'] + tot['insufficient_cov']} 条"
              f"按不一致计，故不作结论")
        print(f"MSU 未表态: 只到类层 {tot['msu_class_only']}，"
              f"仅简单重复 {tot['msu_nonte']}，覆盖不足 {tot['insufficient_cov']}")

        # ---- 超家族层（小样本，只报 n）
        sf = [r for r in ledger if r["msu_superfamily"]]
        print(f"\n超家族层：MSU 给出可映射超家族的 {len(sf)} 条")
        if sf:
            print("  " + " " * 10 + "MSU→" + "".join(f"{s:<7}" for s in ("DTA", "DTC", "DTH", "DTM", "DTT")))
            for ed in ("DTA", "DTC", "DTH", "DTM", "DTT"):
                sub = [r for r in sf if r["ed_superfamily"] == ed]
                if not sub:
                    continue
                print(f"  EDTA {ed:<6}"
                      + "".join(f"{sum(1 for r in sub if r['msu_superfamily'] == s):<7}"
                                for s in ("DTA", "DTC", "DTH", "DTM", "DTT"))
                      + f"{len(sub)}")

    # ---- 边界
    # 两档都报：≥90% 是"MSU 有一条几乎盖住整个元件"的干净子集（LTR 侧 81.1% 那个数
    # 就出自这一档），≥50% 是宽口径。差额全部落在"MSU 只盖了半个元件"的那批里，
    # 它们的边界本来就对不上——两档一起看才知道一致率是被什么拉下来的。
    for pref, share in (("bd", 0.5), ("bd90", 0.9)):
        bd = [r for r in ledger if r[pref + "_ds"] != ""]
        if not bd:
            print(f"\n边界（单条 MSU 行覆盖 ≥{share:.0%}）：无元件满足，不适用")
            continue
        within = [r for r in bd if abs(r[pref + "_ds"]) <= 10 and abs(r[pref + "_de"]) <= 10]
        mode = Counter((r[pref + "_ds"], r[pref + "_de"]) for r in within).most_common(1)
        print(f"\n边界（单条 MSU 行覆盖 ≥{share:.0%} 的 {len(bd)} 条）：两端都在 ±10 bp 的 "
              f"{len(within)} = {len(within) / len(bd) * 100:.1f}%")
        if mode:
            print(f"  众数 (Δstart, Δend) = {mode[0][0]}（{mode[0][1]} 条）")
        d = sorted(r[pref + "_de"] - r[pref + "_ds"] for r in within)
        if d:
            print(f"  Δend − Δstart 中位 {d[len(d) // 2]} bp（负=MSU 的元件更短）")

    if args.tsv:
        cols = ["chrom", "start", "end", "strand", "length", "ed_class",
                "ed_superfamily", "tier", "identity", "tsd", "verdict_code",
                "msu_dominant", "msu_superfamily", "mixed"]
        cols += [c for c, _ in COV_COLS]
        cols += ["n_msu_hits", "msu_note", "bd_ds", "bd_de", "bd90_ds", "bd90_de",
                 "family"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in ledger:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\n台账已写: {args.tsv}（{len(ledger)} 行）")


if __name__ == "__main__":
    main()
