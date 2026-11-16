#!/usr/bin/env python3
"""
模块 1 第一张真实报表：EDTA 结构注释（intact.gff3）的质量事实基线。

回答三个问题：
  1. 不同类别的元件，各自"有多少证据"？（identity / TSD / TIR / motif / strand_evidence）
  2. 结构证据本身是否站得住？（TIR 是否真能构成反向重复，允许 ±10 bp 滑移）
  3. 哪些类别在注释里存在、但几乎没有序列层面的支撑？

输出：
  - 控制台报告
  - intact_quality_baseline.tsv（逐元件指标）
  - intact_quality_summary.tsv（逐类别汇总）
"""
import os
import statistics
from collections import Counter, defaultdict

GFF = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod.EDTA.intact.gff3"
OUT_TSV = "intact_quality_baseline.tsv"
OUT_SUM = "intact_quality_summary.tsv"

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
ELEMENT_TYPES = {
    "CACTA_TIR_transposon", "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon",
    "LTR_retrotransposon", "Mutator_TIR_transposon", "PIF_Harbinger_TIR_transposon",
    "Tc1_Mariner_TIR_transposon", "hAT_TIR_transposon", "helitron",
}
SUB_FEATURES = {"long_terminal_repeat", "target_site_duplication", "repeat_region"}


def rc(s):
    return s.translate(COMP)[::-1]


def attrs(a):
    d = {}
    for kv in a.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def ir_identity(left, right, max_shift=10):
    """右臂滑动对齐 revcomp(左臂)，返回最佳匹配率与最优位移。"""
    target = rc(left)
    if not target or not right:
        return None, None
    best, best_shift = -1.0, None
    n = len(target)
    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            a, b = target[shift:], right[:n - shift] if shift < n else ""
        else:
            a, b = target, right[-shift:]
            a = a[:len(b)]
        if not b:
            continue
        m = min(len(a), len(b))
        if m < 10:
            continue
        ident = sum(1 for i in range(m) if a[i] == b[i]) / m
        if ident > best:
            best, best_shift = ident, shift
    return (best if best >= 0 else None), best_shift


def parse():
    recs = []
    with open(GFF, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype not in ELEMENT_TYPES:
                continue
            d = attrs(f[8])
            start, end = int(f[3]), int(f[4])
            rec = dict(
                ftype=ftype,
                start=start, end=end, len=end - start + 1,
                strand=f[6],
                cls=d.get("classification", ""),
                identity=None, identity_src="",
                tsd=None, tsd_pct=None,
                tir_left=None, tir_right=None,
                ir_ident=None, ir_shift=None,
                ltr_ident=None, motif=None, strand_ev=None,
            )
            # identity 来源
            if d.get("identity") and d["identity"] != "NA":
                try:
                    rec["identity"] = float(d["identity"])
                    rec["identity_src"] = "identity"
                except ValueError:
                    pass
            if d.get("ltr_identity"):
                try:
                    rec["ltr_ident"] = float(d["ltr_identity"])
                except ValueError:
                    pass
            # TSD
            tsd = d.get("TSD") or d.get("tsd")
            if tsd:
                parts = tsd.split("_")
                rec["tsd"] = parts[0] if parts else tsd
                if len(parts) >= 3:
                    try:
                        rec["tsd_pct"] = float(parts[2])
                    except ValueError:
                        pass
            # TIR
            tir = d.get("TIR")
            if tir and "_" in tir:
                p = tir.split("_")
                if len(p) == 2:
                    rec["tir_left"], rec["tir_right"] = p
                    ident, shift = ir_identity(p[0], p[1])
                    rec["ir_ident"], rec["ir_shift"] = ident, shift
            rec["motif"] = d.get("motif")
            rec["strand_ev"] = d.get("strand_evidence")
            recs.append(rec)
    return recs


def main():
    recs = parse()
    print(f"元素级记录: {len(recs)} 条（已排除 LTR/TSD/repeat_region 子特征）")

    # ---- 逐元件输出 ----
    cols = ["ftype", "start", "end", "len", "strand", "cls",
            "identity", "ltr_ident", "tsd", "tsd_pct",
            "tir_left", "tir_right", "ir_ident", "ir_shift", "motif", "strand_ev"]
    with open(OUT_TSV, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"[写出] {OUT_TSV}")

    # ---- 逐类别汇总 ----
    by = defaultdict(list)
    for r in recs:
        by[r["ftype"]].append(r)

    print("\n" + "=" * 104)
    print(f"{'类别':<30}{'n':>4}{'长度中位':>10}{'有identity':>11}{'有TSD':>8}{'有TIR':>7}{'完美IR':>8}")
    print("=" * 104)
    rows = []
    for t in sorted(by, key=lambda x: -len(by[x])):
        g = by[t]
        n = len(g)
        med = statistics.median(r["len"] for r in g)
        has_id = sum(1 for r in g if r["identity"] is not None or r["ltr_ident"] is not None)
        has_tsd = sum(1 for r in g if r["tsd"])
        has_tir = sum(1 for r in g if r["tir_left"])
        perfect = sum(1 for r in g if r["ir_ident"] is not None and r["ir_ident"] >= 0.999)
        print(f"{t:<30}{n:>4}{med:>10.0f}{has_id:>11}{has_tsd:>8}{has_tir:>7}{perfect:>8}")
        rows.append((t, n, med, has_id, has_tsd, has_tir, perfect))

    with open(OUT_SUM, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("type\tn\tlen_median\thas_identity\thas_TSD\thas_TIR\tperfect_IR\n")
        for r in rows:
            fh.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.0f}" for x in r) + "\n")
    print(f"[写出] {OUT_SUM}")

    # ---- 证据缺口 ----
    print("\n=== 证据缺口（按类别）===")
    total = len(recs)
    no_id = [r for r in recs if r["identity"] is None and r["ltr_ident"] is None]
    print(f"  无任何一致性数值的元件: {len(no_id)} / {total} = {len(no_id) / total * 100:.1f}%")
    for t, c in Counter(r["ftype"] for r in no_id).most_common():
        print(f"    {t}: {c}")

    # ---- TIR 反向重复质量 ----
    tir = [r for r in recs if r["ir_ident"] is not None]
    if tir:
        print(f"\n=== TIR 反向重复质量（{len(tir)} 个元件，允许 ±10 bp 滑移）===")
        buckets = Counter()
        for r in tir:
            i = r["ir_ident"]
            buckets["完美 (=1.000)"] += i >= 0.999
            buckets["近乎完美 (0.95–0.999)"] += 0.95 <= i < 0.999
            buckets["较强 (0.85–0.95)"] += 0.85 <= i < 0.95
            buckets["中等 (0.70–0.85)"] += 0.70 <= i < 0.85
            buckets["弱 (<0.70)"] += i < 0.70
        for k, v in buckets.items():
            if v:
                print(f"  {k:<22} {v:>3}  ({v / len(tir) * 100:.1f}%)")
        ident = [r["ir_ident"] for r in tir]
        print(f"  中位数 {statistics.median(ident):.3f}  最小 {min(ident):.3f}")
        shifts = Counter(abs(r["ir_shift"]) for r in tir if r["ir_shift"] is not None)
        print(f"  需要的滑移量分布: {dict(sorted(shifts.items()))}")

    # ---- LTR 证据类型 ----
    ltr = [r for r in recs if r["ftype"].endswith("LTR_retrotransposon")]
    if ltr:
        print(f"\n=== LTR 元件的链方向证据构成（{len(ltr)} 个）===")
        for ev, c in Counter(r["strand_ev"] or "(无)" for r in ltr).most_common():
            print(f"  {ev:<32} {c}")
        idents = [r["ltr_ident"] for r in ltr if r["ltr_ident"] is not None]
        if idents:
            print(f"  LTR 一致性: 中位 {statistics.median(idents):.4f}，"
                  f"范围 {min(idents):.4f}–{max(idents):.4f}")
        no_motif = sum(1 for r in ltr if not r["motif"])
        print(f"  缺 motif 的: {no_motif} / {len(ltr)}")

    # ---- 长度 ----
    print("\n=== 长度分布（各前 3 短/长）===")
    for t in sorted(by, key=lambda x: -len(by[x]))[:4]:
        g = sorted(by[t], key=lambda r: r["len"])
        print(f"  {t:<30} 最短 {g[0]['len']:>6}  中位 {statistics.median(r['len'] for r in g):>6.0f}  最长 {g[-1]['len']:>6}")


if __name__ == "__main__":
    main()
