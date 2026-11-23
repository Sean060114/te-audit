#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把逐元件结构复算与 MSU 对账合起来，重新划 A/B/C。

为什么要重划
------------
`基准_LTR逐元件台账.md` 第七节把 A 集定义成「结构自证 (99.0%) + 跨源超家族一致」，
185 条。但那个 **99.0% 是个批量数字，不是逐元件门槛**：99.0% 出自 5′TSD 一项，
其余两项（末端基序 97.6%、两翼相似度 95.0%）根本没参与筛选。换句话说，A 集
185 条里可能混着结构不成立的元件——这个漏洞直到 `verify_ltr_structure.py` 能
逐条输出才补得上。

本脚本做两件事：

1. 把「结构自证」落实成**每一条都要过**的三项门槛，而不是引用一个比例。
2. 与 `link_msu.py` 的逐元件判定交叉，重新给出 A/B/C，并说明每一条为什么落在
   那一级。

结构自证的三项门槛（与字段校准文档一致）
----------------------------------------
- 5′TSD 与序列**严格**一致（tsd5）
- 3′TSD 允许 1 处错配（tsd3_tol）
- 末端基序两端都在位（motif_ok）
- 两翼 LTR 复算相似度 **≥ 0.90**（ident_ok）

**关于最后一项，不能用 |复算−字段| ≤ 0.01。** 那是「字段可复现性」的判据，不是
「结构成立」的判据。实测 420 条的复算值全部 ≥ 0.9317，一条都不低于 0.90；而
|差| > 0.01 的 20 条**全部是复算值高于字段值**（方向一致，没有一条反向），是
「indel 是否计入分母」的定义偏差，不是结构缺陷。混用这两个判据会把 20 条结构
完好的元件错降一级。

三级怎么划
----------
A：结构三项全过 + 跨源**超家族**一致（agree_superfamily）
B：结构三项全过 + 不冲突，但超家族无双侧旁证（msu_ltr_only / ed_unknown_msu_resolved
   / both_unresolved）
C：结构不全过，或本身有争议/无第二来源（msu_nonltr_label / conflict_superfamily /
   insufficient_cov）

用法
----
  python TE_case/grade_ltr_abc.py \
    --structure TE_case/ltr_structure.tsv \
    --ledger TE_case/ltr_msu_ledger.tsv \
    --tsv TE_case/ltr_abc.tsv
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict

# 判定码 → 所属档位。新增判定码时必须在这里显式归类，不允许静默落进 B。
TIER_OF = {
    "agree_superfamily": "A",
    "msu_ltr_only": "B",
    "ed_unknown_msu_resolved": "B",
    "both_unresolved": "B",
    "msu_nonltr_label": "C",
    "conflict_superfamily": "C",
    "insufficient_cov": "C",
}

IDENT_MIN = 0.90


def load(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def structure_ok(r):
    """返回 (是否全过, 失败项列表)。"""
    bad = []
    if r["tsd5"] != "1":
        bad.append("tsd5")
    if r["tsd3_tol"] != "1":
        bad.append("tsd3")
    if r["motif_ok"] != "1":
        bad.append("motif")
    if not r["ident_calc"] or float(r["ident_calc"]) < IDENT_MIN:
        bad.append("ident")
    return (not bad), bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", required=True, help="verify_ltr_structure.py --tsv 的产物")
    ap.add_argument("--ledger", required=True, help="link_msu.py --kind ltr --tsv 的产物")
    ap.add_argument("--tsv", help="输出逐元件分级表")
    ap.add_argument("--ident-min", type=float, default=IDENT_MIN)
    args = ap.parse_args()

    st = load(args.structure)
    led = load(args.ledger)
    key = lambda r: (r["chrom"], int(r["start"]), int(r["end"]))

    by_st = {key(r): r for r in st}
    by_led = {key(r): r for r in led}
    only_st = set(by_st) - set(by_led)
    only_led = set(by_led) - set(by_st)
    if only_st or only_led:
        sys.exit(f"两份表对不上：结构表多 {len(only_st)} 条、台账多 {len(only_led)} 条。"
                 f"两表必须来自同一次运行的同一批元件。")

    rows, cross = [], defaultdict(Counter)
    for k in by_led:
        s, l = by_st[k], by_led[k]
        sok, bad = structure_ok(s)
        base = TIER_OF.get(l["verdict_code"])
        if base is None:
            sys.exit(f"未知判定码 {l['verdict_code']}，请先在 TIER_OF 里归类")
        # A 要求结构全过；B 也要求结构全过（否则它是 C，不是"结构完好的旁证不足"）
        tier = base if sok or base == "C" else "C"
        cross[l["verdict_code"]][tier] += 1
        rows.append({
            "chrom": l["chrom"], "start": l["start"], "end": l["end"],
            "length": l["length"], "strand": l["strand"],
            "ed_class": l["ed_class"], "verdict_code": l["verdict_code"],
            "structure_ok": int(sok), "structure_fail": "+".join(bad),
            "tsd_len": s["tsd_len"], "ident_field": s["ident_field"],
            "ident_calc": s["ident_calc"], "abc": tier,
        })

    n = len(rows)
    print(f"元件 {n} 条（结构表 {len(st)}、台账 {len(led)}，键完全对齐）")
    print(f"结构三项全过（ident 门槛 {args.ident_min:.2f}）: "
          f"{sum(r['structure_ok'] for r in rows)}/{n} = "
          f"{sum(r['structure_ok'] for r in rows) / n * 100:.1f}%")
    fail = Counter(r["structure_fail"] for r in rows if not r["structure_ok"])
    if fail:
        print("  未全过的失败项组合：")
        for k, v in fail.most_common():
            print(f"    {v:>3} 条  {k}")

    print("\n结构自证 × MSU 判定 交叉表（行=判定码，列=最终档位）")
    print(f"  {'判定码':<24}{'→A':>5}{'→B':>5}{'→C':>5}   合计")
    for code in sorted(cross, key=lambda c: -sum(cross[c].values())):
        c = cross[code]
        tot = sum(c.values())
        print(f"  {code:<24}{c['A']:>5}{c['B']:>5}{c['C']:>5}   {tot:>3}")

    tier_cnt = Counter(r["abc"] for r in rows)
    print(f"\n最终 A/B/C = {tier_cnt['A']} / {tier_cnt['B']} / {tier_cnt['C']}")
    old = {"A": 185, "B": 181, "C": 54}
    print(f"  原口径（结构用批次数 99.0% 代替逐条门槛）= "
          f"{old['A']} / {old['B']} / {old['C']}")
    for t in ("A", "B", "C"):
        d = tier_cnt[t] - old[t]
        if d:
            print(f"  {t} 变化 {d:+d}")

    # 原来记 A 但结构不全过的：就是那条漏洞的具体代价
    dg = [r for r in rows if r["verdict_code"] == "agree_superfamily"
          and not r["structure_ok"]]
    if dg:
        print(f"\n原记 A 但结构不全过的 {len(dg)} 条（应降级）：")
        for r in dg:
            print(f"  {r['chrom'][2:]}:{r['start']:<9}{r['ed_class']:<11}"
                  f"失败 {r['structure_fail']:<14}"
                  f"ident {r['ident_calc'] or '-'}")

    up = [r for r in rows if r["verdict_code"] in
          ("msu_ltr_only", "ed_unknown_msu_resolved") and r["structure_ok"]]
    print(f"\n结构全过但只到 B 的 {len(up)} 条——它们的短板是 MSU 的词汇粒度，"
          f"不是结构（第 11 节已论证，可考虑升 A-）")

    if args.tsv:
        cols = ["chrom", "start", "end", "length", "strand", "ed_class",
                "verdict_code", "structure_ok", "structure_fail", "tsd_len",
                "ident_field", "ident_calc", "abc"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            for r in sorted(rows, key=lambda r: (r["chrom"], int(r["start"]))):
                w.writerow(r)
        print(f"\n逐元件分级已写入 {args.tsv}")


if __name__ == "__main__":
    main()
