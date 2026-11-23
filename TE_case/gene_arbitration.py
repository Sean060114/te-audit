#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 MSU 的**基因注释流程**裁决重复注释流程给出的标签。

为什么这是独立证据
------------------
MSU RGAP7 有两套互不通气的产物：
  repeat_masked.gff3   重复流程：RepeatMasker + 自建文库，给的是"这个区间是什么重复"
  all_models.gff3      基因流程：证据驱动的基因模型，Note 里写功能

后者会给出 `retrotransposon protein`（11,399 个）、`transposon protein`（3,822 个）
这类**转座子自身编码的基因**。它们与重复注释来自不同流程、不同证据，可以用来
裁决"这个位点到底是 LTR 还是 DNA 转座子"。

这正是 LTR 台账里 23 条 `msu_nonltr_label` 需要的东西：EDTA 说 LTR/Gypsy，MSU 的
重复流程说 DNA/TIR（AnacB），两边都只有一套证据。基因流程给出第三票。

先做一次标签可靠性审计（这是本脚本最重要的一段）
------------------------------------------------
`retrotransposon protein` 与 `transposon protein` 这两个 Note **不能直接当裁决票**。
实测：长度 6,156 bp 的基因模型在两种标签下各出现 11 次——同一条模型、同一个家族，
标签取决于它落在哪个位点。而这个家族恰好主导了 LTR 台账里 23 条冲突。

所以脚本先按"模型长度"做审计，把同时出现两种标签的长度列出来；判读结果时必须先看
这一节。标签一致的部分可以用，冲突的部分不能用。

判据
----
对每个元件，看内部有没有**完整**（两端都在元件内）的转座子基因模型：
  retro    Note 含 retrotransposon / reverse transcriptase / integrase / gag / rnase
  dnate    Note 含 transposon（"retrotransposon" 优先判，不会误吞）
  host     其余

注意：**不能用"元件内含完整基因"当指标**。LTR 元件长 6–9 kb，而水稻基因密度是
每 6.7 kb 一个，62% 的随机同长窗口都含完整基因——指标饱和，没有区分力。所以这里
只统计转座子类的基因模型，宿主基因只作对照列。

用法
----
  python TE_case/gene_arbitration.py --ledger TE_case/ltr_msu_ledger.tsv \
    --models C:/Users/Sean/AppData/Local/Temp/msu/all_models.gff3.gz --detail
"""
import argparse
import gzip
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from urllib.parse import unquote

RETRO = re.compile(r"retrotransposon|reverse transcriptase|integrase|\bgag\b|rnase", re.I)
DNATE = re.compile(r"transposon|transposase", re.I)


def classify(note):
    if RETRO.search(note):
        return "retro"
    if DNATE.search(note):
        return "dnate"
    return "host"


def parse_genes(path, chroms):
    """chrom -> (starts, [(start, end, id, note, kind)])，0-based 半开。"""
    rows = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[0] not in chroms:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if e <= s:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            note = unquote(a.get("Note", ""))
            rows[f[0]].append((s, e, a.get("ID", a.get("Name", "")), note, classify(note)))
    out = {}
    for c, v in rows.items():
        v.sort()
        out[c] = ([x[0] for x in v], v)
    return out


def audit_labels(g):
    """按模型长度审计 retro/dnate 标签是否自相矛盾。"""
    by = defaultdict(Counter)
    for _c, v in g.items():
        for a, b, _gid, _note, kind in v[1]:
            if kind != "host":
                by[b - a][kind] += 1
    bad = {L: c for L, c in by.items() if c["retro"] >= 3 and c["dnate"] >= 3}
    tot = sum(c["retro"] + c["dnate"] for c in by.values())
    nbad = sum(c["retro"] + c["dnate"] for c in bad.values())
    print(f"\n标签审计：{len(by)} 个不同长度的转座子基因模型，共 {tot} 个")
    if bad:
        print(f"  同一长度下两种标签各 ≥3 次的有 {len(bad)} 个长度，涉及 {nbad} 个模型"
              f"（{nbad / tot * 100:.1f}%）——这些长度的标签**不可用于裁决**")
        for L, c in sorted(bad.items(), key=lambda x: -(x[1]["retro"] + x[1]["dnate"]))[:6]:
            print(f"    {L:>6} bp  retro {c['retro']:>3}  dnate {c['dnate']:>3}")
    else:
        print("  没有发现同长度双标签的情况")
    return bad


def hits_in(g, c, s, e):
    """返回压在元件内的基因模型列表。"""
    if c not in g:
        return []
    starts, rows = g[c]
    i = bisect_left(starts, s - 300000)      # 最长基因模型 ~200 kb
    out = []
    while i < len(rows) and rows[i][0] < e:
        a, b, gid, note, kind = rows[i]
        i += 1
        ov = min(e, b) - max(s, a)
        if ov > 100:
            out.append({"id": gid, "note": note, "kind": kind, "start": a, "end": b,
                        "bp": b - a, "ov": ov,
                        "complete": a >= s and b <= e,
                        "cs": max(a, s), "ce": min(b, e),
                        "cov": ov / (e - s)})
    return out


def union_cov(hits, L):
    """同一类别多个模型压在同一元件上时取并集占比。直接求和会把重叠算两遍。"""
    if not hits or L <= 0:
        return 0.0
    iv = sorted((x["cs"], x["ce"]) for x in hits)
    tot, cs, ce = 0, iv[0][0], iv[0][1]
    for a, b in iv[1:]:
        if a <= ce:
            ce = max(ce, b)
        else:
            tot += ce - cs
            cs, ce = a, b
    return min(1.0, (tot + ce - cs) / L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--detail", action="store_true", help="逐条打印有转座子基因模型的元件")
    args = ap.parse_args()

    import csv
    led = list(csv.DictReader(open(args.ledger, encoding="utf-8"), delimiter="\t"))
    chroms = {(r["chrom"][2:] if r["chrom"].startswith("Os") else r["chrom"]) for r in led}
    print(f"台账 {len(led)} 条；染色体 {', '.join(sorted(chroms))}")
    g = parse_genes(args.models, chroms)
    n = Counter(x[4] for v in g.values() for x in v[1])
    print(f"基因模型：retro {n['retro']}，dnate {n['dnate']}，host {n['host']}")
    bad = audit_labels(g)

    rows = []
    for r in led:
        c = r["chrom"][2:] if r["chrom"].startswith("Os") else r["chrom"]
        s, e = int(r["start"]) - 1, int(r["end"])
        h = hits_in(g, c, s, e)
        retro = [x for x in h if x["kind"] == "retro"]
        dnate = [x for x in h if x["kind"] == "dnate"]
        rows.append({
            "chrom": r["chrom"], "start": r["start"], "end": r["end"],
            "length": int(r["length"]), "verdict_code": r["verdict_code"],
            "ed_class": r["ed_class"],
            "n_retro": len(retro), "n_dnate": len(dnate),
            "retro_complete": sum(1 for x in retro if x["complete"]),
            "dnate_complete": sum(1 for x in dnate if x["complete"]),
            # 覆盖率用并集：多个模型可能相互重叠，直接求和会超过 1
            "retro_cov": round(union_cov(retro, e - s), 3),
            "dnate_cov": round(union_cov(dnate, e - s), 3),
            "retro_ids": ";".join(x["id"] for x in retro[:3]),
            "dnate_ids": ";".join(x["id"] for x in dnate[:3]),
        })

    print(f"\n{'判定':<26}{'n':>4}{'有 retro':>10}{'有 dnate':>10}"
          f"{'retro 覆盖中位':>15}{'dnate 覆盖中位':>15}")
    by = defaultdict(list)
    for r in rows:
        by[r["verdict_code"]].append(r)
    for v, rs in sorted(by.items(), key=lambda x: -len(x[1])):
        nr = sum(1 for x in rs if x["n_retro"])
        nd = sum(1 for x in rs if x["n_dnate"])
        cr = sorted(x["retro_cov"] for x in rs)[len(rs) // 2]
        cd = sorted(x["dnate_cov"] for x in rs)[len(rs) // 2]
        print(f"{v:<26}{len(rs):>4}{nr:>10}{nd:>10}{cr:>15.3f}{cd:>15.3f}")

    # ---- 裁决：两种转座子基因同时出现时谁覆盖得多
    both = [r for r in rows if r["n_retro"] and r["n_dnate"]]
    if both:
        print(f"\n两种转座子基因模型都有的 {len(both)} 条：")
        for r in both:
            side = "retro" if r["retro_cov"] >= r["dnate_cov"] else "dnate"
            print(f"  {r['chrom']}:{r['start']:<10} {r['verdict_code']:<24} "
                  f"retro {r['retro_cov']:.2f} vs dnate {r['dnate_cov']:.2f} → {side}")

    if args.detail:
        print("\n逐条（只列有转座子基因模型的）：")
        for r in rows:
            if not (r["n_retro"] or r["n_dnate"]):
                continue
            print(f"  {r['chrom']}:{r['start']}-{r['end']} ({r['length']} bp) "
                  f"{r['verdict_code']}  {r['ed_class']}")
            if r["n_retro"]:
                print(f"      retro {r['n_retro']}（完整 {r['retro_complete']}）"
                      f" 覆盖 {r['retro_cov']:.2f}  {r['retro_ids']}")
            if r["n_dnate"]:
                print(f"      dnate {r['n_dnate']}（完整 {r['dnate_complete']}）"
                      f" 覆盖 {r['dnate_cov']:.2f}  {r['dnate_ids']}")


if __name__ == "__main__":
    main()
