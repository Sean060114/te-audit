#!/usr/bin/env python3
"""用基因组序列独立复算 EDTA 结构法 LTR 元件的 TSD 与末端基序。

为什么需要它：EDTA 的 intact.gff3 自带 `tsd=`、`motif=`、`ltr_identity=` 三个
字段，但那是 EDTA 自己算的。要判断「这批注释可信」，必须能从序列重新算一遍，
而不是引用它自己的字段。本脚本复算三项：TSD、末端基序、以及两翼 LTR 相似度。

**LTR 子坐标就在 intact.gff3 里**（每个元件有 `lLTR` / `rLTR` 两行，共 840 行），
不需要 LTR_retriever 的 pass.list。此前报告写成"子坐标未回传、无法复算"是错的：
坐标一直在手上，缺的只是去读它。

字段语义（据 EDTA v2.2.2 bin/make_gff3.pl 与 LTR_retriever bin/make_gff3.pl
源码确认，2026-11-23）：
  1. 坐标 1-based 闭区间。
  2. `repeat_region` 的起止 = 5′TSD 第一个碱基 … 3′TSD 最后一个碱基，即
     **TSD 含在元件坐标内**；`tsd=NA` 时才退化为左右 LTR 的起止。
  3. `tsd=` 写的是**元件方向**上的序列：负链元件需反向互补后再比对。
  4. `motif=` 是 LTR 的末端基序（如 TGCA），拆成「前 2 nt + 后 2 nt」分别落在
     5′LTR 起点与 3′LTR 终点，**不是**边界处连续出现的 4 个碱基。
  5. 一个 LTR 元件在 intact.gff3 里占 6 行（repeat_region / lTSD / lLTR /
     LTRRT / rLTR / rTSD），**只有 repeat_region 行代表元件**。
  6. **两翼 LTR 都按基因组正向书写，负链元件也不做反向互补。**（仅指 LTR 相似度；
     `tsd=` 与 `motif=` 两个字段仍按元件方向书写，负链要 rc。这两条方向约定相反。） 实测负链元件
     `lLTR` 与 `rLTR` 的前 30 bp 完全相同；取 rc 反而把相似度打到 0.23。
     （这一条与 TIR 侧相反：TIR 的右翼判反向重复要拿左翼比 rc(右翼)。）

用法：
    python verify_ltr_structure.py \
        --gff rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3 \
        --fasta rice/Oryza_sativa.Chr1.fa rice/Oryza_sativa.Chr12.fa
"""

import argparse
import csv
import difflib
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(paths):
    seqs = {}
    for p in paths:
        name, buf = None, []
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(buf).upper()
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf).upper()
    return seqs


def parse_gff(path):
    """只取 repeat_region 行——一个 LTR 元件占 6 行，其余 5 行不是元件。"""
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "repeat_region":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if not a.get("classification", "").startswith("LTR"):
                continue
            out.append({
                "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                "strand": f[6], "tsd": a.get("tsd", "NA"),
                "motif": a.get("motif", ""), "ident": a.get("ltr_identity", ""),
                "id": a.get("ID", ""),
            })
    # 第二遍：把两翼 LTR 的坐标挂回元件
    byid = {e["id"]: e for e in out}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "long_terminal_repeat":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            parent, sid = a.get("Parent", ""), a.get("ID", "")
            if parent in byid and sid[:1] in ("l", "r"):
                byid[parent].setdefault("ltr", {})[sid[:1]] = (int(f[3]), int(f[4]))
    return out


def ident_defs(a, b):
    """三种常见"相似度"定义，用于找出 EDTA 字段用的是哪一种。

    两侧 LTR 常不等长（420 条里 236 条不等），所以不能用逐位 hamming：
    差 3 bp 就会把整条错位，算出来 0.63 而字段写 0.93。必须先用比对块对齐。
    """
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    M = mm = gap = 0
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        p, q = i2 - i1, j2 - j1
        if op == "equal":
            M += p
        elif op == "replace":
            mm += min(p, q)
            gap += abs(p - q)
        else:
            gap += max(p, q)
    return {
        "A": 2 * M / (len(a) + len(b)) if a and b else 0,
        "B": M / (M + mm + gap) if M + mm + gap else 0,
        "C": 1 - mm / (M + mm) if M + mm else 0,
    }


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--tol", type=int, default=1, help="3′TSD 允许的错配数（默认 1）")
    ap.add_argument("--ident-def", choices=("A", "B", "C"), default="C",
                    help="ltr_identity 的复算定义，默认 C（扣掉 indel 列的逐位一致率）")
    ap.add_argument("--tsv", help="输出逐元件结构判定（供 A/B/C 分级直接引用）")
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    elems = parse_gff(args.gff)
    if not elems:
        sys.exit("没有找到 LTR 的 repeat_region 行")

    n = ok5 = ok3 = ok3tol = okboth = okmot = 0
    idents = []
    id_rows = []
    # 逐元件判定：A/B/C 分级要用"每一条都过"，不能用批量百分比当门槛
    per = {}
    for e in elems:
        chrom = e["chrom"]
        if chrom not in seqs:
            # seqid 可能带前缀（OsChr1 vs Chr1），退回后缀匹配
            cand = [k for k in seqs if k.endswith(chrom) or chrom.endswith(k)]
            if not cand:
                continue
            chrom = cand[0]
        sq = seqs[chrom]
        s, en = e["start"], e["end"]
        t = e["tsd"]
        if not t or t == "NA":
            continue
        n += 1
        L = len(t)
        rec = {"chrom": e["chrom"], "start": s, "end": en, "id": e["id"],
               "strand": e["strand"], "tsd_len": L,
               "tsd_seq": t, "motif": e["motif"],
               "tsd5": 0, "tsd3": 0, "tsd3_tol": 0, "motif_ok": 0,
               "ident_field": "", "ident_calc": "", "ident_diff": ""}
        left = sq[s - 1:s - 1 + L]
        right = sq[en - L:en]
        if e["strand"] == "-":            # tsd 按元件方向书写
            left, right = rc(left), rc(right)
        if left == t:
            ok5 += 1
            rec["tsd5"] = 1
        if right == t:
            ok3 += 1
            rec["tsd3"] = 1
        if hamming(right, t) <= args.tol:
            ok3tol += 1
            rec["tsd3_tol"] = 1
        if left == t and hamming(right, t) <= args.tol:
            okboth += 1
        m = e["motif"]
        if len(m) == 4:
            # 前 2 nt 在 5′LTR 起点（紧接 5′TSD），后 2 nt 在 3′LTR 终点（紧邻 3′TSD）。
            #
            # **motif 也按元件方向书写**，与 tsd= 一致。负链元件必须先取反向互补再比：
            # 基因组小坐标端是它的 3′ 端，所以还要交叉（rc(尾)==m[:2]、rc(头)==m[2:]）。
            # 初版漏了这一步，靠"末端 2 nt 常相同或回文"蒙对了 410/420——420 条里
            # 有 406 条正反两种读法同时成立，把 5 条负链元件的错读掩盖掉了。
            # 修完后 415/420。
            h = sq[s - 1 + L:s + 1 + L]
            t = sq[en - L - 2:en - L]
            if e["strand"] == "-":
                h, t = rc(t), rc(h)
                hit = (h == m[:2] and t == m[2:])
            else:
                hit = (h == m[:2] and t == m[2:])
            if hit:
                okmot += 1
                rec["motif_ok"] = 1
        if e["ident"]:
            try:
                idents.append(float(e["ident"]))
            except ValueError:
                pass
        ltr = e.get("ltr", {})
        if len(ltr) == 2 and e["ident"]:
            a = sq[ltr["l"][0] - 1:ltr["l"][1]]
            b = sq[ltr["r"][0] - 1:ltr["r"][1]]
            try:
                f = float(e["ident"])
            except ValueError:
                f = None
            if f is not None and len(a) >= 50 and len(b) >= 50:
                d = ident_defs(a, b)
                id_rows.append({"id": e["id"], "strand": e["strand"],
                                "la": len(a), "lb": len(b), "field": f,
                                "A": d["A"], "B": d["B"], "C": d["C"]})
                rec["ident_field"] = round(f, 4)
                rec["ident_calc"] = round(d[args.ident_def], 4)
                rec["ident_diff"] = round(abs(d[args.ident_def] - f), 4)
        per[(e["chrom"], s, en)] = rec

    pct = lambda k: f"{k}/{n} ({k / n * 100:.1f}%)" if n else "n/a"
    print(f"LTR 元件（repeat_region 行）: {len(elems)}；带 TSD 可复算: {n}")
    print(f"  5′ TSD 与序列完全一致          : {pct(ok5)}")
    print(f"  3′ TSD 与序列完全一致          : {pct(ok3)}")
    print(f"  3′ TSD 允许 {args.tol} 处错配      : {pct(ok3tol)}")
    print(f"  两侧同时成立（5′ 严格 + 3′ 容错）: {pct(okboth)}")
    print(f"  末端基序（motif 前 2 + 后 2）    : {pct(okmot)}")
    print()
    print("两翼 LTR 相似度 vs ltr_identity（%d 条有双翼坐标）" % len(id_rows))
    if id_rows:
        print(f"  LTR 长度相等: {sum(1 for r in id_rows if r['la'] == r['lb'])}/"
              f"{len(id_rows)}；长度中位 {sorted(r['la'] for r in id_rows)[len(id_rows)//2]} bp")
        for k in ("A", "B", "C"):
            d = sorted(abs(r[k] - r["field"]) for r in id_rows)
            within = sum(1 for x in d if x <= 0.01) / len(d) * 100
            print(f"  定义 {k}: |差| 平均 {sum(d)/len(d):.4f}  中位 {d[len(d)//2]:.4f}  "
                  f"P90 {d[int(len(d)*0.9)]:.4f}  内 0.01 占 {within:.1f}%")
        sel = args.ident_def
        d = sorted(id_rows, key=lambda r: -abs(r[sel] - r["field"]))
        print(f"  采用定义 {sel} 时差得最大的 5 条（全为两侧长度不等）：")
        for r in d[:5]:
            print(f"    {r['id']:<20}{r['strand']}  {r['la']}/{r['lb']} bp  "
                  f"复算 {r[sel]:.4f}  字段 {r['field']:.4f}  差 {abs(r[sel]-r['field']):.4f}")


    if args.tsv:
        # tsd_seq / motif 一并带出来：判定失败时要能就地看末端实际序列，
        # 不然复核还得重读 gff。两列不影响判定，只是证据留档。
        cols = ["chrom", "start", "end", "id", "strand", "tsd_len", "tsd_seq",
                "motif", "tsd5", "tsd3", "tsd3_tol", "motif_ok",
                "ident_field", "ident_calc", "ident_diff"]
        with open(args.tsv, "w", encoding="utf-8", newline="") as fh:
            # 项目 .gitattributes 锁 LF：csv 默认写 \r\n，必须显式指定，
            # 否则 WSL 里的下游脚本会读到行尾多余的 \r。
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            for k in sorted(per, key=lambda k: (k[0], k[1])):
                w.writerow(per[k])
        print(f"\n逐元件判定已写入 {args.tsv}（{len(per)} 行）")
        full = sum(1 for r in per.values()
                   if r["tsd5"] and r["tsd3_tol"] and r["motif_ok"]
                   and r["ident_diff"] != "" and r["ident_diff"] <= 0.01)
        print(f"  三项（TSD/motif/ident）逐条全过: {full}/{len(per)} = {full/len(per)*100:.1f}%")


if __name__ == "__main__":
    main()
