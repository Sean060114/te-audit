#!/usr/bin/env python3
"""
从 EDTA 注释中筛出「Araport11 未注释」的候选新元件，并按可信度分级。

输入：
  araport11_TE.gff3                      Araport11 分类 TE 注释（TAIR10 坐标）
  chr1.fa.mod.EDTA.TEanno.gff3           EDTA 注释
输出：
  EDTA_candidate_novel.tsv               候选清单（含分级与理由）
  stdout                                 分级汇总

分级规则（基于本次实测的类级精确率）：
  A 级 —— LTR（repeat_region/method=structural，有 motif+tsd）：该类精确率 95%
  B 级 —— 有结构证据（TIR 或 TSD 或 motif）的非 helitron 元件
  C 级 —— ≥1 kb 且 method=structural 的非 helitron 元件
  D 级 —— 其余零重叠元件（helitron 为主，精确率 40.8%，需人工核验）
"""
import sys
import numpy as np

AP_GFF = "araport11_TE.gff3"
ED_GFF = ("//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/"
          "chr1.fa.mod.EDTA.TEanno.gff3")
GENOME_LEN = 30427671
OUT_TSV = "EDTA_candidate_novel.tsv"

SKIP_FEATURES = {"target_site_duplication", "long_terminal_repeat",
                 "primer_binding_site"}

# TEsorter 三字母码 -> 标准名
CODE2NAME = {
    "DTA": "hAT", "DTC": "CACTA", "DTH": "PIF_Harbinger", "DTM": "Mutator",
    "DTT": "Tc1_Mariner", "DTE": "En-Spm", "DTF": "PIF",
    "RIL": "L1", "RIT": "tRNA", "RIX": "Unknown",
    "RLG": "Gypsy", "RLC": "Copia", "RLX": "Unknown",
    "RC": "Helitron", "Helitron": "Helitron",
}


def attr_get(attr, key):
    """从 GFF3 第 9 列取属性值，未找到返回空串。"""
    for kv in attr.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return ""


def norm_class(raw):
    """把 EDTA 的 classification 归一到「超家族/家族」粒度。"""
    if not raw:
        return "Unknown"
    parts = raw.split("/")
    if len(parts) >= 2:
        top, sub = parts[0], parts[1]
        if top in ("DNA", "MITE"):
            return "DNA/" + CODE2NAME.get(sub, sub)
        if top == "LTR":
            return "LTR/" + CODE2NAME.get(sub, sub)
        if top == "LINE":
            return "LINE/" + CODE2NAME.get(sub, sub)
        if top == "SINE":
            return "SINE/" + CODE2NAME.get(sub, sub)
        if top == "RC":
            return "RC/Helitron"
        return raw
    return CODE2NAME.get(raw, raw)


def build_ap_mask():
    mask = np.zeros(GENOME_LEN, dtype=bool)
    for line in open(AP_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != "Chr1":
            continue
        if f[2] != "transposable_element":
            continue
        mask[int(f[3]) - 1:int(f[4])] = True
    return mask


def read_edta():
    rows = []
    for line in open(ED_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] in SKIP_FEATURES:
            continue
        s, e = int(f[3]) - 1, int(f[4])
        a = f[8]
        cls_raw = attr_get(a, "classification")
        rows.append({
            "chrom": f[0], "start": s, "end": e, "len": e - s,
            "feat": f[2],
            "name": attr_get(a, "Name"),
            "cls": norm_class(cls_raw or f[2]),
            "cls_raw": cls_raw or f[2],
            "method": attr_get(a, "method"),
            "ident": attr_get(a, "identity") or attr_get(a, "ltr_identity"),
            "tir": attr_get(a, "TIR"),
            "tsd": attr_get(a, "TSD") or attr_get(a, "tsd"),
            "motif": attr_get(a, "motif"),
            "pbs": attr_get(a, "strand_evidence"),
        })
    return rows


def grade(r):
    """返回 (等级, 理由)。等级越小越可信。"""
    has_struct = bool(r["tir"] or r["tsd"] or r["motif"])
    is_ltr = r["cls"].startswith("LTR/")
    is_rp = r["feat"] == "repeat_region"
    is_rc = r["cls"].startswith("RC/")

    if is_ltr and is_rp:
        return "A", "intact LTR，结构完整（该类精确率 95%）"
    if has_struct and not is_rc:
        return "B", "有 TIR/TSD/motif 结构证据"
    if r["len"] >= 1000 and r["method"] == "structural" and not is_rc:
        return "C", "≥1kb 且 method=structural"
    if is_rc:
        return "D", "helitron（该类精确率仅 40.8%，需人工核验）"
    return "D", "仅同源支持，无结构证据"


def main():
    ap = build_ap_mask()
    rows = read_edta()
    zero = [r for r in rows if not ap[r["start"]:r["end"]].any()]

    for r in zero:
        r["grade"], r["reason"] = grade(r)

    zero.sort(key=lambda r: (r["grade"], -r["len"]))

    with open(OUT_TSV, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(["grade", "chrom", "start", "end", "length",
                            "class", "class_raw", "method", "identity",
                            "has_TIR", "has_TSD", "motif", "feat", "name",
                            "reason"]) + "\n")
        for r in zero:
            fh.write("\t".join([
                r["grade"], r["chrom"], str(r["start"]), str(r["end"]),
                str(r["len"]), r["cls"], r["cls_raw"], r["method"],
                r["ident"], "Y" if r["tir"] else "", "Y" if r["tsd"] else "",
                r["motif"], r["feat"], r["name"], r["reason"],
            ]) + "\n")

    print(f"EDTA 元件 {len(rows)} 个，与 Araport11 零重叠 {len(zero)} 个\n")

    print("=== 分级汇总 ===")
    print(f"{'等级':<6}{'个数':>6}{'bp':>12}   含义")
    print("-" * 78)
    reasons = {}
    for g in "ABCD":
        sub = [r for r in zero if r["grade"] == g]
        if not sub:
            continue
        reasons[g] = sub[0]["reason"]
        print(f"{g:<6}{len(sub):>6}{sum(r['len'] for r in sub):>12,}   {sub[0]['reason']}")

    print("\n=== A/B/C 级按类别拆分 ===")
    import collections
    for g in "ABC":
        sub = [r for r in zero if r["grade"] == g]
        if not sub:
            continue
        c = collections.Counter(r["cls"] for r in sub)
        print(f"  {g} 级 {len(sub)} 个: " +
              ", ".join(f"{k}×{v}" for k, v in c.most_common()))

    print(f"\n清单已写出: {OUT_TSV}")
    print("\n提示：A 级可直接作为「EDTA 新发现的完整 LTR」上报；"
          "D 级（helitron 为主）必须人工核验或直接丢弃。")


if __name__ == "__main__":
    sys.exit(main())
