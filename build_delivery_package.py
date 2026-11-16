#!/usr/bin/env python3
"""
构建可交付的 TE 注释包（TAIR10 Chr1）。

命名统一为 Chr1：
  Araport11 原始注释用 Chr1，EDTA 与 chr1.fa.mod 用 1。
  交付包以 Chr1 为准（拟南芥通行命名，且主注释无需改名，减少出错面）。

产出（delivery/）：
  TAIR10.Chr1.genome.softmasked.fa     软屏蔽基因组（重复区小写，19.79%）
  TAIR10.Chr1.TE.Araport11.gff3        主注释（7,135 个 TE + 7,890 个子片段）
  TAIR10.Chr1.TE.EDTA.gff3             EDTA 原始注释（仅改名，供溯源）
  TAIR10.Chr1.TE.supplement.gff3       43 条高分化 TIR 元件（EDTA 独有）
  TAIR10.Chr1.TE.final.gff3            交付用合并注释（主注释 + 补充）
  TAIR10.Chr1.TE.mask.bed              并集重复区间（BED3，供 maskfasta）
  TAIR10.Chr1.TE.family_summary.tsv    家族级拷贝数/bp/占比
"""
import os
import sys
from collections import defaultdict
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WORK, "delivery")

GENOME = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod"
AP_GFF = os.path.join(WORK, "araport11_TE.gff3")
ED_GFF = ("//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/"
          "chr1.fa.mod.EDTA.TEanno.gff3")
SUPP_TSV = os.path.join(WORK, "EDTA_B_level_candidates.tsv")

SEQID = "Chr1"
GENOME_LEN = 30427671
WRAP = 60
SKIP_ED = {"target_site_duplication", "long_terminal_repeat", "primer_binding_site"}

# EDTA feature type -> 标准分类（写入 TE_class= 便于下游识别）
ED_CLASS = {
    "helitron": "RC/Helitron",
    "Gypsy_LTR_retrotransposon": "LTR/Gypsy",
    "Copia_LTR_retrotransposon": "LTR/Copia",
    "LTR_retrotransposon": "LTR/unknown",
    "Mutator_TIR_transposon": "DNA/Mutator",
    "PIF_Harbinger_TIR_transposon": "DNA/PIF-Harbinger",
    "CACTA_TIR_transposon": "DNA/CACTA",
    "hAT_TIR_transposon": "DNA/hAT",
    "Tc1_Mariner_TIR_transposon": "DNA/Tc1-Mariner",
    "L1_LINE_retrotransposon": "LINE/L1",
    "tRNA_SINE_retrotransposon": "SINE/tRNA",
}


def read_fasta(path):
    name, buf = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip()
            else:
                buf.append(line.strip())
    return name, "".join(buf)


def attrs(a):
    d = {}
    for kv in a.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def merge_iv(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def union_bp(iv):
    return sum(e - s for s, e in merge_iv(iv))


def main():
    os.makedirs(OUT, exist_ok=True)
    _, seq = read_fasta(GENOME)
    assert len(seq) == GENOME_LEN, f"基因组长度异常: {len(seq)}"
    print(f"基因组 {SEQID}: {len(seq):,} bp")

    # ---------- 1. 读取 Araport11（Chr1） ----------
    ap_te, ap_frag = [], []
    for line in open(AP_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != "Chr1":
            continue
        rec = (f[2], int(f[3]), int(f[4]), f[6], f[8])
        if f[2] == "transposable_element":
            ap_te.append(rec)
        elif f[2] == "transposon_fragment":
            ap_frag.append(rec)
    print(f"Araport11: TE {len(ap_te)} 条，子片段 {len(ap_frag)} 条")

    # ---------- 2. 读取 EDTA（仅改名 + 归一化分类） ----------
    ed_recs = []
    ed_strand = {}
    for line in open(ED_GFF, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        d = attrs(f[8])
        cls = d.get("classification") or ED_CLASS.get(f[2], "")
        extra = f";TE_class={cls}" if cls and "TE_class" not in f[8] else ""
        ed_recs.append((f[2], int(f[3]), int(f[4]), f[6], f[8] + extra))
        ed_strand[(int(f[3]), int(f[4]))] = f[6]
    print(f"EDTA: 共 {len(ed_recs)} 行")

    # ---------- 3. 读取补充候选 ----------
    supp = []
    with open(SUPP_TSV, encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(head)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(head):
                continue
            supp.append(p)
    print(f"补充候选: {len(supp)} 条")

    # ---------- 4. 写 Araport11 主注释 ----------
    p_ap = os.path.join(OUT, "TAIR10.Chr1.TE.Araport11.gff3")
    with open(p_ap, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        fh.write("##assembly: TAIR10.1\n")
        fh.write("##primary TE annotation: Araport11 (TAIR), "
                 "classified transposable_element\n")
        fh.write(f"##seqid: renamed Chr1 (unchanged), records: "
                 f"{len(ap_te)} TE + {len(ap_frag)} transposon_fragment\n")
        for t, s, e, st, a in ap_te + ap_frag:
            fh.write(f"{SEQID}\tAraport11\t{t}\t{s}\t{e}\t.\t{st}\t.\t{a}\n")
    print(f"  写出 {os.path.basename(p_ap)}")

    # ---------- 5. 写 EDTA 注释（改名） ----------
    p_ed = os.path.join(OUT, "TAIR10.Chr1.TE.EDTA.gff3")
    with open(p_ed, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        fh.write("##assembly: TAIR10.1\n")
        fh.write("##de novo TE annotation: EDTA v2.2.2 "
                 "(--step anno --anno 1 --species others)\n")
        fh.write("##seqid renamed 1 -> Chr1 (coordinates unchanged)\n")
        fh.write("##added attribute TE_class= normalized from feature type\n")
        for t, s, e, st, a in ed_recs:
            fh.write(f"{SEQID}\tEDTA\t{t}\t{s}\t{e}\t.\t{st}\t.\t{a}\n")
    print(f"  写出 {os.path.basename(p_ed)}")

    # ---------- 6. 写补充注释 ----------
    p_sp = os.path.join(OUT, "TAIR10.Chr1.TE.supplement.gff3")
    sp_recs = []
    with open(p_sp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        fh.write("##assembly: TAIR10.1\n")
        fh.write("##supplementary TE annotation: EDTA de novo calls with no "
                 "overlap with Araport11\n")
        fh.write("##all records are DIVERGENT TIR elements (identity 0.81-0.92), "
                 "NOT intact elements\n")
        fh.write("##qc: TIR_mismatch = mismatches between left TIR revcomp and "
                 "right TIR; perfect_IR_after_shift = 1 only if <=1 mismatch "
                 "after shifting boundaries\n")
        for i, p in enumerate(supp, 1):
            gid = f"EDTAsup_{i:03d}"
            cls = p[idx["class"]]
            s, e = int(p[idx["start"]]), int(p[idx["end"]])
            strand = ed_strand.get((s, e), ".")
            a = (f"ID={gid};Name={gid};classification={cls};"
                 f"evidence=EDTA_de_novo;identity={p[idx['identity']]};"
                 f"TIR_len={p[idx['TIR_len']]};"
                 f"TIR_mismatch={p[idx['TIR_mismatch']]};"
                 f"perfect_IR_after_shift={p[idx['perfect_IR_after_shift']]};"
                 f"shift_bp={p[idx['shift_bp']]};"
                 f"left_TIR={p[idx['left_TIR']]};"
                 f"right_TIR={p[idx['right_TIR']]};"
                 f"left_TSD={p[idx['left_TSD']]}")
            fh.write(f"{SEQID}\tEDTA\ttransposable_element\t{s}\t{e}\t.\t"
                     f"{strand}\t.\t{a}\n")
            sp_recs.append((cls, s, e))
    print(f"  写出 {os.path.basename(p_sp)}  ({len(sp_recs)} 条)")

    # ---------- 7. 合并注释 ----------
    p_fin = os.path.join(OUT, "TAIR10.Chr1.TE.final.gff3")
    with open(p_fin, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##gff-version 3\n")
        fh.write("##assembly: TAIR10.1\n")
        fh.write("##TE annotation for delivery: Araport11 (primary) + "
                 "EDTA divergent TIR supplement\n")
        fh.write(f"##records: {len(ap_te)} Araport11 TE, {len(ap_frag)} "
                 f"Araport11 fragments, {len(sp_recs)} EDTA supplement\n")
        fh.write("##note: supplement records have zero overlap with Araport11\n")
        for t, s, e, st, a in ap_te + ap_frag:
            fh.write(f"{SEQID}\tAraport11\t{t}\t{s}\t{e}\t.\t{st}\t.\t{a}\n")
        for line in open(p_sp, encoding="utf-8"):
            if line.startswith("#"):
                continue
            fh.write(line)
    print(f"  写出 {os.path.basename(p_fin)}")

    # ---------- 8. 掩码 BED ----------
    iv = [(s - 1, e) for _, s, e, _, _ in ap_te]
    iv += [(s - 1, e) for _, s, e, _, _ in ed_recs if _ not in SKIP_ED]
    iv += [(s - 1, e) for _, s, e in sp_recs]
    merged = merge_iv(iv)
    p_bed = os.path.join(OUT, "TAIR10.Chr1.TE.mask.bed")
    with open(p_bed, "w", encoding="utf-8", newline="\n") as fh:
        for s, e in merged:
            fh.write(f"{SEQID}\t{s}\t{e}\n")
    mask_bp = sum(e - s for s, e in merged)
    print(f"  写出 {os.path.basename(p_bed)}  {len(merged):,} 区间, "
          f"{mask_bp:,} bp ({mask_bp/GENOME_LEN*100:.2f}%)")

    # ---------- 9. 软屏蔽基因组 ----------
    mask = np.zeros(GENOME_LEN, dtype=bool)
    for s, e in merged:
        mask[s:e] = True
    lut = str.maketrans("ACGTNacgtn", "acgtnacgtn")
    chars = list(seq.upper())
    for i in np.flatnonzero(mask):
        chars[i] = chars[i].translate(lut)
    out_seq = "".join(chars)
    p_fa = os.path.join(OUT, "TAIR10.Chr1.genome.softmasked.fa")
    with open(p_fa, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f">{SEQID} dna:chromosome chromosome:TAIR10:1:1:"
                 f"{GENOME_LEN}:1 REF softmasked=Araport11_union_EDTA_"
                 f"{mask.sum()/GENOME_LEN*100:.2f}pct\n")
        for i in range(0, GENOME_LEN, WRAP):
            fh.write(out_seq[i:i + WRAP] + "\n")
    low = sum(1 for c in out_seq if c.islower())
    print(f"  写出 {os.path.basename(p_fa)}  小写 {low:,} bp "
          f"({'✓ 与掩码一致' if low == mask.sum() else '✗ 不一致'})")

    # ---------- 10. 家族汇总 ----------
    fam = defaultdict(list)
    for t, s, e, st, a in ap_te:
        d = attrs(a)
        fam[d.get("classification", "unclassified")].append((s - 1, e))
    rows = []
    for cls, ivs in fam.items():
        bp_s = sum(e - s for s, e in ivs)
        bp_u = union_bp(ivs)
        parts = cls.split("/")
        rows.append((cls,
                     parts[0] if parts else "",
                     parts[1] if len(parts) > 1 else "",
                     len(ivs), bp_s, bp_u,
                     bp_u / GENOME_LEN * 100))
    rows.sort(key=lambda r: -r[5])
    p_tab = os.path.join(OUT, "TAIR10.Chr1.TE.family_summary.tsv")
    with open(p_tab, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("classification\tsuperfamily\tfamily\tcopies\tbp_sum\t"
                 "bp_union\tpct_genome\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t"
                     f"{r[6]:.4f}\n")
        tot_s = sum(r[4] for r in rows)
        tot_u = union_bp([(s, e) for ivs in fam.values() for s, e in ivs])
        fh.write(f"TOTAL\t.\t.\t{sum(r[3] for r in rows)}\t{tot_s}\t{tot_u}\t"
                 f"{tot_u/GENOME_LEN*100:.4f}\n")
    print(f"  写出 {os.path.basename(p_tab)}  {len(rows)} 个家族, "
          f"去重 {tot_u:,} bp ({tot_u/GENOME_LEN*100:.2f}%)")

    # ---------- 11. 汇总 ----------
    # sp_recs 存的是 1-based inclusive（GFF 口径），转 0-based 半开再算 bp
    sp_bp = union_bp([(s - 1, e) for _, s, e in sp_recs])
    ap_bp = union_bp([(s - 1, e) for _, s, e, _, _ in ap_te])
    print("\n=== 交付汇总 ===")
    print(f"  Araport11 主注释 : {len(ap_te):>6} TE  {ap_bp:>10,} bp "
          f"({ap_bp/GENOME_LEN*100:.2f}%)")
    print(f"  EDTA 补充        : {len(sp_recs):>6}      {sp_bp:>10,} bp "
          f"({sp_bp/GENOME_LEN*100:.2f}%)")
    print(f"  合计（零重叠可加）: {len(ap_te)+len(sp_recs):>6}      "
          f"{ap_bp+sp_bp:>10,} bp ({(ap_bp+sp_bp)/GENOME_LEN*100:.2f}%)")
    print(f"  屏蔽并集         : {len(merged):>6} 区间 {mask_bp:>10,} bp "
          f"({mask_bp/GENOME_LEN*100:.2f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
