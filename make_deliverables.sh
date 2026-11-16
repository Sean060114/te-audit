#!/usr/bin/env bash
# ============================================================
# 从 EDTA 注释产出下游交付文件（v2）
#
# v1 的问题：家族汇总表按"逐条长度相加"统计，嵌套元件会重复计数，
#             合计 10.10%，与软屏蔽实际值 9.20% 对不上，投稿时会被质疑。
# v2 修正：同时给出
#   bp_redun  —— 逐条相加（含嵌套重复，用于"元件总数"口径）
#   bp_union  —— 同类内合并去重（用于"该类占基因组比例"口径）
#   并额外给出去重后的全局合计，应与软屏蔽比例一致。
# 同时把 TEsorter 的三字母代码映射为标准名称。
#
# 运行：bash make_deliverables.sh
# ============================================================
set -euo pipefail

cd ~/edta_run
source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta

GENOME=chr1.fa.mod
GFF=${GENOME}.EDTA.TEanno.gff3
BED=${GENOME}.EDTA.TE.bed
SOFT=${GENOME}.softmasked.fa
FAM=${GENOME}.EDTA.family.summary.tsv
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

TOTAL=$(awk '/^>/{next} {n+=length($0)} END{print n}' "$GENOME")
echo "基因组总长: $TOTAL bp"

# 三字母代码 → 标准名称（依据 TEanno.sum 的条目数逐类核对得到）
map_name() {
    case "$1" in
        DNA/DTA)        echo "TIR/hAT" ;;
        DNA/DTC)        echo "TIR/CACTA" ;;
        DNA/DTH)        echo "TIR/PIF_Harbinger" ;;
        DNA/DTM)        echo "TIR/Mutator" ;;
        DNA/DTT)        echo "TIR/Tc1_Mariner" ;;
        MITE/DTA)       echo "MITE/hAT" ;;
        MITE/DTH)       echo "MITE/PIF_Harbinger" ;;
        MITE/DTM)       echo "MITE/Mutator" ;;
        DNA/Helitron)   echo "nonTIR/helitron" ;;
        nonDNA/helitron) echo "nonTIR/helitron" ;;
        *)              echo "$1" ;;
    esac
}

# ---------- 1. BED（排除附属特征，避免同一段重复计数）----------
awk -F'\t' -v OFS='\t' '
    $0 ~ /^#/ {next}
    $3=="target_site_duplication" || $3=="long_terminal_repeat" || $3=="primer_binding_site" {next}
    {
        cls=$3
        if (match($9, /classification=[^;]+/)) cls=substr($9, RSTART+15, RLENGTH-15)
        print $1, $4-1, $5, cls, ".", ($7=="-"?"-":"+")
    }' "$GFF" | sort -k1,1 -k2,2n > "$BED"
echo "BED: $(wc -l < "$BED") 条 → $BED"

# ---------- 2. 软屏蔽基因组 ----------
if command -v bedtools >/dev/null 2>&1; then
    bedtools maskfasta -soft -fi "$GENOME" -bed "$BED" -fo "$SOFT"
    MASKED=$(awk '/^>/{next} {n+=length($0); l+=gsub(/[acgt]/,"")} END{printf "%d", l}' "$SOFT")
    echo "软屏蔽基因组 → $SOFT  ($MASKED bp)"
else
    echo "WARNING: bedtools 不在 PATH，跳过软屏蔽"; MASKED=0
fi

# ---------- 3. 家族级汇总 ----------
# 3a. 每条元件一行：cls \t chr \t start0 \t end
awk -F'\t' -v OFS='\t' '
    $0 ~ /^#/ {next}
    $3=="target_site_duplication" || $3=="long_terminal_repeat" || $3=="primer_binding_site" {next}
    {
        cls=$3
        if (match($9, /classification=[^;]+/)) cls=substr($9, RSTART+15, RLENGTH-15)
        print cls, $1, $4-1, $5
    }' "$GFF" > "$TMP/rows.tsv"

# 3b. 条目数与逐条长度
awk -F'\t' -v OFS='\t' '{n[$1]++; r[$1]+=$4-$3} END{for(c in n) print c, n[c], r[c]}' \
    "$TMP/rows.tsv" | sort > "$TMP/count.tsv"

# 3c. 同类内合并去重
sort -k1,1 -k2,2 -k3,3n "$TMP/rows.tsv" | awk -F'\t' -v OFS='\t' '
    {
        if ($1!=pc || $2!=pchr) { if (pc!="") print pc, pu; pc=$1; pchr=$2; pe=$4; pu=$4-$3 }
        else if ($3>pe) { pu += $4-$3; pe=$4 }      # 有间隔：只加本段长度
        else if ($4>pe) { pu += $4-pe; pe=$4 }      # 部分重叠：只加超出部分
    }
    END{ if (pc!="") print pc, pu }' | sort > "$TMP/union.tsv"

# 3d. 全局去重合计（不分类别，应与软屏蔽一致）
sort -k2,2 -k3,3n "$TMP/rows.tsv" | awk -F'\t' '
    { if ($2!=pchr) { if (pchr!="") t+=en-strt; pchr=$2; strt=$3; en=$4 }
      else if ($3>en) { t+=en-strt; strt=$3; en=$4 }
      else if ($4>en) en=$4 }
    END{ if (pchr!="") t+=en-strt; print t }' > "$TMP/grand"

# 3e. 拼表
{
    printf "class_std\tclass_raw\tentries\tbp_redun\tpct_redun\tbp_union\tpct_union\n"
    join -t$'\t' -1 1 -2 1 "$TMP/count.tsv" "$TMP/union.tsv" \
    | while IFS=$'\t' read -r raw n r u; do
        printf "%s\t%s\t%s\t%s\t%.4f\t%s\t%.4f\n" \
            "$(map_name "$raw")" "$raw" "$n" "$r" "$(awk -v a="$r" -v b="$TOTAL" 'BEGIN{print a/b*100}')" \
            "$u" "$(awk -v a="$u" -v b="$TOTAL" 'BEGIN{print a/b*100}')"
      done | sort -t$'\t' -k6,6nr
} > "$FAM"

GRAND=$(cat "$TMP/grand")
echo "家族汇总 → $FAM"
echo ""
column -t -s$'\t' "$FAM" 2>/dev/null || cat "$FAM"
echo ""
echo "-----------------------------------------------"
awk -v g="$GRAND" -v t="$TOTAL" -v m="$MASKED" \
    'BEGIN{printf "全局去重覆盖: %d bp = %.2f%%  （软屏蔽实测 %d bp = %.2f%%）\n", g, g/t*100, m, m/t*100}'
echo "说明：bp_redun 含嵌套重复，合计会大于全局去重值；引用时用 bp_union / 全局去重口径。"
