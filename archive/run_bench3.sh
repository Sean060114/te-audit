#!/usr/bin/env bash
# 第三版：失败后先修复命名，再续跑
#
# 【归因已查明，此前两版均被推翻，2026-11-16】
#   版本一：LTR_retriever.py:70 的 glob(genome+'.mod.*') 拼写缺陷 → 无源码支持。
#   版本二：EDTA issue #599「序列名超过 rmblastn 约 50 字符上限」→ 与 4–6 字符的
#           seqid（>Chr1 / >OsChr1）矛盾，套不上。
# 真实成因：EDTA.pl 已把基因组归一化成 X.fa.mod（必要的一步），但 v2.2.2 的
#   EDTA_raw.pl:440 调用 LTR_retriever 时没传 -convert_seq_name，后者默认 1，
#   于是又归一化了第二次——encode_fasta 的输出路径是 "$genome.mod"，
#   结果就是 X.fa.mod.mod。与序列名长短无关（短名也会写出该文件并被采纳）。
#   v2.3.0:410 起调用里补了 -convert_seq_name 0，故不再出现。
# 完整链条见 ../run_bench_rice.sh 的 repair_names 注释。改名修复仍有效，
# 且不改动环境，是当前唯一安全做法（升级会破坏五份数据的可比性）。
# 已确证的事实（对着 ~/te_bench2 的真实文件核对）：
#   LTR_retriever 的产物是完整的，只是全在双前缀下：
#     col1m.fa.mod.mod.pass.list        627 B
#     col1m.fa.mod.mod.pass.list.gff3   6.8 K
#     col1m.fa.mod.mod.defalse          3.4 K
#     col1m.fa.mod.mod.LTRlib.fa         38 K
#   而单前缀下的 LTR.raw.fa / LTRlib.fa 都是 0 字节。
#
#   失败链：EDTA_raw 读 $genome.pass.list（单前缀）→ 文件不存在
#   → awk 报 "cannot open" → 生成的 LTR.raw.fa 为空 → "Raw LTR results not found"
#
# 因此修复方式是把 X.mod.mod.* 复制成 X.mod.*（纯改名，不改内容）。
# 这些文件的内容只含序列 ID 与坐标，不含文件名，改名是安全的。
#
# 用法：
#   bash archive/run_bench3.sh

set -uo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta

TESTDIR="$CONDA_PREFIX/share/EDTA/test"
OUT="$HOME/te_bench2"          # 复用上一轮已算好的中间结果，不重算
THREADS="${THREADS:-8}"

cd "$OUT" || { echo "目录不存在：$OUT"; exit 1; }

repair_names () {
    local g="$1"
    local ltrdir="$g.mod.EDTA.raw/LTR"
    local n=0
    for f in "$ltrdir/$g.mod.mod."*; do
        [ -e "$f" ] || continue
        local base new
        base=$(basename "$f")
        new=$(printf '%s' "$base" | sed 's/\.mod\.mod\./.mod./')
        [ "$base" = "$new" ] && continue
        cp -f "$f" "$ltrdir/$new" && n=$((n+1))
    done
    echo "  修复：复制 $n 个文件（X.mod.mod.* → X.mod.*）"
    # 验证关键文件
    for k in pass.list defalse LTRlib.fa; do
        if [ -s "$ltrdir/$g.mod.$k" ]; then
            printf "   %-12s %8s 字节 ✅\n" "$k" "$(stat -c%s "$ltrdir/$g.mod.$k")"
        else
            printf "   %-12s 缺失或为空 ❌\n" "$k"
        fi
    done
}

run_genome () {
    local src="$1" tag="$2"
    local g="$tag.fa"

    if [ -s "$g.mod.EDTA.intact.gff3" ]; then
        echo "[跳过] $tag 已完成"; return 0
    fi
    if [ ! -s "$src" ]; then
        echo "[错误] 找不到输入 $src"; return 1
    fi
    [ -s "$g" ] || cp -f "$src" "$g"

    echo "=============================================================="
    echo "[$(date '+%T')] $tag：修复命名后继续"
    echo "=============================================================="
    repair_names "$g"

    echo
    echo "[$(date '+%T')] $tag：重跑 --step all（复用已有中间结果）"
    EDTA.pl --genome "$g" --species others --step all \
            --overwrite 0 --threads "$THREADS" 2>&1 | tail -8

    if [ ! -s "$g.mod.EDTA.raw/LTR/$g.mod.LTR.raw.fa" ]; then
        echo "[结果] $tag 的 LTR 仍未生成，停止。LTR 目录现状："
        ls -lh "$g.mod.EDTA.raw/LTR" 2>/dev/null \
            | awk 'NR>1{printf "         %-58s %s\n",$9,$5}' | tail -20
        return 1
    fi
    echo "[检查] LTR.raw.fa 已非空（$(stat -c%s "$g.mod.EDTA.raw/LTR/$g.mod.LTR.raw.fa") 字节）"

    echo
    echo "[$(date '+%T')] $tag：全基因组注释"
    EDTA.pl --genome "$g" --species others --step anno --anno 1 \
            --threads "$THREADS" 2>&1 | tail -8

    if [ -s "$g.mod.EDTA.intact.gff3" ]; then
        echo "[$(date '+%T')] $tag 完成 ✅ intact.gff3 $(grep -vc '^#' "$g.mod.EDTA.intact.gff3") 行"
    else
        echo "[$(date '+%T')] $tag 未产出 intact.gff3 ❌"
    fi
    echo
}

run_genome "$TESTDIR/Col.test.fa"     col1m
run_genome "$TESTDIR/Alyrata.test.fa" alyrata

echo "=============================================================="
for f in col1m alyrata; do
    if [ -s "$f.fa.mod.EDTA.intact.gff3" ]; then
        printf "  %-8s ✅  %6s 行  raw 分类表 %s 个\n" "$f" \
            "$(grep -vc '^#' "$f.fa.mod.EDTA.intact.gff3")" \
            "$(ls "$f.fa.mod.EDTA.raw"/*/*.rexdb.cls.tsv 2>/dev/null | wc -l)"
    else
        printf "  %-8s ❌ 缺失\n" "$f"
    fi
done
echo "ALL_DONE"
