#!/usr/bin/env bash
# 第三版：失败后先修复命名，再续跑
#
# 【归因已撤回，2026-11-16 二次核查】此前写过两版解释，均被源码推翻：
#   版本一：LTR_retriever.py:70 的 glob(genome+'.mod.*') 拼写缺陷 → 错。
#   版本二：EDTA issue #599 的"序列名超过 rmblastn 约 50 字符上限"→ 对本场景不成立，
#           因为 seqid 是 >Chr1 / >OsChr1（4–6 字符），水稻与 Chr1 却都稳定复现。
# 已对着 v2.2.2 与 LTR_retriever v2.9.9 源码逐条排除 EDTA.pl、EDTA_raw.pl、
# LTR_retriever 的 13 字符截断分支（条件均不满足），详见 run_bench_rice.sh:124-160
# 的完整查证记录。剩余可疑项是 LTR_retriever 内置 RepeatMasker 步骤把已 mask 的
# 输入再 mask 一次，**尚未验证**，验证命令同样记在那里。
# 结论：改名修复已被四个基因组实证有效（col1m / alyrata / 水稻 / Chr1），
# 但它治的是症状，成因待定，不要再引用被推翻的两版解释。
#
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
#   bash /mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/run_bench3.sh

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
