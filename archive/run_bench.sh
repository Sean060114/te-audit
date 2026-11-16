#!/usr/bin/env bash
# 为 te_check 生成第二、三份 EDTA 输出（多基因组 benchmark 素材）
#
# 设计意图：单一变量
#   col1m   与 Chr1 同一基因组（A. thaliana Col-0），但只有 1 Mb  → 隔离"尺度效应"
#   alyrata 近缘种 A. lyrata，2 Mb                              → 隔离"物种效应"
#   若 col1m 的判据表现与 Chr1 一致，而 alyrata 偏离，则偏离可归因于物种。
#
# 参数与 Chr1 保持一致（--species others，不加 --curatedlib，不加 --sensitive），
# 否则差异无法归因。
#
# 用法（在你的 WSL 终端里）：
#   bash /mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/run_bench.sh
#
# 可重复执行：已完成的任务会自动跳过。

set -uo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta

TESTDIR="$CONDA_PREFIX/share/EDTA/test"
OUT="$HOME/te_bench"
THREADS="${THREADS:-8}"

mkdir -p "$OUT"
cd "$OUT" || exit 1

echo "EDTA : $(command -v EDTA.pl)"
echo "测试集: $TESTDIR"
echo "输出  : $OUT"
echo "线程  : $THREADS"
echo

run_one () {
    local src="$1" tag="$2"
    local done_marker="$tag.fa.mod.EDTA.intact.gff3"

    if [ -s "$done_marker" ]; then
        echo "[跳过] $tag 已存在结果：$done_marker"
        return 0
    fi
    if [ ! -s "$src" ]; then
        echo "[错误] 找不到输入：$src"
        return 1
    fi

    cp -f "$src" "$tag.fa"
    echo "=============================================================="
    echo "[$(date '+%F %T')] 开始 $tag"
    echo "=============================================================="

    EDTA.pl --genome "$tag.fa" \
            --species others \
            --step all \
            --anno 1 \
            --overwrite 1 \
            --threads "$THREADS" 2>&1 | tee "$tag.run.log"

    if [ -s "$done_marker" ]; then
        echo "[$(date '+%F %T')] $tag 完成 ✅"
    else
        echo "[$(date '+%F %T')] $tag 未产出 intact.gff3 ❌  见 $tag.run.log"
    fi
    echo
}

run_one "$TESTDIR/Col.test.fa"     col1m
run_one "$TESTDIR/Alyrata.test.fa" alyrata

echo "=============================================================="
echo "结果清单："
for f in col1m alyrata; do
    if [ -s "$f.fa.mod.EDTA.intact.gff3" ]; then
        printf "  %-10s intact.gff3 %8s 行   raw 分类表 %s 个\n" \
            "$f" "$(grep -vc '^#' "$f.fa.mod.EDTA.intact.gff3")" \
            "$(ls "$f.fa.mod.EDTA.raw"/*/*.rexdb.cls.tsv 2>/dev/null | wc -l)"
    else
        printf "  %-10s 缺失\n" "$f"
    fi
done
echo "ALL_DONE"
