#!/usr/bin/env bash
# ############################################################
# 【已废弃 —— 不要运行】2026-11 实测结论：
#   Dfam-RepeatMasker.lib 这个子集里【一条植物序列都没有】
#   （26,292 条几乎全是昆虫/脊椎动物/线虫；曾有 "Oryza" 命中，
#    实为米象 Sitophilus oryzae，是昆虫）
#   → 用它验证拟南芥 TE 注释没有意义，结果不可用
#
# 已改用：Araport11 官方分类 TE 注释（同样基于 TAIR10 坐标）
#   报告：EDTA_vs_Araport11_评估报告.md
#   脚本：compare_araport11.py（在 Windows 侧 python 直接跑）
# ############################################################
# ============================================================
# EDTA 注释交叉验证 v2：用 Dfam 共识库（15.8 MB）跑 RepeatMasker
#
# v1 失败原因：conda 装的 RepeatMasker 4.1.5 里
#   - Dfam.h5 只有 11 字节，内容是字面量 "placeholder"
#   - RepeatMaskerLib.h5 是坏的占位（I/O error）
#   - RepeatMasker.lib 里植物类群条目为 0（全是昆虫/鱼/哺乳动物）
# 所以 RepeatMasker 找不到任何家族，直接退出。
#
# 本版改用 Dfam 官方提供的 RepeatMasker 共识库（FASTA + -lib 方式），
# 兼容 4.1.5，无需升级 RepeatMasker，下载量仅 15.8 MB。
# ============================================================
set -euo pipefail

cd ~/edta_run
source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta

GENOME=chr1.fa.mod
LOG=validation2.log
LIBDIR=dfam_lib

echo "=== [1/4] 检查 / 下载 Dfam 共识库（15.8 MB）===" | tee "$LOG"
mkdir -p "$LIBDIR"
if [ ! -s "$LIBDIR/Dfam-RepeatMasker.lib" ]; then
    echo "库未就绪，调用 get_dfam_lib.sh（带超时与多源回退）..." | tee -a "$LOG"
    bash /mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/get_dfam_lib.sh 2>&1 | tee -a "$LOG"
fi
N=$(grep -c '^>' "$LIBDIR/Dfam-RepeatMasker.lib" 2>/dev/null || true)
N=${N:-0}
echo "库中家族数: $N" | tee -a "$LOG"
if [ "$N" -lt 1000 ]; then
    echo "ERROR: 库不完整（家族数 $N）。请先单独运行 get_dfam_lib.sh 并确认成功。" | tee -a "$LOG"
    exit 1
fi

echo "" | tee -a "$LOG"
echo "=== [2/4] 检查植物类群覆盖 ===" | tee -a "$LOG"
for k in Viridiplantae Embryophyta Arabidopsis Oryza Brassicaceae; do
    printf "  %-16s %s\n" "$k" "$(grep -c "@$k" "$LIBDIR/Dfam-RepeatMasker.lib" || true)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "=== [3/4] RepeatMasker（Dfam 库，约 5-15 分钟）===" | tee -a "$LOG"
rm -rf RM_dfam
mkdir -p RM_dfam
RepeatMasker -pa 8 -qq -no_is -xsmall -gff \
    -lib "$LIBDIR/Dfam-RepeatMasker.lib" \
    -dir RM_dfam "$GENOME" 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== [4/4] 结果汇总 ===" | tee -a "$LOG"
TBL="RM_dfam/${GENOME}.tbl"
if [ -f "$TBL" ]; then
    echo "--- RepeatMasker（Dfam）类别统计 ---" | tee -a "$LOG"
    sed -n '/^  Class/,/^=====/p' "$TBL" | head -40 | tee -a "$LOG"
    echo "" | tee -a "$LOG"
    echo "--- 总计 ---" | tee -a "$LOG"
    grep -E "total inter|Total" "$TBL" | tee -a "$LOG"
else
    echo "WARNING: 未找到 $TBL" | tee -a "$LOG"
fi

echo "" | tee -a "$LOG"
echo "--- 软屏蔽比例 ---" | tee -a "$LOG"
if [ -f "RM_dfam/${GENOME}.masked" ]; then
    awk '/^>/{next} {n+=length($0); l+=gsub(/[acgt]/,"")} END{printf "  Dfam 软屏蔽: %.2f%%\n", l/n*100}' \
        "RM_dfam/${GENOME}.masked" | tee -a "$LOG"
fi
echo "  EDTA 软屏蔽: $(awk '/^>/{next} {n+=length($0); l+=gsub(/[acgt]/,"")} END{printf "%.2f%%", l/n*100}' "$GENOME.softmasked.fa" 2>/dev/null || echo '未生成')" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "完成。请把 validation2.log 和 RM_dfam/${GENOME}.tbl 发我，我做类别级差异分析。"
