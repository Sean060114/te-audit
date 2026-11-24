#!/usr/bin/env bash
# 分步/重试运行 EDTA，规避 LTR_retriever 的双层 .mod 命名问题
#
# 已确认的事实（不是猜测）：
#   1. LTR_retriever.py:70-71 用 glob(genome+'.mod.*') 判断是否已屏蔽，
#      命中就把基因组名再追加一次 '.mod'。
#   2. EDTA_raw 期望的是单层：$genome.LTR.raw.fa / $genome.pass.list。
#      一旦变成双层（X.fa.mod.mod.*），下游全部找不到 → LTR 结果 0 bp → 整个运行失败。
#   3. EDTA.pl v2.2.2 只有 ALL / FILTER / FINAL / ANNO 四个标签，
#      没有 LTR 标签，所以 `--step LTR` 是无效参数（会 goto 失败）。
#   4. Chr1 那次成功，靠的是 `--step all --overwrite 0` 连跑两次：
#      第二次时中间文件已就位，命名收敛到单层。
#
# 因此策略：对每个基因组反复跑 `--step all --overwrite 0`，
# 直到出现单层的 X.fa.mod.pass.list 为止（最多 3 次），再跑 anno。
#
# 用法：
#   bash archive/run_bench2.sh
#
# 输出落在 ~/te_bench2（全新目录，不复用上次失败的中间文件）。

set -uo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta

TESTDIR="$CONDA_PREFIX/share/EDTA/test"
OUT="$HOME/te_bench2"
THREADS="${THREADS:-8}"
MAX_TRY="${MAX_TRY:-3}"

mkdir -p "$OUT"
cd "$OUT" || exit 1

echo "EDTA    : $(command -v EDTA.pl)"
echo "测试集  : $TESTDIR"
echo "输出    : $OUT"
echo "每例上限: $MAX_TRY 次尝试"
echo

run_genome () {
    local src="$1" tag="$2"
    local g="$tag.fa"
    local ltrdir="$g.mod.EDTA.raw/LTR"
    local good="$ltrdir/$g.mod.pass.list"
    local bad="$ltrdir/$g.mod.mod.pass.list"

    if [ -s "$g.mod.EDTA.intact.gff3" ]; then
        echo "[跳过] $tag 已完成"; return 0
    fi
    if [ ! -s "$src" ]; then
        echo "[错误] 找不到输入 $src"; return 1
    fi
    cp -f "$src" "$g"

    local ok=0
    for try in $(seq 1 "$MAX_TRY"); do
        echo "--------------------------------------------------------------"
        echo "[$(date '+%T')] $tag 尝试 $try/$MAX_TRY：--step all"
        echo "--------------------------------------------------------------"
        EDTA.pl --genome "$g" --species others --step all \
                --overwrite 0 --threads "$THREADS" 2>&1 | tail -6

        if [ -s "$good" ]; then
            echo "[检查] ✅ 单层命名已就位：$g.mod.pass.list（$(wc -l < "$good") 条）"
            ok=1; break
        elif [ -s "$bad" ]; then
            echo "[检查] 仍是双层 X.fa.mod.mod.*，再试一次"
        else
            echo "[检查] 未见预期的 pass.list，实际存在的相关文件："
            ls "$ltrdir" 2>/dev/null | grep -E "pass\.list|defalse|LTRlib\.fa$" | sed 's/^/         /'
        fi
    done

    if [ "$ok" -ne 1 ] && [ ! -s "$g.mod.EDTA.intact.gff3" ]; then
        echo "[结果] $tag 未能收敛到单层命名，跳过注释步骤"
        echo
        return 1
    fi

    echo
    echo "[$(date '+%T')] $tag：全基因组注释"
    EDTA.pl --genome "$g" --species others --step anno --anno 1 \
            --threads "$THREADS" 2>&1 | tail -6

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
