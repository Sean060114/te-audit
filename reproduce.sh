#!/usr/bin/env bash
# ==============================================================================
# reproduce.sh —— 从仓库内的 EDTA 输出重建全部分析产物，并与已入库版本逐字节比对
#
# 用法（在仓库根目录）：
#   bash reproduce.sh                  # 台账 + 六数据集表，约 1.5 分钟
#   bash reproduce.sh --with-reports   # 另加四份 te_check JSON，需基因组 FASTA
#
# 退出码：全部一致 → 0；任一不一致 → 1。产物写进 _repro/（已被 .gitignore 忽略），
# 不覆盖任何已入库文件，所以跑失败也不会弄脏工作区。
#
# 为什么值得跑：这些 tsv/json 是报告里每个数字的依据。判据改一行，这里应当失败；
# 若改了判据却仍"全部一致"，说明比对没接上，比数字变化更值得警惕。
# ==============================================================================
set -uo pipefail
cd "$(dirname "$0")" || exit 1

# Windows 上 Python 默认按系统代码页（cp936）写 stdout，中文断言与前缀解析都会失真；
# 文件读写也默认跟随 locale。这里统一成 UTF-8，保证换机器后行为一致。
# 入口脚本自身也做了兜底（见各文件的「控制台编码」段），此处是第二道。
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8

PY=${PYTHON:-python}
OUT=_repro
mkdir -p "$OUT"

PASS=0
FAIL=0
SKIP=0

GFF=rice/edta_out_multi/Oryza_sativa.Chr1_Chr12.fa.mod.EDTA.intact.gff3
MSU=ref/msu/osa1_r7.asm.repeat_masked.gff3.gz
FA1=rice/Oryza_sativa.Chr1.fa
FA12=rice/Oryza_sativa.Chr12.fa

check () {   # check <已入库> <新生成> [显示名]
    local want=$1 got=$2 label=${3:-$(basename "$1")}
    if [ ! -s "$got" ]; then
        printf "  [FAIL] %-26s 未生成（看 %s 下的日志）\n" "$label" "$OUT"
        FAIL=$((FAIL + 1))
    elif diff -q "$want" "$got" >/dev/null 2>&1; then
        printf "  [ OK ] %-26s 逐字节一致\n" "$label"
        PASS=$((PASS + 1))
    else
        printf "  [FAIL] %-26s 有差异\n" "$label"
        printf "         diff %s %s\n" "$want" "$got"
        FAIL=$((FAIL + 1))
    fi
}

step () { printf "\n== %s\n" "$1"; }

# ---------------------------------------------------------------- 0. 冒烟
step "0/3  te_check 能否直接跑（不需要 FASTA）"
if $PY te_check.py --edta-dir chr1_edta_out > "$OUT/smoke.log" 2>&1 \
   && grep -q "^C1 " "$OUT/smoke.log"; then
    printf "  [ OK ] %-26s 输出正常（%s 行）\n" "chr1 无 FASTA 模式" \
        "$(wc -l < "$OUT/smoke.log")"
    PASS=$((PASS + 1))
else
    printf "  [FAIL] %-26s 看 %s\n" "chr1 无 FASTA 模式" "$OUT/smoke.log"
    FAIL=$((FAIL + 1))
fi

# ---------------------------------------------------------------- 1. 结构台账
step "1/3  LTR 结构 → MSU 台账 → A/B/C 分级"
$PY TE_case/verify_ltr_structure.py --gff "$GFF" --fasta "$FA1" "$FA12" \
    --tsv "$OUT/ltr_structure.tsv" > "$OUT/verify_ltr.log" 2>&1
check TE_case/ltr_structure.tsv "$OUT/ltr_structure.tsv"

$PY TE_case/link_msu.py --gff "$GFF" --msu "$MSU" --kind ltr \
    --tsv "$OUT/ltr_msu_ledger.tsv" > "$OUT/link_ltr.log" 2>&1
check TE_case/ltr_msu_ledger.tsv "$OUT/ltr_msu_ledger.tsv"

$PY TE_case/grade_ltr_abc.py --structure "$OUT/ltr_structure.tsv" \
    --ledger "$OUT/ltr_msu_ledger.tsv" --tsv "$OUT/ltr_abc.tsv" > "$OUT/grade.log" 2>&1
check TE_case/ltr_abc.tsv "$OUT/ltr_abc.tsv"

# ---------------------------------------------------------------- 2. TIR 链路
step "2/3  TIR 结构 → 强中档台账 → 家族聚类"
$PY TE_case/verify_tir_structure.py --gff "$GFF" --fasta "$FA1" "$FA12" \
    --tsv "$OUT/tir_structure.tsv" > "$OUT/verify_tir.log" 2>&1
check TE_case/tir_structure.tsv "$OUT/tir_structure.tsv"

$PY TE_case/link_msu.py --gff "$GFF" --msu "$MSU" --kind tir \
    --anchor-tsv "$OUT/tir_structure.tsv" --tier 强,中 \
    --tsv "$OUT/tir_msu_ledger.tsv" > "$OUT/link_tir.log" 2>&1
check TE_case/tir_msu_ledger.tsv "$OUT/tir_msu_ledger.tsv"

$PY TE_case/tir_family_cluster.py --gff "$GFF" --ledger "$OUT/tir_msu_ledger.tsv" \
    --tsv "$OUT/tir_family.tsv" > "$OUT/family.log" 2>&1
check TE_case/tir_family.tsv "$OUT/tir_family.tsv"

# ---------------------------------------------------------------- 3. 横向对比
step "3/3  六数据集横向对比"
$PY bench_multi.py \
    "chr1=chr1_edta_out" \
    "col1m=bench_archive/col1m.fa.mod" \
    "alyrata=bench_archive/alyrata.fa.mod" \
    "rice_pre=rice/edta_out_preflight" \
    "rice_full=rice/edta_out" \
    "rice_multi=rice/edta_out_multi" \
    --tsv "$OUT/bench_six.tsv" > "$OUT/bench.log" 2>&1
check TE_case/bench_six.tsv "$OUT/bench_six.tsv"

# ---------------------------------------------------------------- 4. 四份报告
if [ "${1:-}" = "--with-reports" ]; then
    step "4/4  四份 te_check JSON（需基因组 FASTA）"
    if [ ! -s "$FA1" ] || [ ! -s "$FA12" ] || [ ! -s ref/TAIR10.Chr1.fa ]; then
        printf "  [SKIP] 缺基因组 FASTA，四份报告未比对\n"
        printf "         获取方式见 docs/复现指南.md「二、序列从哪来」\n"
        SKIP=$((SKIP + 4))
    else
        $PY te_check.py --edta-dir chr1_edta_out --fasta ref/TAIR10.Chr1.fa \
            --json "$OUT/chr1_report.json" > "$OUT/rep_chr1.log" 2>&1
        check chr1_report.json "$OUT/chr1_report.json" "chr1_report.json"

        $PY te_check.py --edta-dir rice/edta_out --fasta "$FA1" \
            --json "$OUT/full_report.json" > "$OUT/rep_full.log" 2>&1
        check rice/full_report.json "$OUT/full_report.json" "full_report.json"

        # 预飞用的 3 Mb 序列不入库，按 run_bench_rice.sh 的 PREFLIGHT 规则重建
        {
            echo ">Oryza_sativa.Chr1.pre3000000"
            awk 'NR>1' "$FA1" | tr -d '\n' | cut -c1-3000000 | fold -w 60
        } > "$OUT/pre3000000.fa"
        $PY te_check.py --edta-dir rice/edta_out_preflight --fasta "$OUT/pre3000000.fa" \
            --json "$OUT/preflight_report.json" > "$OUT/rep_pre.log" 2>&1
        check rice/preflight_report.json "$OUT/preflight_report.json" "preflight_report.json"

        $PY te_check.py --edta-dir rice/edta_out_multi --fasta "$FA1" "$FA12" \
            --json "$OUT/multi_report.json" > "$OUT/rep_multi.log" 2>&1
        check rice/multi_report.json "$OUT/multi_report.json" "multi_report.json"
    fi
fi

# ---------------------------------------------------------------- 汇总
printf "\n%s\n" "------------------------------------------------------------"
printf "一致 %d 项   不一致 %d 项   跳过 %d 项\n" "$PASS" "$FAIL" "$SKIP"
if [ "$FAIL" -eq 0 ]; then
    printf "全部可复现。\n"
else
    printf "有产物对不上：先看 diff，再确认这是判据变更还是数据变更。\n"
fi
[ "$SKIP" -gt 0 ] && printf "（跳过项需要基因组 FASTA，见 docs/复现指南.md）\n"
printf "%s\n" "------------------------------------------------------------"
[ "$FAIL" -eq 0 ]
