#!/usr/bin/env bash
# ==============================================================================
# 水稻数据集：IRGSP-1.0 Chr1 的 EDTA 运行（te_check 的第四个数据集）
#
# 为什么是水稻：此前三份数据（TAIR10 Chr1、Col_1M、Alyrata_2M）全部属于
# 十字花科拟南芥属。C1/C3/C4 的结论只在双子叶、TE 含量约 20% 的基因组上成立过。
# 水稻是单子叶、TE 含量约 35–40%、且有 IRGSP-1.0 这样的完整参考组装，
# 是检验"判据能否跨科移植"的最小充分案例。
#
# 参数与前三个数据集的唯一差异是 --species others → Rice。
# 依据（EDTA.pl 源码）：--species 只是转发给 EDTA_raw.pl，用于 TIR 候选鉴定；
# 它不影响 LTR / LINE / SINE / Helitron 的发现，也不影响过滤与注释阶段。
# 单子叶必须显式设为 Rice，否则 TIR_Learner 按双子叶参数跑，TIR 候选会系统性偏少。
#
# 刻意不提供 --curatedlib。EDTA 自带 database/rice7.0.1.liban.{LTR,LINE,SINE,TIR,Helitron}，
# 一旦挂上就不再是 de novo，与前三份数据不可比。想换成"挂库"路径请用
#   USE_CURATED=1 bash run_bench_rice.sh
# 但那份结果不能和前三份放同一张表里比较。
#
# 用法（在 WSL 的终端里执行，不要在 Windows 侧跑）：
#   bash /mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/run_bench_rice.sh
#   THREADS=16 bash ...                       # 指定线程
#   PREFLIGHT=1 bash ...                      # 先截 3 Mb 试跑
#
# 脚本已参数化，可跑任意基因组（默认值 = 水稻 Chr1）：
#   GENOME= TAG= SPECIES= OUTDIR= DEST= THREADS=
# 例如重跑拟南芥 TAIR10 Chr1 作为基准：
#   GENOME=/mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/ref/TAIR10.Chr1.fa \
#   TAG=Chr1 SPECIES=others OUTDIR=$HOME/te_bench_chr1 \
#   DEST=/mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/chr1_edta_out \
#   THREADS=16 bash /mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06/run_bench_rice.sh
#
# 耗时预期：比拟南芥 Chr1（30.4 Mb）慢，主要慢在 TEsorter 与 RepeatMasker——
# 水稻 TE 密度高，进入分类的候选序列数量远多于拟南芥。脚本每一阶段都会打时间戳。
# ==============================================================================
set -uo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate edta || { echo "无法激活 conda 环境 edta"; exit 1; }

WINDIR="/mnt/c/Users/Sean/WorkBuddy/2026-11-13-13-40-06"
# 以下五项均可用环境变量覆盖，默认值即水稻 Chr1 的已验证配置。
SRC="${GENOME:-$WINDIR/rice/Oryza_sativa.Chr1.fa}"
OUT="${OUTDIR:-$HOME/te_bench_rice}"
DEST="${DEST:-$WINDIR/rice/edta_out}"   # 分析用小文件回传到 Windows 侧
THREADS="${THREADS:-8}"
TAG="${TAG:-Oryza_sativa.Chr1}"
SPECIES="${SPECIES:-Rice}"              # Rice | Maize | others
G="$TAG.fa"

# ---------------------------------------------------------------- 预检
echo "=============================================================="
echo " $TAG EDTA 运行（species=$SPECIES）"
echo "=============================================================="
echo " EDTA       : $CONDA_PREFIX"
echo " 线程       : $THREADS"
echo " 输出目录   : $OUT"
echo " 回传目录   : $DEST"
echo

for tool in EDTA.pl LTR_retriever TEsorter RepeatMasker; do
    if command -v "$tool" >/dev/null 2>&1; then
        printf " %-14s ✅ %s\n" "$tool" "$(command -v "$tool")"
    else
        printf " %-14s ❌ 未找到\n" "$tool"
        MISSING=1
    fi
done
if [ -n "${MISSING:-}" ]; then
    echo
    echo "缺工具，停止。请确认 conda 环境 edta 完整。"
    exit 1
fi
echo

if [ ! -s "$SRC" ]; then
    echo "❌ 找不到输入序列：$SRC"
    echo "   该文件应在 Windows 侧预先下载（Ensembl Plants release-57, IRGSP-1.0）。"
    exit 1
fi
printf " 输入序列   : %s\n" "$SRC"
printf " 序列总长   : %s bp\n" \
    "$(awk 'NR>1{n+=length($0)} END{print n}' "$SRC")"
echo

AVAIL=$(df -Pk "$HOME" | awk 'NR==2{printf "%.1f", $4/1048576}')
echo " 家目录可用 : ${AVAIL} GB（43 Mb 基因组的 EDTA 中间产物通常 < 5 GB）"
if [ "${AVAIL%%.*}" -lt 5 ] 2>/dev/null; then
    echo " ⚠ 空间偏紧，继续前建议清理。仍将尝试运行。"
fi
echo

mkdir -p "$OUT" || { echo "无法创建 $OUT"; exit 1; }
cd "$OUT" || exit 1

# ---------------------------------------------------------------- 预飞模式
# 水稻 TE 密度高，LTR 候选多，TEsorter/RepeatMasker 阶段可能比拟南芥慢数倍。
# 全染色体跑几小时后才发现流程问题代价太大，所以支持先截一段跑通全链路。
# 预飞结果不用于判据检验（小片段上建库质量差），只用来确认 EDTA 不中途失败。
#   PREFLIGHT=1 bash run_bench_rice.sh
if [ "${PREFLIGHT:-0}" = "1" ]; then
    SPAN="${PREFLIGHT_SPAN:-3000000}"
    TAG="${TAG}.pre${SPAN}"
    G="$TAG.fa"
    DEST="${DEST}_preflight"                 # 与全量结果分开，避免 find_files 混读
    echo "--------------------------------------------------------------"
    echo " 预飞模式：只取 Chr1 前 ${SPAN} bp（坐标与原染色体一致）"
    echo " 产物前缀 $G —— 不会覆盖全染色体的结果"
    echo " 回传目录改为 $DEST"
    echo "--------------------------------------------------------------"
    echo
    if [ ! -s "$G" ]; then
        {
            echo ">${TAG}"
            awk 'NR>1' "$SRC" | tr -d '\n' | cut -c1-"$SPAN" | fold -w 60
        } > "$G" || { echo "截取失败"; exit 1; }
        printf " 已生成 %s（%s bp）\n\n" "$G" \
            "$(awk 'NR>1{n+=length($0)} END{print n}' "$G")"
    fi
fi

[ -s "$G" ] || cp -f "$SRC" "$G"

# ---------------------------------------------------------------- 双前缀修复
# 现象：LTR 阶段结束后，LTR_retriever 的全部产物（pass.list / defalse / LTRlib.fa）
# 落在 X.fa.mod.mod.* 前缀下，而 EDTA_raw 按 X.fa.mod.* 查找
# → LTR.raw.fa 为空 → "Raw LTR results not found"。
#
# 【成因已查明，2026-11-16 第三次核查，对着实际安装的 LTR_retriever 源码确认】
#
# 一句话：EDTA.pl 已把基因组归一化成 X.fa.mod，LTR_retriever 又归一化了第二次，
#         而它的输出路径正是 "$genome.mod"，于是得到 X.fa.mod.mod。
#
# 完整链条（三处源码，逐行可复现）：
#
#  1) EDTA.pl v2.2.2:365-398 —— 清洗 seqid、非 ATGC 转 N，写出 X.fa.mod，
#     随后无条件 $genome = "$genome.mod"。此后 $genome 就是 X.fa.mod。
#     这一步是必要且正确的。
#
#  2) EDTA_raw.pl v2.2.2:440 —— 调用 LTR_retriever 时没传 -convert_seq_name：
#       LTR_retriever -genome $genome -inharvest ... -noanno ... -repeatmasker ...
#     而 LTR_retriever 该变量默认为 1，于是第二次归一化被打开。
#
#  3) LTR_retriever Init 段（实际安装版本 402-415 行）——
#       if ($convert_seq_name){
#           `perl $seqid_codec encode_fasta $genome $genome.mod.seqid.map $genome.mod --prefix _J`;
#           if (-s "$genome.mod"){ $index="$genome.mod"; $genome="$genome.mod"; ... }
#       }
#     encode_fasta 的输出路径是 "$genome.mod"，此时 $genome 已是 X.fa.mod，
#     于是写出 X.fa.mod.mod；紧随的判断只检查文件存在，而 encode_fasta
#     无论是否需要编码都会写出该文件（日志里 "Genome seq IDs cleaned
#     (no encoding needed)" 走的就是这条），所以采纳是无条件的。
#     后续 pass.list / defalse / LTRlib.fa 全部带双前缀。
#
# 为什么与序列名长短无关：对 >Chr1 这种 4 字符 ID，encode_fasta 只是"清洗"（无变化），
# 但文件照写、照样被采纳。所以水稻(>OsChr1) 与 Chr1(>Chr1) 都复现。
# EDTA issue #599 那句"序列名超过 rmblastn 约 50 字符上限"为何套不上：那是长名
# 确实需要编码时暴露的**同一段代码**，成因同源、触发条件不同。
#
# 为什么 v2.3.0 不再出现（已核实 v2.3.0:410 与 master:419）：调用里补了
#   -convert_seq_name 0
# 即明确关掉这次多余的二次归一化。所以"升级可根治"方向对，但此前给的理由（长头）错。
#
# 为什么暂不升级、也不给 conda 里的 EDTA_raw.pl 打这一行补丁：五份对比数据全在
# v2.2.2 下跑的，任何改变环境行为的操作都会让新旧不可比。改名修复不改环境。
#
# 处置：把 X.fa.mod.mod.* 复制成 X.fa.mod.*（纯改名，不改内容）。这些文件只含
# 序列 ID 与坐标，不含文件名，改名不改语义。已验证有效：水稻 9 个文件、
# Chr1 9 个文件，改后建库一次通过。
repair_names () {
    local ltrdir="$G.mod.EDTA.raw/LTR"
    [ -d "$ltrdir" ] || { echo "  （尚无 LTR 目录，跳过修复）"; return 0; }
    local n=0 f base new
    for f in "$ltrdir/$G.mod.mod."*; do
        [ -e "$f" ] || continue
        base=$(basename "$f")
        new=$(printf '%s' "$base" | sed 's/\.mod\.mod\./.mod./')
        [ "$base" = "$new" ] && continue
        cp -f "$f" "$ltrdir/$new" && n=$((n+1))
    done
    echo "  修复：复制 $n 个文件（X.mod.mod.* → X.mod.*）"
    local k
    for k in pass.list defalse LTRlib.fa; do
        if [ -s "$ltrdir/$G.mod.$k" ]; then
            printf "   %-12s %8s 字节 ✅\n" "$k" "$(stat -c%s "$ltrdir/$G.mod.$k")"
        else
            printf "   %-12s 缺失或为空\n" "$k"
        fi
    done
}

# ---------------------------------------------------------------- 阶段一：建库
if [ -s "$G.mod.EDTA.intact.gff3" ]; then
    echo "[跳过] intact.gff3 已存在，直接进入回传阶段"
else
    echo "=============================================================="
    echo "[$(date '+%F %T')] 阶段一：--step all（过滤 + 建库）"
    echo "=============================================================="

    CURATED_ARGS=()
    if [ "${USE_CURATED:-0}" = "1" ]; then
        DB="$CONDA_PREFIX/share/EDTA/database"
        # EDTA 不含单文件水稻库，而是按类别拆成 rice7.0.1.liban.{LTR,LINE,SINE,TIR,Helitron}
        if compgen -G "$DB/rice7.0.1.liban.*" >/dev/null; then
            COMB="$OUT/rice7.0.1.liban.fa"
            [ -s "$COMB" ] || cat "$DB"/rice7.0.1.liban.* > "$COMB"
            CURATED_ARGS=(--curatedlib "$COMB")
            echo "[注意] 已挂 curatedlib（$(grep -c '^>' "$COMB") 条），"
            echo "       本份结果不可与前三份 de novo 数据同表比较"
        else
            echo "[警告] 未找到 $DB/rice7.0.1.liban.*，退回纯 de novo"
        fi
    fi
    echo

    EDTA.pl --genome "$G" --species "$SPECIES" --step all \
            --overwrite 0 --threads "$THREADS" ${CURATED_ARGS[@]+"${CURATED_ARGS[@]}"} 2>&1 | tail -15

    # 首次运行若踩到双前缀缺陷，修复后重跑一次（复用中间结果）
    if [ ! -s "$G.mod.EDTA.raw/LTR/$G.mod.LTR.raw.fa" ]; then
        echo
        echo "[$(date '+%T')] LTR.raw.fa 为空，尝试修复双前缀后重跑"
        repair_names
        echo
        EDTA.pl --genome "$G" --species "$SPECIES" --step all \
                --overwrite 0 --threads "$THREADS" ${CURATED_ARGS[@]+"${CURATED_ARGS[@]}"} 2>&1 | tail -15
    fi

    if [ ! -s "$G.mod.EDTA.raw/LTR/$G.mod.LTR.raw.fa" ]; then
        echo
        echo "❌ LTR 仍未生成。LTR 目录现状："
        ls -lh "$G.mod.EDTA.raw/LTR" 2>/dev/null \
            | awk 'NR>1{printf "   %-58s %s\n",$9,$5}' | tail -25
        echo
        echo "已停止。中间结果保留在 $OUT，修复后可重跑本脚本（--overwrite 0 会复用）。"
        exit 1
    fi
    printf "\n[$(date '+%T')] 建库完成，LTR.raw.fa %s 字节\n" \
        "$(stat -c%s "$G.mod.EDTA.raw/LTR/$G.mod.LTR.raw.fa")"

    # ------------------------------------------------------------ 阶段二：注释
    echo
    echo "=============================================================="
    echo "[$(date '+%F %T')] 阶段二：--step anno（全基因组注释）"
    echo "=============================================================="
    EDTA.pl --genome "$G" --species "$SPECIES" --step anno --anno 1 \
            --threads "$THREADS" 2>&1 | tail -12
fi

if [ ! -s "$G.mod.EDTA.intact.gff3" ]; then
    echo
    echo "❌ 未产出 intact.gff3，停止。"
    exit 1
fi

# ---------------------------------------------------------------- 回传
echo
echo "=============================================================="
echo "[$(date '+%F %T')] 回传分析用文件到 Windows 侧"
echo "=============================================================="
mkdir -p "$DEST"
for suffix in EDTA.intact.gff3 EDTA.TEanno.gff3 EDTA.TEanno.sum EDTA.TElib.fa; do
    f="$G.mod.$suffix"
    if [ -s "$f" ]; then
        cp -f "$f" "$DEST/" && printf "  %-26s %10s 字节 ✅\n" "$suffix" "$(stat -c%s "$f")"
    fi
done
# 分类表需保留 <prefix>.EDTA.raw/<lib>/ 两层结构，te_check 按该模式 glob
if [ -d "$G.mod.EDTA.raw" ]; then
    for f in "$G.mod.EDTA.raw"/*/*.rexdb.cls.tsv; do
        [ -e "$f" ] || continue
        sub=$(basename "$(dirname "$f")")
        mkdir -p "$DEST/$G.mod.EDTA.raw/$sub"
        cp -f "$f" "$DEST/$G.mod.EDTA.raw/$sub/"
        printf "  分类表 %-20s → %s\n" "$sub" "$(basename "$f")"
    done
fi
# WSL 侧路径 → Windows 可见路径，供 Windows 侧 Python 直接读
WIN_HOME=$(wslpath -w "$HOME"   2>/dev/null | sed 's#\\#/#g')
WWIN=$(wslpath -w "$WINDIR" 2>/dev/null | sed 's#\\#/#g')

# 把物种信息一并留档，防止日后忘记这份数据用的什么参数
cat > "$DEST/RUN_INFO.txt" <<EOF
tag        : $TAG
species    : $SPECIES
parameters : --step all, then --step anno --anno 1; --overwrite 0; --threads $THREADS
curatedlib : ${USE_CURATED:-0}
input      : Ensembl Plants release-57, Oryza_sativa.IRGSP-1.0.dna.chromosome.1
length     : $(awk 'NR>1{n+=length($0)} END{print n}' "$SRC") bp
edta       : $CONDA_PREFIX
outdir     : $OUT
win_path   : $WWIN/rice/$(basename "$DEST")
finished   : $(date '+%F %T')
EOF

# ---------------------------------------------------------------- 归档前三个数据集
# EDTA 产物里真正用于分析的只有 intact.gff3、TEanno.sum、TElib.fa 和各分类表，
# 都是小文件。镜像到 Windows 侧之后，四数据集对比不再依赖"WSL 此刻是否在运行"。
ARCH="$WINDIR/bench_archive"
archived=0
for d in "$HOME/te_bench" "$HOME/te_bench2"; do
    [ -d "$d" ] || continue
    for gff in "$d"/*.EDTA.intact.gff3; do
        [ -s "$gff" ] || continue
        pre=$(basename "$gff" .EDTA.intact.gff3)
        dst="$ARCH/$pre"
        mkdir -p "$dst"
        cp -f "$gff" "$dst/" 2>/dev/null
        for sfx in EDTA.TEanno.sum EDTA.TElib.fa; do
            [ -s "$d/$pre.$sfx" ] && cp -f "$d/$pre.$sfx" "$dst/"
        done
        if [ -d "$d/$pre.EDTA.raw" ]; then
            for c in "$d/$pre.EDTA.raw"/*/*.rexdb.cls.tsv; do
                [ -e "$c" ] || continue
                sub=$(basename "$(dirname "$c")")
                mkdir -p "$dst/$pre.EDTA.raw/$sub"
                cp -f "$c" "$dst/$pre.EDTA.raw/$sub/"
            done
        fi
        printf "  归档 %-44s → %s\n" "$pre" "$dst"
        archived=$((archived+1))
    done
done
BENCH_ARGS=()
if [ "$archived" -gt 0 ]; then
    echo
    echo "  已归档 $archived 份历史运行到 $ARCH"
    for d in "$ARCH"/*/; do
        [ -d "$d" ] || continue
        name=$(basename "${d%/}")
        label=$(printf '%s' "$name" | sed 's/\.fa\.mod$//; s/\.fa$//')
        BENCH_ARGS+=("\"${label}=$WWIN/bench_archive/${name}\"")
        printf "    %-10s %s\n" "$label" "$WWIN/bench_archive/$name"
    done
else
    echo "  （未找到 te_bench / te_bench2，跳过归档）"
fi

echo
echo "=============================================================="
printf " intact 非注释行 : %s\n" "$(grep -vc '^#' "$G.mod.EDTA.intact.gff3")"
printf " 分类表数量      : %s\n" "$(ls "$G.mod.EDTA.raw"/*/*.rexdb.cls.tsv 2>/dev/null | wc -l)"
printf " 回传目录        : %s\n" "$DEST"
echo "=============================================================="
echo "RICE_DONE"

# 生成可直接执行的对比命令。用 Windows 侧归档路径而不是 WSL 路径，
# 这样即使关掉 WSL 也能复现对比。
# 标签与路径都从本次实际取值派生。早先这里写死成 "rice=$WWIN/rice/$(basename "$DEST")"，
# 跑 Chr1 时会拼出 ../rice/chr1_edta_out（不存在），Tag 和 DEST 一旦改动就失真。
THIS_WIN=$(wslpath -w "$DEST" 2>/dev/null | sed 's#\\#/#g')
[ -n "$THIS_WIN" ] || THIS_WIN="$DEST"
LABEL="${LABEL:-$(printf '%s' "$TAG" | sed 's/\.fa$//')}"
if [ ${#BENCH_ARGS[@]} -gt 0 ]; then
    CMD="cd \"$WWIN\"; python bench_multi.py ${BENCH_ARGS[*]} \"${LABEL}=$THIS_WIN\""
    mkdir -p "$ARCH"
    printf '%s\n' "$CMD" > "$ARCH/run_bench_cmd.txt"
    echo
    echo "下一步：在 Windows 侧（PowerShell 或 Git Bash）执行"
    echo
    echo "  $CMD"
    echo
    echo "  （同一条命令已写入 $ARCH/run_bench_cmd.txt）"
else
    echo
    echo "下一步：python bench_multi.py \"${LABEL}=$THIS_WIN\""
fi
echo
