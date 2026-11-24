#!/usr/bin/env bash
# ############################################################
# 【已废弃】下载的库经实测不含任何植物序列，已改用 Araport11 参照。
# 详见 run_validation.sh 顶部说明与 EDTA_vs_Araport11_评估报告.md
# ############################################################
# ============================================================
# 稳健下载 Dfam 共识库（15.8 MB）
#
# 场景：wget 卡住不写盘（dfam.org 从国内访问常卡在 TLS/IPv6）
# 对策：强制 IPv4 + 连接超时 + 有限重试 + 多版本回退 + 断点续传
#
# 运行：bash get_dfam_lib.sh
# ============================================================
set -uo pipefail

cd ~/edta_run
mkdir -p dfam_lib
cd dfam_lib

LOG=../get_dfam_lib.log
: > "$LOG"

say() { echo "$@" | tee -a "$LOG"; }

say "=== 1. 连通性测试（强制 IPv4，10 秒超时）==="
if curl -4 -sSI --connect-timeout 10 https://www.dfam.org/ > /tmp/dfam_head.txt 2>&1; then
    head -3 /tmp/dfam_head.txt | tee -a "$LOG"
    say "IPv4 可达"
else
    say "IPv4 连接失败，内容如下："
    cat /tmp/dfam_head.txt | tee -a "$LOG"
    say ""
    say "→ 若提示 Could not resolve host，说明 DNS 有问题"
    say "→ 若提示 Connection timed out，说明网络层不通，需要代理"
    say "→ 若提示 SSL/TLS 相关，试：export CURL_SSL_BACKEND=openssl"
fi

say ""
say "=== 2. 下载（按版本从新到旧回退）==="
URLS=(
    "https://www.dfam.org/releases/Dfam_3.9/families/Dfam-RepeatMasker.lib.gz"
    "https://www.dfam.org/releases/Dfam_3.8/families/Dfam-RepeatMasker.lib.gz"
    "https://www.dfam.org/releases/Dfam_3.7/families/Dfam-RepeatMasker.lib.gz"
)

OK=0
for u in "${URLS[@]}"; do
    say "--- 尝试: $u"
    if curl -4 -L --fail --retry 3 --retry-delay 2 \
             --connect-timeout 15 --speed-time 60 --speed-limit 1024 \
             -C - -o Dfam-RepeatMasker.lib.gz \
             --progress-bar "$u" 2>&1 | tee -a "$LOG"; then
        say "下载成功"
        OK=1
        break
    fi
    say "本地址失败，换下一个"
done

if [ "$OK" -ne 1 ] || [ ! -s Dfam-RepeatMasker.lib.gz ]; then
    say ""
    say "ERROR: 三个地址都失败。请检查网络，或手动用浏览器下载后放到:"
    say "  ~/edta_run/dfam_lib/Dfam-RepeatMasker.lib.gz"
    exit 1
fi

say ""
say "=== 3. 解压与校验 ==="
BYTES=$(stat -c%s Dfam-RepeatMasker.lib.gz)
say "压缩包大小: $BYTES 字节（应约 16,000,000）"

if ! gunzip -kf Dfam-RepeatMasker.lib.gz 2>>"$LOG"; then
    say "ERROR: 解压失败，文件可能不完整。删掉重跑本脚本（会自动续传）。"
    exit 1
fi

N=$(grep -c '^>' Dfam-RepeatMasker.lib)
say "家族数: $N"

say ""
say "=== 4. 植物类群覆盖检查 ==="
for k in Viridiplantae Embryophyta Arabidopsis Oryza Brassicaceae; do
    printf "  %-16s %s\n" "$k" "$(grep -c "@$k" Dfam-RepeatMasker.lib || true)" | tee -a "$LOG"
done

VP=$(grep -c "@Viridiplantae" Dfam-RepeatMasker.lib || true)
say ""
if [ "$VP" -gt 0 ]; then
    say "OK: 库中有 $VP 条植物类群记录，可以用于交叉验证。"
    say "下一步：cd ~/edta_run && bash /path/to/te-audit/run_validation.sh"
else
    say "WARN: 未发现 @Viridiplantae 标注。Dfam 的类群标注格式可能不同，"
    say "      可以先用下面的命令看看实际标注长什么样："
    say "      grep -o '@[A-Za-z ]*' Dfam-RepeatMasker.lib | sort | uniq -c | sort -rn | head -20"
fi
