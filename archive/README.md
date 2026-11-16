# archive/ —— 已被取代的运行脚本

保留这些脚本不是为了复用，而是因为**每一条失败链都是判据检验结论的一部分**。
报告里"C4 的分母算错"这类教训，源头就藏在这几版脚本的注释里。

| 脚本 | 做过什么 | 结局 | 被谁取代 |
|---|---|---|---|
| `run_bench.sh` | 第一次跑 Chr1 + col1m + alyrata，`--step all --anno 1` 一次跑完 | 碰到 LTR 阶段双 `.mod` 前缀失败 | → run_bench2.sh |
| `run_bench2.sh` | 重试策略：失败后重跑同一阶段 | **原理上不可能成功**——`EDTA_raw.pl` 用 `-s "$genome.LTRlib.fa"` 判断是否跳过，而该文件因改名而永远不存在，于是每轮都重跑 retriever、每轮都产出双前缀 | → run_bench3.sh |
| `run_bench3.sh` | 失败后先改名（`X.mod.mod.*` → `X.mod.*`）再续跑 | 一次通过，产出 col1m / alyrata 两份可比数据 | → `../run_bench_rice.sh` |

`../run_bench_rice.sh` 已把 `run_bench3.sh` 的改名修复内置，并参数化到任意基因组
（`GENOME= TAG= SPECIES= OUTDIR= DEST= THREADS=`），默认值仍是水稻 Chr1。

## 一处必须说明的坑

这三个脚本的注释里，关于 `.mod.mod` 成因的说明**都不可信**：

- `run_bench3.sh` v1：写成 `LTR_retriever.py:70` 的 `glob` 拼写缺陷 → 无源码支持
- 之后的修订：写成 EDTA issue #599 的「序列名超过 rmblastn 约 50 字符」→ 与本场景
  4–6 字符的 seqid 矛盾

两版都已被源码推翻，现行记录（含排除过程与未验证的假设）在
`../run_bench_rice.sh` 的 `repair_names` 注释、以及
`../TE_case/报表_水稻跨科判据检验.md` 第 5.3 节。

**改造名修复本身是有效的**（四个基因组实证），错的只是对原因的猜测。
