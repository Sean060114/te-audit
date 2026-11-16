# EDTA 中 1 个完整 LTR 被静默丢弃的溯源报告

- 运行目录：`~/edta_run`，基因组 `chr1.fa.mod`
- 现象：`chr1.fa.mod.pass.list` 有 **41** 条完整 LTR，但 `chr1.fa.mod.LTR.intact.raw.fa` 只有 **40** 条，`EDTA.TElib.fa` 中也没有这个元件
- 涉事元件：`1:17528667..17536138`（7,472 bp，`LTR/unknown`，`ltr_identity=0.9866`）

---

## 一、结论（先看这段）

丢序列的不是提取步骤，也不是 TEsorter 过滤，而是 **`bin/rename_LTR_skim.pl` 第 51 行的静默 `next`**：

```perl
my $anno = "unknown";
next unless defined $anno{$id};   # ← 精确字符串匹配失败就整条丢掉，不报错、不记录
```

`%anno` 的键来自 `chr1.fa.mod.defalse` 第 1 列，名字写作**升序** `1:17528667..17536138`；
而待查的 fasta 头部名字来自 `pass.list`，写作**降序** `1:17536138..17528667`。
两者只是坐标顺序反了，字符串对不上 → 该元件被丢弃。

因为 `$genome.LTR.intact.fa.ori.rmlist`（本该记录“被排除者”的名单）最终是 **0 字节**，
下游没有任何日志能发现少了这一条，所以表现为“41 进、40 出”且无告警。

---

## 二、逐步计数定位

`EDTA_raw.pl` 第 445–462 行的实际链路：

| 步骤 | 命令 | 条数 |
|---|---|---|
| 1 | `awk ... pass.list \| call_seq_by_list.pl - -C` → `LTR.intact.fa.ori` | 41 → |
| 2 | `rename_LTR_skim.pl LTR.intact.fa.ori defalse` → 覆盖 `.ori` | → **40** ← 丢在这里 |
| 3 | `mdust` → `.ori.dusted` | 40 |
| 4 | `cleanup_tandem` → `.ori.dusted.cln` | 40 |
| 5 | `TEsorter` + `cleanup_misclas` → `LTR.intact.raw.fa` | 40 |

配套证据：

- `pass.list` 非注释行 = **41**
- `.ori` 的 `>` 头 = **40**，且用 `17528667|17536138` 检索**无任何命中**
- 该元件全序列在 `.ori` 中既非正链也非反链存在，也无可匹配的前 40 bp 片段 —— 说明是整条丢弃，不是被裁剪
- `.ori.rmlist` = **0 字节**（`output_by_list -ex` 没有排除任何东西，因为 `.ori` 已经不含它了）

> 注：`LTR.intact.fa.ori` 在流程中被 `mv` 覆盖了两次，磁盘上看到的是 **改名后** 的结果，
> 所以第 2 步才是真正断点。`call_seq_by_list.pl` 本身对反序坐标有正确处理
> （第 110–115 行：`if ($start>$stop)` 自动交换并置 `$strand="-"`），不是它的问题。

---

## 三、为什么只有这一条中招

`pass.list` 41 条中，25 条升序、16 条降序。按字符串精确比对 `defalse` 第一列：

```
pass.list 有、defalse 无（精确匹配）：
    1:17536138..17528667   降序        ← 唯一一条
```

其余 15 条降序元件都能匹配上，因为 `defalse` 里它们**也是降序**写的
（如 `1:4411120..4406006`、`1:22105329..22100644`）。`defalse` 共 246 个唯一名字，其中 54 条是降序。

差别在链方向的判定依据：

| 元件 | domain_strand | strand | strand_evidence | defalse 写法 |
|---|---|---|---|---|
| `1:4411120..4406006` | `-` | `-` | `domain` | 降序 ✓ 一致 |
| `1:11152835..11148091` | `-` | `-` | `domain+PBS+PPT` | 降序 ✓ 一致 |
| **`1:17536138..17528667`** | **`?`** | `-` | **`PBS`** | **升序 ✗ 不一致** |

它是 41 条中**唯一一条蛋白结构域完全没打上（`domain_strand=?`）、链方向纯靠 PBS 推断**的元件
（GFF3 里也留着 `strand_evidence=PBS`）。

`defalse` 中它还出现了**两条重复候选**：

```
1254: 1:17528667..17536138  pass  motif:TGCA  TSD:ATCTC  0.9866  ?  unknown  514243
1261: 1:17528667..17536138  pass  motif:TGCA  TSD:ATCTC  0.9867  ?  unknown  513263
（LTR coordinates: 17528667, 17529195, 17535610, 17536138 —— 长度 7473 / 7472 差 1 bp）
```

LTRharvest（`17528667 17536139 7473`）与 LTR_FINDER（`17528667 17536138 7472`）各报一次；
`pass.list` 保留了其中一条并在链判定后被写成降序，`defalse` 两条都保持升序。
即：**坐标书写方向在 `defalse` 与 `pass.list` 之间失去同步**，而这正好撞上精确字符串查表。

---

## 四、实际影响

**已出现的不一致：注释在，但库和序列里没有。**

- `chr1.fa.mod.LTR.intact.raw.gff3` 与最终 `chr1.fa.mod.EDTA.TEanno.gff3` **都保留了**这条注释：

  ```
  repeat_region          17528662  17536143  ID=repeat_region_34;Name=1:17536138..17528667;classification=LTR/unknown
  LTR_retrotransposon    17528667  17536138  ID=LTRRT_34;... ltr_identity=0.9866;method=structural;motif=TGCA;tsd=GAGAT;strand_evidence=PBS
  long_terminal_repeat   17528667  17529195  ID=lLTR_34   (529 bp)
  long_terminal_repeat   17535610  17536138  ID=rLTR_34   (529 bp)
  primer_binding_site    17535589  17535607  ID=PBS_34    (19 bp, Lys-CTT)
  ```

- 但 `LTR.intact.raw.fa` 无此序列 → `TEsorter` 从未见过它 → 家族分类停留在 `LTR/unknown`
- `EDTA.TElib.fa`（253 个家族）中用全序列检索：**完全相同的 0 个，完整包含 7,472 bp 的 0 个**

因此：

1. **基因组屏蔽/注释不受影响** —— `TEanno.gff3` 照样会遮蔽这 7.5 kb，`RepeatMasker` 走注释结果时能覆盖；
2. **全长 LTR 库缺一条** —— 建非冗余库、家族代表性序列、LTR 插入年龄分布（5' / 3' LTR 差异）、
   插入多态性分析都会漏掉这个 7,472 bp 的完整元件；
3. **家族统计偏差** —— 该元件没有成为任何家族的成员，`unknown` 类被低估 1 条。

---

## 五、修复建议

### 方案 A：打补丁（推荐，一次性解决同类问题）

在 `rename_LTR_skim.pl` 的查表处加坐标顺序回退：

```perl
my $anno = "unknown";
# 原：next unless defined $anno{$id};
unless (defined $anno{$id}) {
    # 坐标顺序可能被链判定翻转，尝试交换后重查
    if ($id =~ /^(.*):([0-9]+)\.\.([0-9]+)$/) {
        my $alt = "$1:$3..$2";
        $id = $alt if defined $anno{$alt};
    }
}
next unless defined $anno{$id};   # 仍找不到才真丢
```

补丁后 **必须重跑 LTR 模块**（`--overwrite 1 --sensitive 1` 或删掉 `*.EDTA.raw/LTR` 后重跑），
否则 `EDTA.TElib.fa` 不会更新。

### 方案 B：只补这一条（快速验证用）

```bash
# 1. 用升序坐标从基因组切出序列（注意是负链，需反向互补）
# 2. 以 >1:17536138..17528667#LTR/unknown 为头追加进 chr1.fa.mod.LTR.intact.raw.fa
# 3. 重新执行 TEsorter + cleanup_misclas + 重建 TElib
```

不建议手改 `EDTA.TElib.fa`，因为家族编号（`TE_xxxxxxxx`）与 `TEanno.gff3` 的映射会错位。

### 附：给上游报 issue 的要点

- EDTA 版本、LTR_retriever 版本
- 触发条件：`defalse` 与 `pass.list` 坐标书写方向不一致（本例由 `strand_evidence=PBS` 的
  元件引起，`domain_strand=?`）
- 后果：`rename_LTR_skim.pl:51` 静默丢弃，且 `LTR.intact.fa.ori.rmlist` 为空，无任何日志
- 建议：查表前统一坐标顺序，或丢弃时写入 `rmlist` 并 `warn`

---

## 六、可复用的自查脚本

对任意 EDTA 运行目录，快速找出所有“进了 pass.list 却没进 raw.fa”的元件：

```python
# 用法：改 L 为本机 <genome>.EDTA.raw/LTR/ 路径
pl  = [l.split("\t")[0] for l in open(L+"chr1.fa.mod.pass.list")
       if l.strip() and not l.startswith("#")]
ori = [l[1:].split("#")[0].strip() for l in open(L+"chr1.fa.mod.LTR.intact.raw.fa")
       if l.startswith(">")]
norm = lambda n: n.split(":")[0] + ":" + "..".join(
    str(x) for x in sorted(int(v) for v in n.split(":")[1].split("..")))
ori_n = {norm(n) for n in ori}
# 注意：raw.fa 已经过 TEsorter 过滤，差额里既有本报告的“静默丢弃”，
# 也有真实的 false-positive 过滤。区分方法：看 .ori 阶段是否已缺（本报告即属此类）。
print("pass.list 有、raw.fa 无：", [n for n in pl if norm(n) not in ori_n])
```
