# 提交清单：#542 别名 PR

逐步照抄即可。所有命令均已在本机实测（git 2.55 / Windows），非推测。

> **已定策略（2026-11-16）**：本次只提主补丁（21 个别名）。clade 那 30 条不进这个 PR，
> 待维护者回应后再单独开一个 PR（见第 9 节）。理由是两者独立提交，clade 即使被拒也不影响主补丁。
>
> **本次实操参数**：fork = `https://github.com/newbiesean0114/EDTA`（已核实为 `oushujun/EDTA` 的 fork，
> 默认分支 master）；本地克隆目录 = `D:\github\EDTA`；分支名 `add-issue542-aliases`。

---

## 0. 基准：补丁对着哪个版本生成的

| 项 | 值 |
|---|---|
| 仓库 | `oushujun/EDTA`，分支 `master` |
| 目标文件 | `bin/TE_Sequence_Ontology.txt` |
| 基准文件 MD5 | `366ebd57098a2a875e5907df65c0ea26` |
| 基准文件字节数 | 15,738 |

已实测：本地 `data/TE_Sequence_Ontology.txt` 与上游 `master` 的 `bin/TE_Sequence_Ontology.txt`
**逐字节一致**（同 MD5）。补丁也已在上游原始文件上实跑应用成功。

> 如果提交时该文件已变动，先重跑 `python make_alias_pr.py` 再提交，不要硬打。

---

## 1. 先改一个配置（不做这步大概率失败）

`bin/TE_Sequence_Ontology.txt` 是 **LF** 行尾、无 BOM。若 `core.autocrlf` 生效，
`git clone` 会把它转成 CRLF，补丁将无法应用或悄悄改坏整个文件的行尾。

```bash
git config --global core.autocrlf false
```

本机验证时出现过 `warning: LF will be replaced by CRLF`，说明确实存在该风险，务必先设。

---

## 2. Fork 并克隆

网页点 Fork，然后：

```bash
git clone https://github.com/<你的用户名>/EDTA.git
cd EDTA
git remote add upstream https://github.com/oushujun/EDTA.git

# 确认基准未变（应输出 366ebd57098a2a875e5907df65c0ea26）
md5sum bin/TE_Sequence_Ontology.txt
```

---

## 3. 建分支

```bash
git checkout -b add-issue542-aliases
```

---

## 4. 应用主补丁

把 `<素材目录>` 换成你的 `pr_issue542` 实际路径：

```bash
git apply -p1 <素材目录>/aliases_patch_main.diff
```

预期输出：无（静默成功）。随后：

```bash
git diff --stat
#  应为：1 file changed, 11 insertions(+), 11 deletions(-)
```

---

## 5. 三项自查（都要看到预期值）

**a) 记录数与别名数**

```bash
python -c "
rec=al=0
for l in open('bin/TE_Sequence_Ontology.txt',encoding='utf-8'):
    if l.startswith('#') or not l.strip(): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<3: continue
    rec+=1; al+=len([a for a in f[2].split(',') if a.strip()])
print('记录',rec,'别名',al)
"
#  应为：记录 98 别名 797
```

**b) 只改了第三列**

```bash
diff <(git diff -U0 | grep '^-' | grep -v '^---' | cut -c2- | cut -f1,2 | sort) \
     <(git diff -U0 | grep '^+' | grep -v '^+++' | cut -c2- | cut -f1,2 | sort) \
  && echo "OK: 前两列完全未变"
#  应输出：OK: 前两列完全未变
```

> 这条已做反向测试：故意改坏第一列时它会报差异，因此不是恒真。

**c) 行尾与 BOM 未被破坏**

```bash
python -c "
d=open('bin/TE_Sequence_Ontology.txt','rb').read()
print('CRLF',d.count(b'\r\n'),'LF',d.count(b'\n'),'BOM',d[:3]==b'\xef\xbb\xbf')
"
#  应为：CRLF 0 LF 156 BOM False
```

---

## 6. 提交

```bash
git add bin/TE_Sequence_Ontology.txt
git commit -F- <<'MSG'
Add missing TE name aliases reported in #542

Extend the alias column (third field) of 11 existing rows so that names
reported in #542 resolve on input. No SO entry is added and no SO name or
ID is modified, so bin/gff2bed.pl needs no corresponding update.

Alias count 776 -> 797; record count unchanged at 98.
MSG
```

---

## 7. 推送并开 PR

```bash
git push -u origin add-issue542-aliases
```

然后在网页上点 `Compare & pull request`，标题填：

```
Add missing TE name aliases reported in #542
```

正文直接粘贴 `PR_BODY.md` 全文。

已安装 `gh` 的话可省掉网页操作：

```bash
gh pr create --repo oushujun/EDTA --base master \
  --head <你的用户名>:add-issue542-aliases \
  --title "Add missing TE name aliases reported in #542" \
  --body-file <素材目录>/PR_BODY.md
```

---

## 8. 三件不要做的事

1. **不要对文件做格式化或全文件重排。** 会把 11 行 diff 变成上百行，维护者无法 review。
2. **不要动 `bin/gff2bed.pl`。** 该文件仅在新**增** SO 条目行时才需同步；本次只追加别名。
   这点已在 PR 正文中主动说明，避免维护者先要求补改再打回。
3. **不要顺手把 `aliases_patch_clade.diff` 一起提交。** 维护者曾明确表示 clade 级命名
   不是 RepeatMasker 库的惯例；未经他同意塞进同一个 PR 有被直接关闭的风险。

---

## 9. 若维护者同意收 clade 那 30 条

在已应用主补丁的分支上追加一个独立 commit：

```bash
git apply -p1 <素材目录>/aliases_patch_clade.diff
#  别名 797 -> 827，记录仍为 98
git commit -am "Also accept clade-level LTR aliases, stripped to subfamily
(currently discarded by users, see #542)"
```

这样他能只合第一个 commit、拒第二个。

---

## 10. 回滚

```bash
git checkout -- bin/TE_Sequence_Ontology.txt     # 丢弃未提交的改动
# 或整支作废：
git checkout master && git branch -D add-issue542-aliases
```

---

## 附：三份补丁的实测数据

| 补丁 | 别名数 | 别名总数 | 记录数 | 改动行 | 新增命中（本工具） |
|---|---|---|---|---|---|
| `aliases_patch_main.diff` | 21 | 776 → 797 | 98 | 11 | +5 |
| `aliases_patch_main_plus.diff` | 25 | 776 → 801 | 98 | 11 | +9 |
| `aliases_patch_clade.diff` | 30 | 776 → 806 | 98 | 2 | +0 |

"新增命中"指 `normalize.py` 打补丁前已能通过兜底路径猜中的数量之外的增量；
对只按表精确查名的消费者（EDTA 自身、TEsorter），增量等于别名数本身。
