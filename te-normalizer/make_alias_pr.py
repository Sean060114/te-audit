#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_alias_pr.py —— 从 aliases_issue542.csv 生成 EDTA #542 的别名 PR 素材

输入：
  data/aliases_issue542.csv          73 个被报告名称的比对结果（61 个缺失）
  data/TE_Sequence_Ontology.txt      EDTA 官方 SO 别名表（PR 目标文件）

输出（默认写到 pr_issue542/）：
  aliases_triage.tsv                 61 行完整分诊表：目标 SO、决策、理由、风险
  aliases_patch_core.diff            50 条高置信度（结构后缀 + LTR 亚家族）
  aliases_patch_full.diff            56 条（core + 6 条判据型）
  TE_Sequence_Ontology.patched.txt   full 版的完整结果文件，便于本地 diff 复核
  PR_BODY.md                         可直接粘贴的英文 PR 描述

设计约束（重要）：
  1. 只往**已有行的第三列追加别名**，绝不新增 SO 条目行。
     因此无需改动 EDTA/bin/gff2bed.pl（该文件只在新增 SO 条目时才需同步）。
  2. 输出必须逐字节保持原文件风格：制表符分隔、别名逗号分隔且逗号后无空格、LF 换行、UTF-8 无 BOM。
  3. 决策写死在本文件的 DECISIONS 字典里，任何人可复核与修改。
"""
import csv
import difflib
import os
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CSV_IN = os.path.join(DATA, "aliases_issue542.csv")
SO_IN = os.path.join(DATA, "TE_Sequence_Ontology.txt")
OUT_DIR = os.path.join(HERE, "pr_issue542")

# 结构后缀：剥掉这些末段后，父名应为表内已有别名
QUALIFIERS = {"mite", "nmite", "unknown", "na", "mixture", "nltr", "n/a"}

# ---- 疑难项的显式决策 ----------------------------------------------------
# key: 原始名
#   target : 目标 SO 标准名；None 表示不进补丁
#   tier   : core / full  —— core 进两份补丁，full 只进 full 版
#   conf   : 置信度
#   reason / risk : 写进分诊表与 PR 正文
DECISIONS = {
    "TIR/Novosib": dict(
        target="Novosib_TIR_transposon", tier="core", conf="高",
        reason="表中已有 Novosib_TIR_transposon 行且已有 DNA/Novosib 别名，补 TIR/Novosib 属同义写法",
        risk="极低；Novosib 为已描述的 TIR 家族，归属无争议",
    ),
    "nMITE": dict(
        target="MITE", tier="full", conf="中高",
        reason="nMITE = non-autonomous MITE，而 MITE 定义本身即非自主元件，SO:0000338 即其类名",
        risk="若报告方另指某一特定家族的非自主拷贝，此别名会丢一层信息；但作为名称归一化可接受",
    ),
    "nMITE/NA": dict(
        target="MITE", tier="full", conf="中高",
        reason="同上，nMITE 的结构后缀写法",
        risk="同 nMITE",
    ),
    "Unknown//ClassI": dict(
        target="retrotransposon", tier="full", conf="中",
        reason="双斜杠是中间字段为空的解析产物；表内 retrotransposon 行已有 Class_I 别名",
        risk="把字面双斜杠串写进官方表略显宽松；建议在 PR 中说明这是为兼容输出工具",
    ),
    "Unknown//ClassII": dict(
        target="DNA_transposon", tier="full", conf="中",
        reason="同上；表内 DNA_transposon 行已有 Class_II 别名",
        risk="同 Unknown//ClassI",
    ),
    "ARTEFACT": dict(
        target=None, tier="core", conf="—",
        reason="维护者 2026-03-02 在 #542 原话：\"ARTEFACT is ambiguous to me ... "
               "the current default is probably good for now\"——即维持 repeat_fragment 兜底，不新增别名",
        risk="若纳入等于替维护者做他明确表示暂不做的分类决定；映射结果虽相同，但不该由 PR 推动",
        en="`ARTEFACT` — per your reply in #542 (2026-03-02): *\"ARTEFACT is "
           "ambiguous to me, and I am unsure which TE classification fits it better than "
           "repeat_fragment, so the current default is probably good for now.\"* The existing "
           "`repeat_fragment` fallback is therefore left as-is.",
    ),
    "LTR/PIF_Harbinger": dict(
        target=None, tier="core", conf="—",
        reason="PIF_Harbinger 属 TIR 超家族，报成 LTR 前缀是报告方分类错误，不应作为合法别名固化",
        risk="若强行加入会合法化一条错误分类，后续统计会把 TIR 元件算进 LTR",
        en="`LTR/PIF_Harbinger` — `PIF_Harbinger` is a TIR superfamily; the `LTR/` "
           "prefix looks like a misclassification on the reporting side. Encoding it "
           "as an alias would legitimise a wrong call and skew LTR vs TIR counts.",
    ),
    "mixture": dict(
        target=None, tier="core", conf="—",
        reason="mixture 是 TEsorter 对嵌合元件的标记，不是分类名；固化后会把嵌合件统一压成兜底类，掩盖真实归属",
        risk="保持缺失则报告方仍能看到原始标记",
        en="`mixture` — a TEsorter marker for chimeric elements rather than a "
           "taxonomic name. Mapping it to `repeat_fragment` would silently flatten "
           "every chimeric call into the fallback class and hide its real affiliation.",
    ),
    "mixture/NA": dict(
        target=None, tier="core", conf="—",
        reason="同上",
        risk="同上",
        en="`mixture/NA` — same as `mixture`.",
    ),
    "nLTR": dict(
        target=None, tier="core", conf="—",
        reason="nLTR 有歧义：可读作 non-LTR，也可读作 nested LTR，二者指向完全不同的 SO 分支",
        risk="建议在 PR 中请报告方澄清后单独处理，避免猜错方向",
        en="`nLTR` — ambiguous between *non-LTR* and *nested LTR*, which sit in "
           "completely different branches of the ontology. Worth confirming with the "
           "reporter before assigning a target.",
    ),
    "nLTR/NA": dict(
        target=None, tier="core", conf="—",
        reason="同上",
        risk="同上",
        en="`nLTR/NA` — same ambiguity as `nLTR`.",
    ),
}


def parse_so(path):
    """返回 (原始行列表, 记录列表)。记录只收非注释非空行。"""
    with open(path, encoding="utf-8", newline="") as fh:
        raw = fh.readlines()
    recs = []
    for i, line in enumerate(raw):
        s = line.rstrip("\n")
        if not s.strip() or s.lstrip().startswith("#"):
            continue
        f = s.split("\t")
        if len(f) < 3:
            continue
        recs.append(dict(
            idx=i,
            name=f[0],
            so_id=f[1],
            aliases=[a.strip() for a in f[2].split(",") if a.strip()],
            extra=f[3:],
        ))
    return raw, recs


def build_index(recs):
    """alias(小写) -> 记录"""
    idx = {}
    for r in recs:
        for a in r["aliases"]:
            idx.setdefault(a.lower(), r)
    return idx


def strip_qualifiers(name):
    parts = [p for p in name.split("/") if p]
    while len(parts) > 1 and parts[-1].lower() in QUALIFIERS:
        parts.pop()
    return "/".join(parts)


def resolve(name, idx):
    """返回 (目标记录 or None, 分档, 依据说明)"""
    base = strip_qualifiers(name)
    hit = idx.get(base.lower())
    if hit:
        tier = "结构后缀" if base != name else "直命中"
        return hit, tier, f"剥离结构后缀后命中已有别名 {base}"
    parts = base.split("/")
    if len(parts) >= 2:
        sf = "/".join(parts[:2])
        hit = idx.get(sf.lower())
        if hit:
            return hit, "亚家族", f"按前两级 {sf} 归入所属超家族"
    return None, "需人工", ""


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    # 清理上一版命名残留，避免旧补丁与新分组混淆
    for stale in ("aliases_patch_core.diff", "aliases_patch_full.diff"):
        p = os.path.join(OUT_DIR, stale)
        if os.path.exists(p):
            os.remove(p)
    raw, recs = parse_so(SO_IN)
    idx = build_index(recs)

    rows = [r for r in csv.DictReader(open(CSV_IN, encoding="utf-8-sig"))
            if "缺失" in r["状态"]]
    if not rows:
        sys.exit("未在 CSV 中找到缺失行")

    # 逐条判定
    plan = []          # 分诊表行
    by_target = OrderedDict()   # SO 名 -> [(name, tier)]
    for r in rows:
        name = r["被报告的名称"].strip()
        who = r["报告人@日期"].strip()
        dec = DECISIONS.get(name)

        if dec is not None:
            tgt_name = dec["target"]
            rec = next((x for x in recs if x["name"] == tgt_name), None) if tgt_name else None
            if tgt_name and rec is None:
                sys.exit(f"决策里的目标 SO 名不存在: {tgt_name}")
            tier, conf = dec["tier"], dec["conf"]
            reason, risk = dec["reason"], dec["risk"]
            bucket = "疑难-显式决策"
        else:
            rec, bucket, why = resolve(name, idx)
            if rec is None:
                sys.exit(f"未判定且无法自动映射: {name}")
            tgt_name = rec["name"]
            tier = "core"
            conf = "高"
            reason = why
            risk = "低；该父名已在表内，仅补结构后缀写法"
            if bucket == "亚家族":
                reason += "（EDTA 表内已有 LTR/Gypsy-Cigr、LINE/L1-Tx1 等亚家族别名先例）"
            # clade 级名称单独成组：维护者 2026-03-02 已表态「不是 RepeatMasker 库命名惯例」，
            # 故不并入主补丁，单独出一份备选补丁供讨论。
            if name.startswith(("LTR/Copia/", "LTR/Gypsy/")):
                tier = "clade"
                risk = "维护者已表示 clade 级命名不属于 RepeatMasker 库惯例，可能被判为超出范围"

        if rec is not None:
            existing = {a.lower() for a in rec["aliases"]}
            if name.lower() in existing:
                sys.exit(f"{name} 已存在于 {rec['name']}，分诊逻辑有误")
            by_target.setdefault(rec["name"], []).append((name, tier))

        plan.append(dict(
            name=name, who=who, bucket=bucket,
            target=tgt_name or "", so_id=(rec["so_id"] if rec else ""),
            decision="纳入" if rec is not None else "不纳入",
            tier=tier, conf=conf, reason=reason, risk=risk,
            en=(dec.get("en", "") if dec else ""),
        ))

    # ---- 生成打补丁后的文件 ----
    def emit(tier_filter):
        out = list(raw)
        n_alias = 0
        for r in recs:
            adds = [n for n, t in by_target.get(r["name"], []) if t in tier_filter]
            if not adds:
                continue
            existing = {a.lower() for a in r["aliases"]}
            fresh = [a for a in adds if a.lower() not in existing]
            if not fresh:
                continue
            fields = [r["name"], r["so_id"], ",".join(r["aliases"] + fresh)] + r["extra"]
            out[r["idx"]] = "\t".join(fields) + "\n"
            n_alias += len(fresh)
        return out, n_alias

    main, n_main = emit({"core"})
    plus, n_plus = emit({"core", "full"})
    clade, n_clade = emit({"clade"})
    allp, n_all = emit({"core", "full", "clade"})

    def write_patch(lines, tag, n_add, note):
        diff = difflib.unified_diff(
            raw, lines,
            fromfile="a/bin/TE_Sequence_Ontology.txt",
            tofile="b/bin/TE_Sequence_Ontology.txt",
            n=1,
        )
        body = "".join(diff)
        p = os.path.join(OUT_DIR, f"aliases_patch_{tag}.diff")
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"# {note}：新增 {n_add} 个别名\n")
            fh.write("# 仅改动已有行的第三列，不新增任何 SO 条目行\n")
            fh.write(body)
        return p, body

    p_main, d_main = write_patch(main, "main", n_main, "主补丁（建议提交）")
    p_plus, d_plus = write_patch(plus, "main_plus", n_plus, "主补丁 + 判据型")
    p_clade, d_clade = write_patch(clade, "clade", n_clade, "备选：clade 级别名（不并入主补丁）")

    for fname, lines in (("TE_Sequence_Ontology.patched.txt", plus),
                         ("TE_Sequence_Ontology.with_clade.txt", allp)):
        with open(os.path.join(OUT_DIR, fname), "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(lines)

    # ---- 分诊表 ----
    tier_label = {"core": "主补丁", "full": "判据型", "clade": "仅备选(clade)"}
    order = {"core": 0, "full": 1, "clade": 2}
    plan.sort(key=lambda x: (order.get(x["tier"], 0) if x["decision"] == "纳入" else 3,
                             x["bucket"], x["target"], x["name"]))
    tri = os.path.join(OUT_DIR, "aliases_triage.tsv")
    with open(tri, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["序号", "被报告名称", "报告人", "分档", "目标SO标准名",
                            "SO ID", "决策", "补丁分组", "置信度", "理由", "风险"]) + "\n")
        for i, p in enumerate(plan, 1):
            grp = tier_label.get(p["tier"], "—") if p["decision"] == "纳入" else "—"
            fh.write("\t".join([str(i), p["name"], p["who"], p["bucket"], p["target"],
                                p["so_id"], p["decision"], grp, p["conf"],
                                p["reason"], p["risk"]]) + "\n")

    # ---- 校验：打补丁后每个"纳入"的名称都必须能命中 ----
    for fname, wanted in (("TE_Sequence_Ontology.patched.txt", {"core", "full"}),
                          ("TE_Sequence_Ontology.with_clade.txt", {"core", "full", "clade"})):
        _, recs2 = parse_so(os.path.join(OUT_DIR, fname))
        idx2 = build_index(recs2)
        miss = [p["name"] for p in plan
                if p["decision"] == "纳入" and p["tier"] in wanted
                and p["name"].lower() not in idx2]
        if miss:
            sys.exit(f"{fname} 校验失败，以下名称仍未命中: {miss}")
        assert len(recs2) == len(recs), "补丁改变了记录数（说明新增了行），违反设计约束"

    # ---- PR 正文 ----
    n_in = sum(1 for p in plan if p["decision"] == "纳入")
    n_out = sum(1 for p in plan if p["decision"] == "不纳入")
    n_criteria = sum(1 for p in plan if p["decision"] == "纳入" and p["tier"] == "full")
    rows_main = len([r for r in recs
                     if any(t == "core" for _, t in by_target.get(r["name"], []))])
    pr = os.path.join(OUT_DIR, "PR_BODY.md")
    with open(pr, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"""# Add missing TE name aliases reported in #542

Follow-up to #542. Of the 73 names reported in that thread, 12 were already
resolvable in `TE_Sequence_Ontology.txt`. This PR adds aliases for {n_main} of the
remaining 61 — those that are structural variants of superfamilies already present in
the table, plus `TIR/Novosib`, whose target row already exists.

## What changed

- **Only the third column (alias list) of {rows_main} existing rows is extended.**
  No new SO entry line is added, and no SO name or SO ID is modified.
- Therefore `bin/gff2bed.pl` needs **no** corresponding update: as noted in the
  header of `TE_Sequence_Ontology.txt`, that file only has to change when a *new*
  SO entry line is introduced.
- Total: {n_main} aliases.

## Deliberately kept out: clade-level names ({n_clade})

In #542 you noted that clade-level classifications are *"not the convention of
RepeatMasker library naming"* and suggested stripping them from the library before
feeding it to EDTA2. This PR therefore leaves them untouched: names that carry a clade
segment (`LTR/Copia/Ale`, `LTR/Gypsy/Athila`, ...) continue to behave exactly as they do
today.

If you would rather have the table discard the clade segment automatically, that change
affects {n_clade} names across two rows (`Copia_LTR_retrotransposon`,
`Gypsy_LTR_retrotransposon`). I can open it as a separate PR so the two decisions stay
independent.

## Aliases added, grouped by target term

| Target SO name | SO ID | Aliases added |
|---|---|---|
""")
        for r in recs:
            adds = [n for n, t in by_target.get(r["name"], []) if t == "core"]
            if adds:
                fh.write(f"| `{r['name']}` | {r['so_id']} | {len(adds)}: "
                         + ", ".join(f"`{a}`" for a in adds) + " |\n")
        fh.write(f"""
## Deliberately **not** included ({n_out} names)

{chr(10).join('- ' + (p['en'] or p['reason']) for p in plan if p['decision'] == '不纳入')}

`TIR/Novosib` is included because a `Novosib_TIR_transposon` row already exists; only
the `TIR/` spelling was missing.

Note that `LTR/mixture` is included while the bare `mixture` is not: `LTR/mixture`
still names the class (an LTR retrotransposon whose clade is indeterminate), so what
remains after stripping is well defined, whereas `mixture` on its own carries no class
information to fall back on.

## Not included here: {n_criteria} judgement calls

`nMITE`, `nMITE/NA`, `Unknown//ClassI` and `Unknown//ClassII` have defensible targets,
but assigning them rests on a judgement (e.g. whether the `n` prefix marks a
non-autonomous copy, or which class an unspecified element belongs to) rather than on a
mechanical strip. They are left out of this PR; I can add them, or explain the reasoning
for each, if you find it useful.

## Verification

- Every added alias was checked to be absent from its target row beforehand
  (case-insensitive), so no duplicate is introduced.
- After patching, all {n_main} names resolve to an SO term via the alias index, and the
  record count is unchanged ({len(recs)}) — i.e. no SO entry was added.
- Tab separators, comma-without-space alias formatting, LF line endings and the absence
  of a BOM match the original file outside the edited cells.

A per-name breakdown of all 61 reported names (what resolves, what does not, and why) is
available on request.
""")

    # ---- 终端汇报 ----
    print(f"缺失行 {len(rows)} 条：纳入 {n_in}，不纳入 {n_out}")
    print(f"  主补丁 core        +{n_main}  （改 {rows_main} 行）")
    print(f"  主补丁+判据型      +{n_plus}")
    print(f"  备选 clade         +{n_clade}  （独立补丁）")
    print(f"  全量               +{n_all}")
    print(f"记录总数保持不变（{len(recs)}），未新增 SO 条目")
    print()
    for p in (tri, p_main, p_plus, p_clade, pr):
        print(f"  {os.path.relpath(p, HERE)}")
    print()
    print("=== 主补丁预览 ===")
    print(d_main)


if __name__ == "__main__":
    main()
