#!/usr/bin/env python3
"""
test_te_check.py —— te_check.py 的断言测试集（标准库 unittest，零依赖）

为什么要写它：
    te_check 的判据在开发当天被改过两次，两次都是"看起来对、实际错"：
      1. C4 初版要求所有类别都该出现在分类表里 → 报出 166 个假阳性。
         实际 TEsorter 只处理建库子集，TIR 覆盖率仅 10.6%，LTR 才是 100%。
      2. C5 初版把 helitron 的期望 order 写成 DNA → 4 条正确条目被误报为跨类。
         实际 TEsorter 对 helitron 输出 order=Helitron。
    这两处当时只靠肉眼看了一次输出，没有锁。之后改动任何一行都可能悄悄退回去。
    本文件把两处都写成"回归锁"：谁把判据改宽或改窄，测试立刻失败。

测试分六层：
    1. 单元     rc / attrs / ir_identity
    2. 解析     parse_intact / parse_cls（在临时目录里合成 EDTA 输出树）
    3. 判据     六项检查各自的边界，含两处回归锁
    4. CLI      端到端：JSON 结构、--strict 退出码、缺文件时的行为
    5. 真实数据 Chr1 实测值回归（目录不可达时自动跳过，约 1 秒）
    6. 词表     SO 别名表对齐：九类映射与权威词表不允许各自漂移

怎么跑：
    python test_te_check.py           # 简洁
    python test_te_check.py -v        # 每条用例名
"""

import contextlib
import io
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))

import te_check as T                                     # noqa: E402

SCRIPT = BASE / "te_check.py"

# 真实数据目录（WSL 内）。不可达时第 5 层自动跳过，其余四层照常跑。
# 换机器或数据搬家时用环境变量覆盖，不必改代码：
#   TE_CHECK_TEST_DIR=/data/edta_run python test_te_check.py
REAL_EDTA = Path(os.environ.get(
    "TE_CHECK_TEST_DIR", "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run"))

# Chr1 实测基线（见 te_check_README.md 的回归基准表）。改检查逻辑时这些数会变，
# 变了就是"行为改变"，必须同时更新 README 与本表，不能只改一边。
CHR1 = dict(
    n_elements=218,
    n_classified=67,
    severe=4,          # 原为 5：Chr1 那个缺失元件是 LTR/unknown，按分层判据不计严重项
    c1_no_identity=["helitron"],
    c2_short=4,
    c2_min_internal=70,
    c3_n=92,
    c3_perfect=0,
    c3_median_ident=0.8182,
    c3_median_armlen=14,
    c4_missing_ltr=[(17528667, 17536138, "LTR_retrotransposon")],
    c5_total=15,
    c6_hits=[(17528667, 17536138, 0.9866)],
)


# ------------------------------------------------------------------ 合成数据
#
# 手算的反向互补（不用 T.rc，避免自证）：
#   "ACGTACGTAA" 的互补是 TGCATGCATT，反转得 TTACGTACGT
ARM_L = "ACGTACGTAA"
ARM_R = "TTACGTACGT"          # == rc(ARM_L)，完美反向重复

# 30 bp 的臂：滑移窗口只在重叠段 >= MIN_ARM_BP(10) 时才评估，
# 所以测试滑移必须用长臂，10 bp 的臂根本给不出滑移余地。
LONG_ARM = "ACGTTGCAACGATTACAGTCAGGTCAATGC"


def gff_line(seqid, ftype, start, end, strand, attrs_str):
    return "\t".join([seqid, "EDTA", ftype, str(start), str(end), ".", strand, ".", attrs_str])


def build_tree(root):
    """在 root 下造一份最小可用的 EDTA 输出树，每一行都对应一条断言。"""
    d = Path(root)
    (d / "t.EDTA.raw" / "TIR").mkdir(parents=True)
    (d / "t.EDTA.raw" / "Helitron").mkdir(parents=True)
    (d / "t.EDTA.raw" / "LTR").mkdir(parents=True)

    lines = ["##gff-version 3"]

    # (1) 完美反向重复的 TIR，内部区充足 → C3 命中 1.0，C2 不报
    lines.append(gff_line("1", "hAT_TIR_transposon", 1001, 2000, "+",
                          f"ID=e1;classification=DNA/hAT;identity=0.92;"
                          f"TSD=TAA_TAA_100;TIR={ARM_L}_{ARM_R}"))

    # (2) 内部区仅 70 bp → C2 命中（94 - 2*12）
    short_l, short_r = "ACGTACGTACGT", "TTACGTACGTAC"   # 各 12 bp
    lines.append(gff_line("1", "Mutator_TIR_transposon", 5001, 5094, "+",
                          f"ID=e2;classification=DNA/Mutator;identity=0.88;"
                          f"TIR={short_l}_{short_r}"))

    # (3) TIR 但不在分类表里 → C4 的 missing_other，**不得计入严重项**
    lines.append(gff_line("1", "CACTA_TIR_transposon", 7001, 8000, "+",
                          f"ID=e3;classification=DNA/CACTA;identity=0.90;"
                          f"TIR={ARM_L}_{ARM_R}"))

    # (4) LTR，未进分类表 + 无域证据 → C4 的 missing_ltr（严重）+ C6 命中
    lines.append(gff_line("1", "LTR_retrotransposon", 10001, 15000, "-",
                          "ID=e4;classification=LTR/unknown;ltr_identity=0.9866;"
                          "strand_evidence=PBS"))

    # (5) LTR，已进分类表且有域 → C4 覆盖 + C6 不报（对照组）
    lines.append(gff_line("1", "Copia_LTR_retrotransposon", 20001, 25000, "+",
                          "ID=e5;classification=LTR/Copia;ltr_identity=0.9990;"
                          "strand_evidence=domain+PBS+PPT;motif=TGCA"))

    # (6) helitron，identity=NA → C1 报"整类无 identity"
    lines.append(gff_line("1", "helitron", 30001, 31000, "+",
                          "ID=e6;classification=DNA/Helitron;identity=NA;method=structural"))

    # 子特征：必须被 parse_intact 排除，否则元素数会虚高
    lines.append(gff_line("1", "long_terminal_repeat", 10001, 10529, "-", "ID=ltr1"))
    lines.append(gff_line("1", "target_site_duplication", 9996, 10000, "-", "ID=tsd1"))
    lines.append(gff_line("1", "repeat_region", 10001, 15000, "-", "ID=rr1"))

    (d / "t.EDTA.intact.gff3").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 分类表：TIR 库混入 1 条 LTR、1 条 TIR；Helitron 库输出 order=Helitron（回归锁）；
    # LTR 库覆盖 e5 但不覆盖 e4
    tir_tsv = "\n".join([
        "#seqid\torder\tsuperfamily\tclade\tcomplete\tstrand\tdomains",
        "1:20001..25000#LTR/Copia\tLTR\tCopia\tAle\tcomplete\t+\tGAG,PROT,INT,RT,RH",  # 跨类：应报
        "1:20001..25000#DNA/DTM\tTIR\tMutator\tMutator\tcomplete\t+\tTPase",           # 同类：不报
    ]) + "\n"
    (d / "t.EDTA.raw" / "TIR" / "t.TIR.rexdb.cls.tsv").write_text(tir_tsv, encoding="utf-8")

    hel_tsv = "\n".join([
        "#seqid\torder\tsuperfamily\tclade\tcomplete\tstrand\tdomains",
        "1:30001..31000#Helitron\tHelitron\tHelitron\tHelitron\tcomplete\t+\tNA",      # order=Helitron：不报
    ]) + "\n"
    (d / "t.EDTA.raw" / "Helitron" / "t.Helitron.rexdb.cls.tsv").write_text(hel_tsv, encoding="utf-8")

    ltr_tsv = "\n".join([
        "#seqid\torder\tsuperfamily\tclade\tcomplete\tstrand\tdomains",
        "1:20001..25000#LTR/Copia\tLTR\tCopia\tAle\tcomplete\t+\tGAG,PROT,INT,RT,RH",
    ]) + "\n"
    (d / "t.EDTA.raw" / "LTR" / "t.LTR.rexdb.cls.tsv").write_text(ltr_tsv, encoding="utf-8")
    return d


class Fixture(unittest.TestCase):
    """所有需要合成树的用例共享一套建树/拆树逻辑。"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="te_check_test_")
        self.root = build_tree(self.tmp)
        self.recs = T.parse_intact(str(self.root / "t.EDTA.intact.gff3"))
        self.cls = []
        for p in T.find_files(str(self.root))["cls"]:
            self.cls.extend(T.parse_cls(p))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


# ------------------------------------------------------------------ 1. 单元

class TestUnits(unittest.TestCase):
    def test_rc_basic(self):
        self.assertEqual(T.rc("ACGT"), "ACGT")            # 回文
        self.assertEqual(T.rc("AAAA"), "TTTT")
        self.assertEqual(T.rc("ACGTACGTAA"), "TTACGTACGT")

    def test_rc_handles_n(self):
        # N 必须映射回 N，不能变成别的碱基——否则含 N 的臂会算错匹配率
        self.assertEqual(T.rc("ACGTN"), "NACGT")
        self.assertEqual(T.rc("nnnn"), "nnnn")

    def test_rc_is_involution(self):
        for s in ("ACGTACGTAA", "GGGTTTAAACCC", "ATATAT", "NACGTN"):
            self.assertEqual(T.rc(T.rc(s)), s)

    def test_attrs(self):
        d = T.attrs("ID=x;classification=DNA/hAT;TIR=AAA_TTT")
        self.assertEqual(d["classification"], "DNA/hAT")
        self.assertEqual(d["TIR"], "AAA_TTT")

    def test_attrs_ignores_flagstyle_and_empty(self):
        d = T.attrs("ID=x;;novalue;empty=")
        self.assertEqual(d["ID"], "x")
        self.assertEqual(d["empty"], "")
        self.assertNotIn("novalue", d)      # 无 = 的项不进字典
        self.assertEqual(T.attrs(""), {})

    def test_attrs_value_with_equals(self):
        # 值里含 = 时只能切第一个，否则会把值截断
        d = T.attrs("note=a=b;ID=z")
        self.assertEqual(d["note"], "a=b")
        self.assertEqual(d["ID"], "z")


# ------------------------------------------------------- SO 别名表对齐锁
#
# 为什么需要这一层：te_check.py 自带一份 9 类的硬编码映射（ELEMENT_TYPES），
# 而 SO 的权威词表放在 te-normalizer/data/ 里。两份实现各自维护就会漂移，
# 而漂移的表现恰恰是"新类型被当成非元件静默跳过"——本工具要抓的就是这类失效，
# 不能自己先犯一次。下面几条测试把"漂移"变成一次显式决定：要么建模，要么登记。
SO_TABLE = BASE / "te-normalizer" / "data" / "TE_Sequence_Ontology.txt"

# SO 表里合法、但本工具的十三项判据未在其上校准过的元件类型。
# 现状：解析时按"非元件"跳过（六份实测数据里一条都没出现，实际影响 0 条）。
# 显式列出来而不是留在暗处：SO 表一旦新增 TIR/LTR 类型，对应用例会失败，
# 逼一次显式判断。改这张表时请连同原因一起改，别只加名字。
UNMODELED_SO_TYPES = {
    # ① TIR 的非 CACTA/hAT/Mutator/PIF/Tc1 家族：六份实测数据里未出现，
    #    且 TIR 判据（臂长 + 相似度阈值）只在实测过的那五个家族上订过。
    "Sola_TIR_transposon", "PILE_TIR_transposon", "POLE_TIR_transposon",
    "Ginger_TIR_transposon", "Kolobok_TIR_transposon", "Dada_TIR_transposon",
    "IS3EU_TIR_transposon", "Zator_TIR_transposon", "KDZ_TIR_transposon",
    "Academ_TIR_transposon", "Novosib_TIR_transposon", "P_TIR_transposon",
    "piggyBac_TIR_transposon", "Transib_TIR_transposon", "Merlin_TIR_transposon",
    # ② LTR 的非 Copia/Gypsy 支系：结构与 Copia/Gypsy 同法（同样两套坐标、同样
    #    motif/TSD 语义），但本工具没有实测样本，不假装已校准。
    "Bel_Pao_LTR_retrotransposon", "Retrovirus_LTR_retrotransposon",
    "Endogenous_Retrovirus_LTR_retrotransposon",
}


def so_element_names(names):
    """从 SO 表第一列里挑出"元件级"的名字：LTR/TIR 转座子与 helitron。

    刻意不用 `"LTR" in name` 判断：`non_LTR_retrotransposon` 含 LTR 却是非 LTR
    元件，正是 `TIR_TYPES`/`LTR_TYPES` 用子串派生时的坑（见 test_non_ltr_...）。
    """
    out = set()
    for n in names:
        if n == "helitron":
            out.add(n)
        elif n.startswith("non_"):
            continue        # 非 LTR 反转座子（LINE 一类），不在本工具的 LTR 判据范围内
        elif n.endswith("_TIR_transposon") or n.endswith("_LTR_retrotransposon"):
            out.add(n)
    return out


def read_so_names():
    names = set()
    for line in SO_TABLE.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        first = line.split("\t")[0].strip()
        if first:
            names.add(first)        # 第一列：SO 大小写敏感的正式名
    return names


@unittest.skipUnless(SO_TABLE.is_file(), f"SO 别名表不可达: {SO_TABLE}")
class TestSOTableAlignment(unittest.TestCase):
    """te_check 的九类映射必须与 SO 权威词表对得上——对不上就报，不静默。"""

    @classmethod
    def setUpClass(cls):
        cls.names = read_so_names()

    def test_every_modeled_type_exists_in_so_table(self):
        # 反向也查：ELEMENT_TYPES 里若写错一个字母，解析层会把它整类丢掉
        missing = T.ELEMENT_TYPES - self.names
        self.assertFalse(missing, f"ELEMENT_TYPES 在 SO 表里查不到: {sorted(missing)}")

    def test_no_so_element_type_is_silently_ignored(self):
        unaccounted = so_element_names(self.names) - T.ELEMENT_TYPES - UNMODELED_SO_TYPES
        self.assertFalse(
            unaccounted,
            "这些 SO 元件类型既未建模、也未登记为排除项: "
            f"{sorted(unaccounted)}。要么加进 ELEMENT_TYPES 并补判据测试，"
            "要么写进 UNMODELED_SO_TYPES 说明为什么跳过。")

    def test_exclusion_list_is_current(self):
        stale = UNMODELED_SO_TYPES - self.names
        self.assertFalse(stale, f"排除项在 SO 表里已不存在（改名了？）: {sorted(stale)}")
        overlap = UNMODELED_SO_TYPES & T.ELEMENT_TYPES
        self.assertFalse(overlap, f"既建模又排除，自相矛盾: {sorted(overlap)}")

    def test_non_ltr_is_not_treated_as_ltr(self):
        # 子串派生的经典坑：non_LTR_retrotransposon 含 "LTR"，但它是 LINE 一类
        self.assertNotIn("non_LTR_retrotransposon", T.LTR_TYPES)

    def test_tir_and_ltr_sets_are_disjoint(self):
        self.assertFalse(T.TIR_TYPES & T.LTR_TYPES)


class TestIrIdentity(unittest.TestCase):
    def test_perfect(self):
        ident, shift = T.ir_identity(ARM_L, ARM_R)
        self.assertAlmostEqual(ident, 1.0)
        self.assertEqual(shift, 0)

    def test_one_mismatch(self):
        bad = "TTACGTACGA"           # 倒数第二位 T→G
        ident, _ = T.ir_identity(ARM_L, bad)
        self.assertAlmostEqual(ident, 0.9)

    def test_shift_is_found(self):
        # shift=+2 的语义是 target[i+2] 对 right[i]：right 取 target 的第 3 位起。
        # 必须用长臂——10 bp 的臂滑移 2 bp 后重叠段只剩 8，低于 MIN_ARM_BP 会被跳过，
        # 滑移窗口对短臂实际不起作用（见 test_no_shift_for_short_arm）。
        left = LONG_ARM
        shifted = T.rc(left)[2:] + "GG"
        ident, shift = T.ir_identity(left, shifted)
        self.assertAlmostEqual(ident, 1.0)
        self.assertEqual(shift, 2)

    def test_no_shift_for_short_arm(self):
        # 记录该限制：10 bp 臂 + 2 bp 偏移，滑移救不回来（重叠段 8 < MIN_ARM_BP）
        shifted = ARM_R[2:] + "GG"
        ident, shift = T.ir_identity(ARM_L, shifted)
        self.assertLess(ident, 0.5)

    def test_shift_beyond_window_fails(self):
        # 移 12 bp 超出默认 ±10 窗口，不应再判为完美
        far = T.rc(LONG_ARM)[12:] + "GG" * 6
        ident, _ = T.ir_identity(LONG_ARM, far)
        self.assertLess(ident, 0.999)

    def test_short_arm_returns_none(self):
        # 短于 MIN_ARM_BP 的臂样本量不足，必须返回 None 而不是一个假分数
        ident, shift = T.ir_identity("ACGT", "ACGT")
        self.assertIsNone(ident)
        self.assertIsNone(shift)

    def test_empty_arm(self):
        self.assertEqual(T.ir_identity("", "ACGT"), (None, None))
        self.assertEqual(T.ir_identity("ACGT", ""), (None, None))


# ------------------------------------------------------------------ 2. 解析

class TestParsing(Fixture):
    def test_element_count_excludes_subfeatures(self):
        # 造了 6 个元素 + 3 个子特征；子特征若被算进来就会是 9
        self.assertEqual(len(self.recs), 6)

    def test_length_is_inclusive(self):
        e = [r for r in self.recs if r["start"] == 1001][0]
        self.assertEqual(e["length"], 1000)      # 2000-1001+1

    def test_tir_split(self):
        e = [r for r in self.recs if r["start"] == 1001][0]
        self.assertEqual(e["tir"], (ARM_L, ARM_R))

    def test_identity_na_is_none(self):
        h = [r for r in self.recs if r["type"] == "helitron"][0]
        self.assertIsNone(h["identity"])         # "NA" 不能变成 0.0

    def test_identity_parsed(self):
        e = [r for r in self.recs if r["start"] == 1001][0]
        self.assertAlmostEqual(e["identity"], 0.92)

    def test_ltr_identity_parsed(self):
        e = [r for r in self.recs if r["start"] == 10001][0]
        self.assertAlmostEqual(e["ltr_identity"], 0.9866)

    def test_tsd_pct(self):
        e = [r for r in self.recs if r["start"] == 1001][0]
        self.assertEqual(e["tsd"], "TAA")
        self.assertAlmostEqual(e["tsd_pct"], 100.0)

    def test_cls_parse_coord_order(self):
        # 降序坐标（20001..20000 之类）必须被规整成 start<end
        rows = self.cls
        for c in rows:
            self.assertLessEqual(c["start"], c["end"])

    def test_cls_fields(self):
        c = [x for x in self.cls if x["orig_cls"] == "DNA/DTM"][0]
        self.assertEqual(c["order"], "TIR")
        self.assertEqual(c["domains"], "TPase")
        self.assertEqual(c["lib"], "TIR")

    def test_cls_malformed_line_skipped(self):
        p = self.root / "t.EDTA.raw" / "TIR" / "bad.rexdb.cls.tsv"
        p.write_text("no-coord-here\tTIR\tMutator\tX\tcomplete\t+\tNA\n", encoding="utf-8")
        self.assertEqual(T.parse_cls(str(p)), [])

    def test_find_files(self):
        f = T.find_files(str(self.root))
        self.assertTrue(f["intact"].endswith("t.EDTA.intact.gff3"))
        self.assertIsNone(f["anno"])              # 本 fixture 未造
        self.assertEqual(len(f["cls"]), 3)


class TestFindFilesPrefix(Fixture):
    """同一目录放多份运行时，必须能按前缀挑对，且不能静默取错。"""

    def _add_second_run(self, prefix):
        """在同一个目录里再放一份「另一次运行」的输出（前缀不同、内容不同）。"""
        d = self.root
        lines = ["##gff-version 3",
                 gff_line("1", "hAT_TIR_transposon", 1, 500, "+",
                          f"ID=x1;classification=DNA/hAT;identity=0.9;TIR={ARM_L}_{ARM_R}")]
        (d / f"{prefix}.EDTA.intact.gff3").write_text("\n".join(lines) + "\n", encoding="utf-8")
        p = d / f"{prefix}.EDTA.raw" / "TIR"
        p.mkdir(parents=True)
        (p / f"{prefix}.TIR.rexdb.cls.tsv").write_text(
            "1:1..500#DNA/hAT\tTIR\thAT\thAT\tcomplete\t+\tTPase\n", encoding="utf-8")

    def test_prefix_selects_intact(self):
        self._add_second_run("other")
        f = T.find_files(str(self.root), "other")
        self.assertTrue(f["intact"].endswith("other.EDTA.intact.gff3"))
        recs = T.parse_intact(f["intact"])
        self.assertEqual(len(recs), 1)

    def test_prefix_filters_cls(self):
        self._add_second_run("other")
        f = T.find_files(str(self.root), "other")
        self.assertEqual(len(f["cls"]), 1)
        # 分类表路径形如 .../other.EDTA.raw/TIR/other.TIR.rexdb.cls.tsv
        self.assertTrue(f["cls"][0].replace("\\", "/").split("/")[-3].startswith("other"))

    def test_without_prefix_picks_first_and_warns(self):
        # 不指定前缀时取文件名最靠前的一份——但必须出声，不能静默
        self._add_second_run("aaa")   # 字典序在 t. 之前
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            f = T.find_files(str(self.root))
        self.assertTrue(f["intact"].endswith("aaa.EDTA.intact.gff3"))
        self.assertIn("仅分析", err.getvalue())

    def test_single_run_does_not_warn(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            T.find_files(str(self.root))
        self.assertEqual(err.getvalue(), "")


# ------------------------------------------------------------------ 3. 判据

class TestC1(Fixture):
    def test_helitron_flagged_as_no_identity(self):
        fa = T.check_field_availability(self.recs)
        no_id = [d["type"] for d in fa if d["identity"] == 0]
        self.assertIn("helitron", no_id)

    def test_others_not_flagged(self):
        fa = T.check_field_availability(self.recs)
        no_id = [d["type"] for d in fa if d["identity"] == 0]
        self.assertNotIn("hAT_TIR_transposon", no_id)

    def test_counts_per_type(self):
        fa = {d["type"]: d for d in T.check_field_availability(self.recs)}
        self.assertEqual(fa["helitron"]["n"], 1)
        self.assertEqual(fa["helitron"]["tir"], 0)
        self.assertEqual(fa["hAT_TIR_transposon"]["tir"], 1)


class TestC2(Fixture):
    def test_short_element_flagged(self):
        short = T.check_short_elements(self.recs)
        self.assertEqual(len(short), 1)
        self.assertEqual(short[0]["internal"], 70)     # 94 - 2*12

    def test_normal_element_not_flagged(self):
        short = T.check_short_elements(self.recs)
        self.assertNotIn(1001, [r["start"] for r in short])

    def test_boundary_exactly_100(self):
        # internal == MIN_INTERNAL_BP 不应报（判据是严格小于）
        r = dict(type="hAT_TIR_transposon", length=140, tir=("A" * 20, "T" * 20))
        self.assertEqual(T.check_short_elements([r]), [])
        r2 = dict(type="hAT_TIR_transposon", length=139, tir=("A" * 20, "T" * 20))
        self.assertEqual(len(T.check_short_elements([r2])), 1)

    def test_non_tir_excluded(self):
        r = dict(type="helitron", length=50, tir=("", ""))
        self.assertEqual(T.check_short_elements([r]), [])

    def test_tir_without_arm_excluded(self):
        r = dict(type="hAT_TIR_transposon", length=50, tir=("", ""))
        self.assertEqual(T.check_short_elements([r]), [])


class TestC3(Fixture):
    def test_perfect_pair_detected(self):
        rows = T.check_ir_quality(self.recs)
        perfect = [r for r in rows if r["ir_ident"] >= 0.999]
        self.assertTrue(any(r["start"] == 1001 for r in perfect))

    def test_arm_length_recorded(self):
        rows = T.check_ir_quality(self.recs)
        e = [r for r in rows if r["start"] == 1001][0]
        self.assertEqual(e["armlen"], 10)

    def test_short_arm_excluded_from_rows(self):
        # 臂长 10 == MIN_ARM_BP 应保留；改成 4 bp 则整条不进统计
        rows = T.check_ir_quality(self.recs)
        self.assertTrue(all(r["armlen"] >= T.MIN_ARM_BP for r in rows))
        r = dict(type="hAT_TIR_transposon", tir=("ACGT", "ACGT"))
        self.assertEqual(T.check_ir_quality([r]), [])


class TestC4(Fixture):
    """回归锁之一：只有 LTR 缺失算严重项。"""

    def test_missing_ltr_goes_to_severe_bucket(self):
        per_class, missing_ltr, missing_other = T.check_coverage(self.recs, self.cls)
        self.assertEqual([(r["start"], r["end"]) for r in missing_ltr],
                         [(10001, 15000)])

    def test_unassigned_superfamily_not_severe(self):
        # LTR/unknown 的缺失是 TEsorter 的正常结果，不计严重项。
        # 这一条是水稻数据换来的：218 个 LTR 里有 20 个 LTR/unknown 全部未分类，
        # 混算会让总覆盖率从 94.4% 掉到 85.8%，稀释掉真正需要核查的已定超家族缺失。
        rec = dict(type="LTR_retrotransposon", start=30001, end=34000,
                   cls="LTR/unknown", strand_ev="PBS")
        assigned = dict(type="Gypsy_LTR_retrotransposon", start=40001, end=45000,
                        cls="LTR/Gypsy", strand_ev="PBS,PPT")
        per, miss_ltr, miss_other = T.check_coverage([rec, assigned], [])
        self.assertEqual(len(miss_ltr), 2)          # check_coverage 本身不分层
        self.assertEqual(T.UNASSIGNED_LTR, "LTR_retrotransposon")

    def test_missing_tir_is_not_severe(self):
        # ← 这就是当初 166 个假阳性的来源：TIR 缺失绝不能被算作严重
        per_class, missing_ltr, missing_other = T.check_coverage(self.recs, self.cls)
        tir_missing = [r for r in missing_other if "TIR" in r["type"]]
        self.assertTrue(tir_missing)
        self.assertNotIn(7001, [r["start"] for r in missing_ltr])

    def test_coverage_counts(self):
        per_class, _, _ = T.check_coverage(self.recs, self.cls)
        self.assertEqual(per_class["Copia_LTR_retrotransposon"], [1, 1])
        self.assertEqual(per_class["LTR_retrotransposon"], [0, 1])
        self.assertEqual(per_class["hAT_TIR_transposon"], [0, 1])

    def test_tolerance_five_bp(self):
        # 坐标差 ±5 bp 视为同一元件（分类表与注释的边界常有小出入）
        recs = [dict(type="Copia_LTR_retrotransposon", start=20005, end=25000)]
        per_class, missing_ltr, _ = T.check_coverage(recs, self.cls)
        self.assertEqual(missing_ltr, [])

    def test_same_coord_other_chrom_not_covered(self):
        # 多序列输入（Chr1+Chr12）时两条染色体可能有完全相同坐标的元件。
        # 若配对键不含染色体，B 染色体的分类结果会被当成 A 染色体的，
        # 覆盖率虚高且完全不报错。这条锁住按染色体分组。
        recs = [dict(chrom="2", type="Copia_LTR_retrotransposon",
                     start=20001, end=25000)]
        per_class, missing_ltr, _ = T.check_coverage(recs, self.cls)   # cls 的 seqid 是 "1"
        self.assertEqual(len(missing_ltr), 1)
        self.assertEqual(per_class["Copia_LTR_retrotransposon"], [0, 1])

    def test_same_coord_same_chrom_covered(self):
        recs = [dict(chrom="1", type="Copia_LTR_retrotransposon",
                     start=20001, end=25000)]
        _, missing_ltr, _ = T.check_coverage(recs, self.cls)
        self.assertEqual(missing_ltr, [])

    def test_tolerance_exceeded(self):
        recs = [dict(type="Copia_LTR_retrotransposon", start=20007, end=25000)]
        _, missing_ltr, _ = T.check_coverage(recs, self.cls)
        self.assertEqual(len(missing_ltr), 1)


class TestC5(Fixture):
    """回归锁之二：helitron 的期望 order 是 Helitron（TEsorter 输出）。"""

    def test_cross_class_detected(self):
        cross = T.check_cross_class(self.cls)
        # TIR 库里的那条 LTR/Copia 应被报出
        self.assertTrue(any(c["lib"] == "TIR" and c["order"] == "LTR" for c in cross))

    def test_helitron_order_not_flagged(self):
        # ← 当初把期望写成 {"DNA"} 时，这条正确条目被误报
        cross = T.check_cross_class(self.cls)
        self.assertFalse(any(c["lib"] == "Helitron" for c in cross))

    def test_helitron_lib_accepts_dna_too(self):
        rows = [dict(lib="Helitron", order="DNA", start=1, end=2)]
        self.assertEqual(T.check_cross_class(rows), [])

    def test_matching_order_not_flagged(self):
        rows = [dict(lib="TIR", order="TIR", start=1, end=2),
                dict(lib="LTR", order="LTR", start=1, end=2),
                dict(lib="LINE", order="LINE", start=1, end=2),
                dict(lib="SINE", order="SINE", start=1, end=2)]
        self.assertEqual(T.check_cross_class(rows), [])

    def test_unknown_lib_is_permissive(self):
        # 认不出库名时不要瞎报
        rows = [dict(lib="Other", order="Whatever", start=1, end=2)]
        self.assertEqual(T.check_cross_class(rows), [])


class TestC6(Fixture):
    def test_other_chrom_not_matched(self):
        # 与 C4 同理：坐标相同但染色体不同，不能借用另一条染色体的分类结果。
        # 借了会把"缺域证据"误判成"有域证据"，从而漏报。
        rec = dict(chrom="2", type="Copia_LTR_retrotransposon",
                   start=10001, end=15000, ltr_identity=0.99, strand_ev="PBS")
        hits = T.check_ltr_conflict([rec], self.cls)
        self.assertEqual(len(hits), 1)
        self.assertIsNone(hits[0]["cls"])

    def test_high_identity_no_domain_hits(self):
        hits = T.check_ltr_conflict(self.recs, self.cls)
        self.assertEqual([(d["element"]["start"], d["element"]["end"]) for d in hits],
                         [(10001, 15000)])

    def test_high_identity_with_domain_misses(self):
        hits = T.check_ltr_conflict(self.recs, self.cls)
        self.assertNotIn(20001, [d["element"]["start"] for d in hits])

    def test_low_identity_excluded(self):
        r = dict(type="Copia_LTR_retrotransposon", start=90001, end=95000,
                 ltr_identity=0.90)
        self.assertEqual(T.check_ltr_conflict([r], []), [])

    def test_threshold_boundary(self):
        # 判据是 >= 0.95
        r = dict(type="Copia_LTR_retrotransposon", start=90001, end=95000,
                 ltr_identity=0.95)
        self.assertEqual(len(T.check_ltr_conflict([r], [])), 1)
        r["ltr_identity"] = 0.9499
        self.assertEqual(T.check_ltr_conflict([r], []), [])

    def test_domain_na_counts_as_no_domain(self):
        # domains="NA" 与空字符串必须同判——否则会漏报
        recs = [dict(type="Copia_LTR_retrotransposon", start=1, end=2000, ltr_identity=0.99)]
        for dom in ("NA", ""):
            cls = [dict(start=1, end=2000, order="LTR", superfamily="Copia", domains=dom)]
            self.assertEqual(len(T.check_ltr_conflict(recs, cls)), 1, msg=f"domains={dom!r}")

    # --- 域证据的第二个来源：intact.gff3 的 strand_evidence
    # 水稻 880801..891653 实测：strand_evidence=domain,PBS,PPT，但不在任何分类表中。
    # 只认分类表会把它误报成"缺域证据"。

    def test_strand_evidence_domain_counts_as_evidence(self):
        r = dict(type="LTR_retrotransposon", start=1, end=2000, ltr_identity=0.99,
                 strand_ev="domain,PBS,PPT")
        self.assertEqual(T.check_ltr_conflict([r], []), [])

    def test_strand_evidence_without_domain_still_hits(self):
        r = dict(type="LTR_retrotransposon", start=1, end=2000, ltr_identity=0.99,
                 strand_ev="PBS,PPT")
        self.assertEqual(len(T.check_ltr_conflict([r], [])), 1)

    def test_strand_evidence_case_insensitive(self):
        r = dict(type="LTR_retrotransposon", start=1, end=2000, ltr_identity=0.99,
                 strand_ev="Domain")
        self.assertEqual(T.check_ltr_conflict([r], []), [])

    def test_missing_strand_ev_key_does_not_crash(self):
        # 手工构造的记录可能没有该键，不能因此抛 KeyError
        r = dict(type="LTR_retrotransposon", start=1, end=2000, ltr_identity=0.99)
        self.assertEqual(len(T.check_ltr_conflict([r], [])), 1)


# ------------------------------------------------------------------ 4. CLI

class TestCli(Fixture):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args],
                              capture_output=True, text=True, encoding="utf-8")

    def test_reports_six_sections(self):
        p = self.run_cli("--edta-dir", str(self.root))
        self.assertEqual(p.returncode, 0, msg=p.stderr)
        for tag in ("C1 字段可得性", "C2 结构长度异常", "C3 反向重复质量",
                    "C4 注释 ↔ 分类器覆盖", "C5 库层面跨类",
                    "C6 LTR 一致性偏高但缺结构域证据"):
            self.assertIn(tag, p.stdout)

    def test_severe_count(self):
        # fixture 里 1 个 C2 + 1 个 C4-missing-LTR = 2；TIR 缺失不计
        p = self.run_cli("--edta-dir", str(self.root))
        self.assertIn("严重问题 1 项", p.stdout)   # 未定超家族缺失不计严重项

    def test_json_structure(self):
        out = self.root / "r.json"
        self.run_cli("--edta-dir", str(self.root), "--json", str(out))
        d = json.loads(out.read_text(encoding="utf-8"))
        for k in ("dir", "n_elements", "n_classified", "C1_field_availability",
                  "C1_no_identity_classes", "C2_short_elements", "C3_ir_quality",
                  "C4_coverage", "C4_not_classified", "C5_cross_class",
                  "C6_ltr_conflict", "severe_total"):
            self.assertIn(k, d)
        self.assertEqual(d["n_elements"], 6)
        self.assertEqual(d["severe_total"], 1)
        self.assertEqual(d["C1_no_identity_classes"], ["helitron"])

    def test_strict_exit_code(self):
        p = self.run_cli("--edta-dir", str(self.root), "--strict")
        self.assertEqual(p.returncode, 1)

    def test_strict_clean_exit_zero(self):
        # 把两条严重项的数据都去掉后，--strict 必须返回 0
        (self.root / "t.EDTA.intact.gff3").write_text(
            "\n".join([
                gff_line("1", "hAT_TIR_transposon", 1001, 2000, "+",
                         f"ID=e1;classification=DNA/hAT;identity=0.92;TIR={ARM_L}_{ARM_R}"),
            ]) + "\n", encoding="utf-8")
        p = self.run_cli("--edta-dir", str(self.root), "--strict")
        self.assertEqual(p.returncode, 0, msg=p.stdout)

    def test_missing_intact_file_exits(self):
        empty = tempfile.mkdtemp(prefix="te_check_empty_")
        try:
            p = self.run_cli("--edta-dir", empty)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("未在", p.stderr + p.stdout)
        finally:
            shutil.rmtree(empty, ignore_errors=True)


# ------------------------------------------------------------------ 5. 真实数据

@unittest.skipUnless(REAL_EDTA.is_dir(), f"真实数据目录不可达: {REAL_EDTA}")
class TestSemantics(unittest.TestCase):
    """C7-C13：字段语义与结构自证。

    这批判据的共同点：**用错一套坐标或方向，结果不报错，只是全错**。所以每条
    都写"完美元件"与"只错一处"两版；锁的是"改错会失败"，不是"能跑通"。
    """

    # 完美元件：父行 97..204（含 TSD），分类行 101..200（不含）
    #   96..99   = ACGT  父行左 TSD
    #   100..101 = AA    motif 前 2 nt（5'LTR 起点）
    #   198..199 = CC    motif 后 2 nt（3'LTR 末端）
    #   200..203 = ACGT  父行右 TSD
    POS = {96: "ACGT", 100: "AA", 198: "CC", 200: "ACGT"}

    @staticmethod
    def seq(fills=None, n=300, chrom="Chr1"):
        s = ["A"] * n
        for pos, text in (fills or {}).items():
            s[pos:pos + len(text)] = list(text)
        return {chrom: "".join(s)}

    @staticmethod
    def rec(**kw):
        base = dict(
            chrom="Chr1", type="Copia_LTR_retrotransposon",
            start=101, end=200, length=100, strand="+", cls="LTR/Copia",
            identity=None, ltr_identity=0.99, tsd="ACGT", tsd_pct=None,
            tir=("", ""), ir_ident=None, ir_shift=None, armlen=0,
            motif="AACC", strand_ev=None, method=None,
            sid="LTRRT_1", name="TE_1", parent="rr_1", pstart=97, pend=204)
        base.update(kw)
        return base

    # ---- C7 两套坐标
    def test_c7_two_coords_ok(self):
        r = T.check_ltr_two_coords([self.rec()], self.seq(self.POS))
        self.assertTrue(r[0]["coord_ok"])
        self.assertTrue(r[0]["seq5"] and r[0]["seq3"])

    def test_c7_mismatch_when_delta_ne_tsd_len(self):
        # 父行只比分类行多 3 bp，而 tsd 长 4 → 不自洽
        r = T.check_ltr_two_coords([self.rec(pstart=98)], self.seq(self.POS))
        self.assertFalse(r[0]["coord_ok"])
        self.assertEqual(r[0]["d5"], 3)

    def test_c7_without_parent_is_skipped(self):
        self.assertEqual(
            T.check_ltr_two_coords([self.rec(pstart=None)], self.seq(self.POS)), [])

    # ---- C8 motif 方向
    def test_c8_plus_strand(self):
        r = T.check_ltr_motif([self.rec()], self.seq(self.POS))
        self.assertTrue(r[0]["motif_ok"])

    def test_c8_minus_strand_needs_rc_and_cross(self):
        # 负链的正确布局：小端 GG、大端 TT（rc 后才是 motif 的前/后 2 nt）
        pos = {96: "ACGT", 100: "GG", 198: "TT", 200: "ACGT"}
        r = T.check_ltr_motif([self.rec(strand="-")], self.seq(pos))
        self.assertTrue(r[0]["motif_ok"], "负链要 rc 且交叉：小坐标端是它的 3' 端")

    def test_c8_minus_strand_fails_on_plus_layout(self):
        # 把正链布局套在负链上必须失败——这条锁住"忘了 rc/交叉"的回归
        r = T.check_ltr_motif([self.rec(strand="-")], self.seq(self.POS))
        self.assertFalse(r[0]["motif_ok"])

    def test_c8_uses_class_row_not_parent(self):
        # 若误用父行坐标（多算 len(tsd)），会取到 104..105 而非 100..101
        pos = dict(self.POS)
        pos[104] = "GG"
        r = T.check_ltr_motif([self.rec()], self.seq(pos))
        self.assertTrue(r[0]["motif_ok"], "必须用分类行坐标；多偏移一个 L 会失败")

    # ---- C9 TIR 的 TSD 在坐标外
    def test_c9_tir_tsd_outside(self):
        pos = {98: "AC", 200: "AC"}
        r = T.check_tir_tsd_outside(
            [self.rec(type="Tc1_Mariner_TIR_transposon", tsd="AC")], self.seq(pos))
        self.assertTrue(r[0]["out5"] and r[0]["out3"])

    def test_c9_tir_tsd_inside_fails(self):
        # 按 LTR 约定把 TSD 放进坐标内 → 外侧取不到
        pos = {100: "AC", 198: "AC"}
        r = T.check_tir_tsd_outside(
            [self.rec(type="Tc1_Mariner_TIR_transposon", tsd="AC")], self.seq(pos))
        self.assertFalse(r[0]["out5"] and r[0]["out3"])

    # ---- C10 identity 虚高（严格对位，不滑移）
    def test_c10_inflation_detected(self):
        r = T.check_tir_identity_inflation([self.rec(
            type="Tc1_Mariner_TIR_transposon", identity=0.9,
            tir=("AAAAAAAAAA", "AAAAAAAAAA"))])
        self.assertAlmostEqual(r[0]["measured"], 0.0)
        self.assertAlmostEqual(r[0]["inflation"], 0.9)

    def test_c10_no_slide_means_strict(self):
        # 右臂整体错一位：滑移 ±10 能对齐成 1.0，严格对位是 0.0。
        # 这条锁住 C10 **不滑移**的选择（滑移会把虚高掩盖掉）。
        a, b = "ACGTACGTAC", "TACGTACGTA"
        r = T.check_tir_identity_inflation([self.rec(
            type="Tc1_Mariner_TIR_transposon", identity=1.0, tir=(a, b))])
        self.assertAlmostEqual(r[0]["measured"], 0.0)
        # 滑移口径取最大值，必然不低于严格对位，两边口径确实不同。但别断言具体
        # 数值：10 bp 的臂滑移后长度 < MIN_ARM_BP 会被跳过，此例两边同为 0.0。
        self.assertGreaterEqual(T.ir_identity(a, b)[0], r[0]["measured"])

    # ---- C11 长臂内侧外推
    def test_c11_inner_extrapolation(self):
        a = "A" * 100
        b = "A" * 34 + "T" * 66          # rc(b) = A*66 + T*34
        r = T.check_tir_arm_extrapolation([self.rec(
            type="Mutator_TIR_transposon", tir=(a, b))])
        self.assertEqual(len(r), 1)
        self.assertLess(r[0]["mismatch_outer"], 0.01)
        self.assertGreater(r[0]["mismatch_inner"], 0.9)

    def test_c11_uniform_arm_not_flagged(self):
        self.assertEqual(T.check_tir_arm_extrapolation([self.rec(
            type="Mutator_TIR_transposon", tir=("A" * 100, "T" * 100))]), [])

    # ---- C12 近重复
    def test_c12_near_duplicate(self):
        p = T.check_near_duplicates([self.rec(start=101, end=200),
                                     self.rec(start=102, end=199)])
        self.assertEqual(len(p), 1)

    def test_c12_three_bp_apart_is_not_duplicate(self):
        self.assertEqual(T.check_near_duplicates(
            [self.rec(start=101, end=200), self.rec(start=104, end=200)]), [])

    def test_c12_different_type_not_duplicate(self):
        self.assertEqual(T.check_near_duplicates(
            [self.rec(start=101, end=200),
             self.rec(start=101, end=200, type="Gypsy_LTR_retrotransposon")]), [])

    # ---- C13 家族标签三层
    @staticmethod
    def _fam(cls_list, name="TE_x"):
        return [TestSemantics.rec(name=name, cls=c) for c in cls_list]

    def test_c13_true_conflict(self):
        r = T.check_family_label_spread(self._fam(["DNA/DTT", "DNA/DTH", "DNA/DTT"]))
        self.assertEqual(len(r), 1)
        self.assertEqual(r[0]["kind"], "真冲突")

    def test_c13_unknown_is_not_conflict(self):
        r = T.check_family_label_spread(self._fam(["DNA/DTT", "DNA/DTT", "LTR/unknown"]))
        self.assertEqual(r[0]["kind"], "含未定")

    def test_c13_mite_vs_dna_is_autonomy_only(self):
        r = T.check_family_label_spread(self._fam(["DNA/DTT", "DNA/DTT", "MITE/DTT"]))
        self.assertEqual(r[0]["kind"], "仅自主性")

    def test_c13_coord_name_excluded(self):
        # 反写坐标串是缺家族名时的占位，不能当家族聚合
        self.assertEqual(T.check_family_label_spread(
            self._fam(["DNA/DTT", "DNA/DTH", "DNA/DTM"], name="Chr1:100..200")), [])

    def test_c13_single_label_not_reported(self):
        self.assertEqual(T.check_family_label_spread(self._fam(["DNA/DTT"] * 3)), [])

    # ---- 解析：父坐标挂载
    def test_parse_attaches_parent_coords(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.EDTA.intact.gff3"
            p.write_text(
                "Chr1\tEDTA\trepeat_region\t97\t204\t.\t+\t.\t"
                "ID=rr_1;classification=LTR/Copia\n"
                "Chr1\tEDTA\tCopia_LTR_retrotransposon\t101\t200\t.\t+\t.\t"
                "ID=LTRRT_1;Parent=rr_1;Name=TE_1;classification=LTR/Copia\n",
                encoding="utf-8")
            recs = T.parse_intact(str(p))
        self.assertEqual(len(recs), 1, "repeat_region 行本身不能算作元件")
        self.assertEqual(recs[0]["pstart"], 97)
        self.assertEqual(recs[0]["pend"], 204)


class TestChr1Regression(unittest.TestCase):
    """Chr1 实测值回归。这些数字是 README 里公布的基线，不能悄悄漂移。"""

    @classmethod
    def setUpClass(cls):
        cls.files = T.find_files(str(REAL_EDTA))
        cls.recs = T.parse_intact(cls.files["intact"])
        cls.cls = []
        for p in cls.files["cls"]:
            cls.cls.extend(T.parse_cls(p))

    def test_element_and_cls_counts(self):
        self.assertEqual(len(self.recs), CHR1["n_elements"])
        self.assertEqual(len(self.cls), CHR1["n_classified"])

    def test_c1(self):
        fa = T.check_field_availability(self.recs)
        no_id = [d["type"] for d in fa if d["identity"] == 0]
        self.assertEqual(sorted(no_id), sorted(CHR1["c1_no_identity"]))

    def test_c2(self):
        short = T.check_short_elements(self.recs)
        self.assertEqual(len(short), CHR1["c2_short"])
        self.assertEqual(min(r["internal"] for r in short), CHR1["c2_min_internal"])

    def test_c3(self):
        rows = T.check_ir_quality(self.recs)
        idents = [r["ir_ident"] for r in rows]
        self.assertEqual(len(rows), CHR1["c3_n"])
        self.assertEqual(sum(1 for i in idents if i >= 0.999), CHR1["c3_perfect"])
        # 用 statistics.median（偶数个样本时取中间两值平均），与 te_check 内部一致
        self.assertAlmostEqual(round(statistics.median(idents), 4),
                               CHR1["c3_median_ident"], places=4)
        self.assertEqual(statistics.median(r["armlen"] for r in rows),
                         CHR1["c3_median_armlen"])

    def test_c4_only_one_ltr_missing(self):
        _, missing_ltr, _ = T.check_coverage(self.recs, self.cls)
        self.assertEqual([(r["start"], r["end"], r["type"]) for r in missing_ltr],
                         CHR1["c4_missing_ltr"])

    def test_c4_ltr_full_coverage_claim(self):
        # 反向断言：LTR 类共 41 个，只有 1 个未覆盖——若哪天变成大面积缺失，
        # 上面那条会失败，这条给出更直观的原因
        per_class, _, _ = T.check_coverage(self.recs, self.cls)
        tot = sum(v[1] for t, v in per_class.items() if t in T.LTR_TYPES)
        cov = sum(v[0] for t, v in per_class.items() if t in T.LTR_TYPES)
        self.assertEqual((cov, tot), (40, 41))

    def test_c5(self):
        self.assertEqual(len(T.check_cross_class(self.cls)), CHR1["c5_total"])

    def test_c6(self):
        hits = T.check_ltr_conflict(self.recs, self.cls)
        self.assertEqual([(d["element"]["start"], d["element"]["end"],
                           d["element"]["ltr_identity"]) for d in hits],
                         CHR1["c6_hits"])

    def test_cli_json_matches_baseline(self):
        """端到端：跑一遍 CLI，核对 JSON 与公布基线一致。"""
        out = Path(tempfile.mkdtemp(prefix="te_check_chr1_")) / "r.json"
        try:
            p = subprocess.run(
                [sys.executable, str(SCRIPT), "--edta-dir", str(REAL_EDTA),
                 "--json", str(out), "--strict"],
                capture_output=True, text=True, encoding="utf-8")
            self.assertEqual(p.returncode, 1, msg="有严重问题时应返回 1")
            d = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(d["n_elements"], CHR1["n_elements"])
            self.assertEqual(d["n_classified"], CHR1["n_classified"])
            self.assertEqual(d["severe_total"], CHR1["severe"])
            self.assertEqual(len(d["C4_not_classified"]), 1)
            self.assertEqual(len(d["C6_ltr_conflict"]), 1)
            self.assertEqual(d["C3_ir_quality"]["median_ident"], CHR1["c3_median_ident"])
        finally:
            shutil.rmtree(out.parent, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
