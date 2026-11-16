"""
tests.py —— normalize.py 的断言测试集（标准库 unittest，零依赖）

为什么要写它：
    在此之前，改一次代码只能靠肉眼看 7 个文件的输出摘要。今天写第一版断言时，
    立刻抓到一个真 bug：describe_so("non_LTR_retrotransposon") 的 superfamily
    仍然是 "non"——因为它以 _LTR_retrotransposon 结尾，先被 LTR 分支截走了，
    之前放在 LINE 分支里的"修复"是死代码。肉眼看不出来，断言一行就露。

测试分五层：
    1. 单元：classify_format / lookup / describe_so / relate / read_names
    2. 不变量：所有非未命中结果，命中键必须真能在别名表里查到（防"自说自话"）
    3. 数据：105 个 SO 名必须全部归入某个类别，不能有 unknown
    4. 回归：7 个真实测试文件的首行摘要必须与记录一致（改动前后对比）
    5. CLI：--verbose 的输出格式、CSV 列数、默认输出不变

怎么跑：
    python tests.py                 # 简洁模式
    python tests.py -v              # 每条用例名
    python -m unittest tests -v     # 同上
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))

import normalize as N                                     # noqa: E402

ALIAS, SO_NAMES, CONFLICTS = N.load_te_so(str(BASE / "data" / "TE_Sequence_Ontology.txt"))
EX = BASE / "examples"

# 角色词表：describe_so 只允许返回这些值
VALID_ROLES = {"element", "fragment", "sub_feature", "region",
               "simple_repeat", "organelle", "RNA_gene", "unknown"}


def names_of(fname):
    """读 examples 下的名字清单（默认模式）"""
    return N.read_names(str(EX / fname))


def counts_of(fname):
    """按 lookup 状态统计一个文件，返回 {状态: 数量}（只保留非零项）"""
    tally = {}
    for nm in names_of(fname):
        st = N.lookup(nm, ALIAS, SO_NAMES)["status"]
        tally[st] = tally.get(st, 0) + 1
    return {k: v for k, v in tally.items() if v}


# ---------------------------------------------------------------------------
# 第 1 层：单元测试
# ---------------------------------------------------------------------------
class TestLoadTeSo(unittest.TestCase):
    """加载器的边界：别名拆分、空格分行兜底、别名冲突记账。"""

    def _load(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False,
                                         encoding="utf-8", newline="\n") as f:
            f.write(text)
            path = f.name
        try:
            return N.load_te_so(path)
        finally:
            Path(path).unlink()

    def test_split_by_comma_and_space(self):
        alias, names, _ = self._load("SO_A\t0001\ta1, a2 a3\n")   # 逗号与空格混用
        self.assertEqual(set(alias), {"a1", "a2", "a3"})
        self.assertEqual(alias["a2"], ("SO_A", "0001"))
        self.assertIn("SO_A", names)

    def test_space_separated_line_fallback(self):
        # 整行用空格而非制表符：按 \t 只能切出 1 列，必须靠 split() 兜底救回
        alias, _, _ = self._load("SO_B 0002 b1  b2\n")
        self.assertEqual(alias["b1"], ("SO_B", "0002"))
        self.assertEqual(alias["b2"], ("SO_B", "0002"))

    def test_empty_alias_dropped(self):
        # 连续逗号会产生空串，若不滤掉，"" 就成了一个合法键
        alias, _, _ = self._load("SO_A\t0001\ta1,,a2\n")
        self.assertNotIn("", alias)
        self.assertEqual(set(alias), {"a1", "a2"})

    def test_bad_line_skipped(self):
        alias, _, _ = self._load("SO_A\nSO_B\t0002\tb\n")      # 第一行不足 3 列
        self.assertNotIn("SO_A", alias)
        self.assertEqual(alias["b"], ("SO_B", "0002"))

    def test_alias_conflict_is_reported(self):
        # 同一别名指向不同 SO 名：后者覆盖前者（保留原策略），但必须记账
        alias, _, conflicts = self._load("SO_A\t0001\tshared\nSO_B\t0002\tshared\n")
        self.assertEqual(alias["shared"], ("SO_B", "0002"))
        self.assertIn("shared", conflicts)
        self.assertEqual([n for n, _ in conflicts["shared"]], ["SO_A", "SO_B"])

    def test_same_target_repeated_is_not_conflict(self):
        _, _, conflicts = self._load("SO_A\t0001\tx\nSO_A\t0001\tx\n")
        self.assertEqual(conflicts, {})

    def test_three_way_conflict_lists_all(self):
        _, _, conflicts = self._load(
            "SO_A\t0001\ts\nSO_B\t0002\ts\nSO_C\t0003\ts\n")
        self.assertEqual([n for n, _ in conflicts["s"]], ["SO_A", "SO_B", "SO_C"])

    def test_alias_equal_to_other_so_name_not_reported(self):
        # 别名与另一行的标准名同名、但只作为别名出现一次：不是重复定义。
        # 它确实会影响查表（lookup 路径1 会命中别名），但那是 lookup 的语义问题。
        _, _, conflicts = self._load("SO_A\t0001\tSO_B\nSO_B\t0002\tz\n")
        self.assertEqual(conflicts, {})

    def test_real_table_has_no_conflict(self):
        # 回归锁：随仓库分发的 TE_SO 表当前零冲突。若哪天变成非空，
        # 要么表被改出了重复定义，要么该去修它——都不该静默通过。
        self.assertEqual(CONFLICTS, {})
        self.assertGreater(len(ALIAS), 700)


class TestClassifyFormat(unittest.TestCase):
    """命名格式判定：standard / clade_style / malformed"""

    def test_standard_two_levels(self):
        self.assertEqual(N.classify_format("DNA/hAT"), "standard")
        self.assertEqual(N.classify_format("LTR/Copia"), "standard")

    def test_clade_three_levels(self):
        # 两种层级顺序都算 clade_style：RepeatMasker 风格与 TEsorter 风格
        self.assertEqual(N.classify_format("LTR/Copia/SIRE"), "clade_style")
        self.assertEqual(N.classify_format("ClassI/LTR/Gypsy"), "clade_style")

    def test_malformed(self):
        self.assertEqual(N.classify_format("Unknown//ClassI"), "malformed")  # 空级
        self.assertEqual(N.classify_format("A/B/C/D"), "malformed")          # 超过三级


class TestLookup(unittest.TestCase):
    """五条查找路径，重点验证 hit（表里实际命中的键）"""

    def test_path1_exact(self):
        r = N.lookup("DNA/hAT", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched")
        self.assertEqual(r["so_id"], "SO:0002279")
        self.assertEqual(r["hit"], "DNA/hAT")          # 精确命中时 hit 就是输入

    def test_path2_case_insensitive(self):
        r = N.lookup("dna/hat", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched_ci")
        self.assertEqual(r["hit"], "DNA/hAT")           # 表里写法与输入不同
        self.assertNotEqual(r["hit"], "dna/hat")

    def test_path2_real_case(self):
        # 真实案例：用户写了 sine/alu，表里是 SINE/Alu
        r = N.lookup("sine/alu", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched_ci")
        self.assertEqual(r["hit"], "SINE/Alu")

    def test_path3_strip_repeatmasker_style(self):
        # 超家族在前：砍掉最后一级
        r = N.lookup("LTR/Copia/SIRE", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched_after_strip")
        self.assertEqual(r["hit"], "LTR/Copia")

    def test_path3_strip_tesorter_style(self):
        # 超家族在最后：取最后两级。只砍最后一级会得到 ClassI/LTR，永远查不到。
        r = N.lookup("ClassI/LTR/Gypsy", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched_after_strip")
        self.assertEqual(r["hit"], "LTR/Gypsy")

    def test_path4_separator(self):
        r = N.lookup("DNA/CMC_Transib", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "matched_norm")
        self.assertEqual(r["hit"], "DNA/CMC-Transib")

    def test_path5_unknown_never_discarded(self):
        r = N.lookup("ZZZ_not_a_te_name", ALIAS, SO_NAMES)
        self.assertEqual(r["status"], "unknown")
        self.assertEqual(r["hit"], "")                  # 未命中没有命中键
        self.assertEqual(r["so_name"], "")
        # 关键：必须给建议，不能静默丢弃
        self.assertTrue(r["suggestion"], "未命中时必须给出相似建议")


class TestDescribeSo(unittest.TestCase):
    """层级归类：类别(order) / 角色(role) / 超家族(superfamily)"""

    CASES = [
        # SO 名                        order       role          superfamily
        ("helitron",                   "Class II", "element",    "helitron"),
        ("hAT_TIR_transposon",         "Class II", "element",    "hAT"),
        ("Mutator_TIR_transposon",     "Class II", "element",    "Mutator"),
        ("Copia_LTR_retrotransposon",  "Class I",  "element",    "Copia"),
        ("Gypsy_LTR_retrotransposon",  "Class I",  "element",    "Gypsy"),
        ("L1_LINE_retrotransposon",    "Class I",  "element",    "L1"),
        ("long_terminal_repeat",       "-",        "sub_feature", "long_terminal_repeat"),
        ("repeat_region",              "-",        "region",     "repeat_region"),
        ("repeat_fragment",            "-",        "fragment",   "repeat_fragment"),
        ("cytosolic_18S_rRNA",         "-",        "RNA_gene",   "cytosolic_18S_rRNA"),
    ]

    def test_layer_cases(self):
        for so, order, role, sf in self.CASES:
            with self.subTest(so=so):
                d = N.describe_so(so)
                self.assertEqual(d["order"], order)
                self.assertEqual(d["role"], role)
                self.assertEqual(d["superfamily"], sf)

    def test_non_ltr_superfamily(self):
        """回归用例：non_LTR_retrotransposon 曾被截成假超家族 'non'。

        它同时以 _LTR_retrotransposon 结尾，若 LTR 分支先执行就会被截走，
        所以这个断言守住"分支顺序"这个隐性约定。
        """
        d = N.describe_so("non_LTR_retrotransposon")
        self.assertEqual(d["order"], "Class I")
        self.assertEqual(d["superfamily"], "non_LTR")

    def test_all_so_names_classified(self):
        """105 个 SO 标准名必须全部归入某个类别，0 个 unknown"""
        bad = [s for s in SO_NAMES if N.describe_so(s)["order"] == "unknown"]
        self.assertEqual(bad, [], f"有 {len(bad)} 个 SO 名未归类：{bad[:5]}")

    def test_role_vocabulary(self):
        """role 只能取词表内的值，防止以后随手加出拼写不一致的新词"""
        got = {N.describe_so(s)["role"] for s in SO_NAMES}
        self.assertTrue(got <= VALID_ROLES, f"出现了词表外的 role：{got - VALID_ROLES}")


class TestRelate(unittest.TestCase):
    """两个名字之间的关系（--pairs 模式的核心）"""

    def judge(self, a, b):
        return N.relate(N.lookup(a, ALIAS, SO_NAMES), N.lookup(b, ALIAS, SO_NAMES))[0]

    def test_equivalent(self):
        # DNA/DTC 与 DNA/CACTA 是同一概念两种写法
        self.assertEqual(self.judge("DNA/DTC", "DNA/CACTA"), "等价")

    def test_same_order_diff_family(self):
        # 都是 Class II，但 hAT 与 Mutator 不同科 —— EDTA 只记录不处置的那一类
        self.assertEqual(self.judge("DNA/hAT", "DNA/Mutator"), "同目不同科")

    def test_cross_class(self):
        self.assertEqual(self.judge("LTR/Gypsy", "DNA/Tc1_Mariner"), "跨类冲突")

    def test_role_mismatch(self):
        # 元件 vs 零件，不可直接比较
        self.assertEqual(self.judge("LTR/Copia", "long_terminal_repeat"), "角色不同")

    def test_undecidable(self):
        self.assertEqual(self.judge("ZZZ_not_a_te_name", "DNA/hAT"), "无法判定")


class TestRankPairs(unittest.TestCase):
    """加权排序（--pairs --weight）：镜像对必须合并，方向必须保留"""

    def rec(self, a, b, w):
        ra = N.lookup(a, ALIAS, SO_NAMES)
        rb = N.lookup(b, ALIAS, SO_NAMES)
        return {"a": a, "a_so": ra["so_name"], "a_id": ra["so_id"],
                "a_order": N.describe_so(ra["so_name"])["order"],
                "b": b, "b_so": rb["so_name"], "b_id": rb["so_id"],
                "b_order": N.describe_so(rb["so_name"])["order"],
                "verdict": N.relate(ra, rb)[0], "why": "", "w": w}

    def test_mirror_pairs_merged(self):
        # Tc1→PIF 7 与 PIF→Tc1 1 是同一对分歧的两个方向，合并后应为 8
        ranked = N.rank_pairs([self.rec("DNA/Tc1_Mariner", "DNA/PIF_Harbinger", 7),
                               self.rec("DNA/PIF_Harbinger", "DNA/Tc1_Mariner", 1)])
        self.assertEqual(len(ranked), 1)
        self.assertEqual(ranked[0]["weight"], 8)

    def test_direction_kept(self):
        # 合并后仍要能看出谁把谁吸走：7 在前、1 在后
        ranked = N.rank_pairs([self.rec("DNA/Tc1_Mariner", "DNA/PIF_Harbinger", 7),
                               self.rec("DNA/PIF_Harbinger", "DNA/Tc1_Mariner", 1)])
        self.assertTrue(ranked[0]["dirs"].startswith("Tc1_Mariner→PIF_Harbinger 7"))

    def test_sorted_by_weight(self):
        recs = [self.rec("DNA/CACTA", "DNA/Mutator", 1),
                self.rec("LTR/Copia", "DNA/Mutator", 5),
                self.rec("DNA/hAT", "DNA/Mutator", 3)]
        weights = [r["weight"] for r in N.rank_pairs(recs)]
        self.assertEqual(weights, sorted(weights, reverse=True))

    def test_equivalent_pair_single_group(self):
        # 等价对也必须合成一组，否则同一件事会被数两次
        ranked = N.rank_pairs([self.rec("DNA/DTC", "DNA/CACTA", 2),
                               self.rec("DNA/CACTA", "DNA/DTC", 2)])
        self.assertEqual(len(ranked), 1)
        self.assertEqual(ranked[0]["weight"], 4)


class TestReadNames(unittest.TestCase):
    """read_names 的三种真实行格式"""

    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix=".txt", text=True)
        os.close(fd)
        Path(path).write_text(text, encoding="utf-8")
        self.addCleanup(os.remove, path)
        return path

    def test_name_then_count(self):
        p = self._write("DNA/hAT 104231\n")
        self.assertEqual(N.read_names(p), ["DNA/hAT"])

    def test_count_then_name(self):
        # uniq -c 的输出：数量在前，不能把数字当名字
        p = self._write("1461 DNA/hAT-hAT19\n")
        self.assertEqual(N.read_names(p), ["DNA/hAT-hAT19"])

    def test_pairs_with_count(self):
        # EDTA 混淆矩阵：两个名字 + 次数。第三个数是权重，不是名字，必须原样带出来
        p = self._write("DNA/Tc1_Mariner  DNA/PIF_Harbinger  7\n")
        self.assertEqual(N.read_names(p, pairs=True),
                         [("DNA/Tc1_Mariner", "DNA/PIF_Harbinger", 7)])

    def test_pairs_without_count(self):
        # 没有次数列的配对文件：权重记 1，不能记 0（0 会让统计凭空消失）
        p = self._write("DNA/hAT\tDNA/Mutator\n")
        self.assertEqual(N.read_names(p, pairs=True), [("DNA/hAT", "DNA/Mutator", 1)])

    def test_pairs_count_in_front(self):
        # 次数写在行首也要认（uniq -c 风格）
        p = self._write("3  LTR/Gypsy  DNA/Mutator\n")
        self.assertEqual(N.read_names(p, pairs=True), [("LTR/Gypsy", "DNA/Mutator", 3)])

    def test_skip_comment_and_blank(self):
        p = self._write("# 注释\n\nDNA/hAT\n")
        self.assertEqual(N.read_names(p), ["DNA/hAT"])

    # --- 解析失败必须可见：静默丢行比读错更难查
    def test_all_digits_line_skipped(self):
        # 整行都是数字＝没有名字。旧实现用 tokens[-1:] 兜底，会把数字当名字造出来
        p = self._write("12345\n")
        issues = []
        self.assertEqual(N.read_names(p, issues=issues), [])
        self.assertEqual(len(issues), 1)
        self.assertIn("没有名字", issues[0][2])

    def test_pairs_single_name_reported(self):
        # 配对模式遇到单名字行：整行丢弃，但必须报出来（旧实现让它静默消失）
        p = self._write("DNA/A\n")
        issues = []
        self.assertEqual(N.read_names(p, pairs=True, issues=issues), [])
        self.assertEqual(len(issues), 1)
        self.assertIn("只有 1 个", issues[0][2])

    def test_single_mode_extra_name_reported(self):
        # 单列模式遇到两个名字：只取第一个，第二个的去向必须可见
        p = self._write("DNA/A  DNA/B\n")
        issues = []
        self.assertEqual(N.read_names(p, issues=issues), ["DNA/A"])
        self.assertEqual(len(issues), 1)
        self.assertIn("DNA/B", issues[0][2])

    def test_pairs_extra_name_reported(self):
        p = self._write("DNA/A  DNA/B  DNA/C\n")
        issues = []
        self.assertEqual(N.read_names(p, pairs=True, issues=issues),
                         [("DNA/A", "DNA/B", 1)])
        self.assertEqual(len(issues), 1)

    def test_normal_formats_produce_no_issue(self):
        # 三种正常格式都不能产生 issue，否则警告变噪音，真问题反而被忽略
        p = self._write("DNA/hAT 104231\n104231 DNA/hAT\nDNA/hAT\n")
        issues = []
        self.assertEqual(len(N.read_names(p, issues=issues)), 3)
        self.assertEqual(issues, [])

    def test_issue_carries_lineno_and_text(self):
        p = self._write("# 注释\nDNA/A\n\n12345\n")
        issues = []
        N.read_names(p, issues=issues)
        self.assertEqual(issues[0][0], 4)              # 1-based，注释与空行不计数
        self.assertEqual(issues[0][1], "12345")


# ---------------------------------------------------------------------------
# 第 2 层：不变量（比"输出对不对"更重要的性质）
# ---------------------------------------------------------------------------
class TestInvariants(unittest.TestCase):
    """对所有真实数据都成立的硬性要求"""

    FIXTURES = ["test_names.txt", "test_names_counts.txt", "edta_chr1_TElib_names.txt",
                "matrix_labels.txt", "gff3_col3.txt", "gff3_class.txt", "tesorter_style.txt"]

    def test_hit_is_a_real_alias(self):
        """凡是非未命中的结果，hit 必须真能在别名表里查到，且指向同一个 SO。

        这条防的是"程序自己编一个命中名"——如果 hit 只是原样回显输入，
        这一条就会立刻失败（因为输入本身并不在表里）。
        """
        checked = 0
        for f in self.FIXTURES:
            for nm in names_of(f):
                r = N.lookup(nm, ALIAS, SO_NAMES)
                if r["status"] == "unknown":
                    continue
                self.assertIn(r["hit"], ALIAS, f"{nm} 的命中键 {r['hit']} 不在别名表里")
                self.assertEqual(ALIAS[r["hit"]], (r["so_name"], r["so_id"]))
                if r["status"] == "matched":
                    self.assertEqual(r["hit"], nm)
                else:
                    self.assertNotEqual(r["hit"], nm,
                                        f"{nm} 走了 {r['status']}，命中键却与输入相同")
                checked += 1
        self.assertGreater(checked, 50, "测试数据太少，不变量没被真正验证")

    def test_reflexivity(self):
        """输出必须能被自己再吃回去：产出的 SO 名要能作为输入再次命中"""
        for f in self.FIXTURES:
            for nm in names_of(f):
                so = N.lookup(nm, ALIAS, SO_NAMES)["so_name"]
                if not so:
                    continue
                back = N.lookup(so, ALIAS, SO_NAMES)
                self.assertEqual(back["status"], "matched",
                                 f"输出的 {so} 无法被自己再查回来")


# ---------------------------------------------------------------------------
# 第 3 层：回归（数字是记录下来的，不是算出来的）
# ---------------------------------------------------------------------------
class TestFixtureRegression(unittest.TestCase):
    """7 个真实文件的首行摘要。改动前后对比，任何行为漂移都会在这里暴露。"""

    EXPECTED = {
        # 文件名                         精确  大小写  剥级  连字符  未知
        "test_names.txt":               (4,   2,     1,    1,     2),
        "test_names_counts.txt":        (4,   0,     1,    0,     0),
        "edta_chr1_TElib_names.txt":    (14,  0,     0,    0,     0),
        "matrix_labels.txt":            (12,  0,     0,    0,     0),   # EDTA 混淆矩阵标签
        "gff3_col3.txt":                (15,  0,     0,    0,     0),   # gff3 第 3 列 SO 名
        "gff3_class.txt":               (14,  0,     0,    0,     0),   # gff3 classification
        "tesorter_style.txt":           (0,   0,     5,    0,     0),   # TEsorter 层级命名
    }

    def test_counts(self):
        for fname, exp in self.EXPECTED.items():
            with self.subTest(file=fname):
                c = counts_of(fname)
                got = (c.get("matched", 0), c.get("matched_ci", 0),
                       c.get("matched_after_strip", 0), c.get("matched_norm", 0),
                       c.get("unknown", 0))
                self.assertEqual(got, exp)


# ---------------------------------------------------------------------------
# 第 4 层：CLI（--verbose 与 CSV 输出）
# ---------------------------------------------------------------------------
class TestCli(unittest.TestCase):

    def run_cli(self, *args):
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        return subprocess.run([sys.executable, str(BASE / "normalize.py"), *args],
                              cwd=str(BASE), capture_output=True, text=True,
                              encoding="utf-8", env=env)

    def test_default_output_unchanged(self):
        r = self.run_cli(str(EX / "test_names.txt"))
        self.assertEqual(r.returncode, 0)
        first = r.stdout.splitlines()[0]
        self.assertEqual(first, "共 10 个名称：精确命中 4 ｜ 大小写命中 2 ｜ "
                                "剥级后命中 1 ｜ 连字符归一 1 ｜ 未知 2")

    def test_verbose_shows_hit_key(self):
        r = self.run_cli(str(EX / "test_names.txt"), "--verbose")
        self.assertEqual(r.returncode, 0)
        out = r.stdout
        self.assertIn("命中键", out)                 # 表头必须有这一列
        self.assertIn("路径2 大小写", out)
        self.assertIn("路径3 剥级", out)
        self.assertIn("路径4 连字符", out)
        # sine/alu 这一行必须同时出现原始名与真实命中键
        line = [l for l in out.splitlines() if "sine/alu" in l]
        self.assertTrue(line, "verbose 输出里找不到 sine/alu 的明细行")
        self.assertIn("SINE/Alu", line[0])

    def test_verbose_all_exact_no_empty_table(self):
        """全部精确命中时不该画出空表格"""
        r = self.run_cli(str(EX / "edta_chr1_TElib_names.txt"), "--verbose")
        self.assertEqual(r.returncode, 0)
        self.assertIn("全部走路径1", r.stdout)
        self.assertNotIn("SO 标准名 (SO ID)", r.stdout)

    def test_csv_columns(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        self.addCleanup(os.remove, path)
        r = self.run_cli(str(EX / "test_names.txt"), "--out", path)
        self.assertEqual(r.returncode, 0)
        header = Path(path).read_text(encoding="utf-8-sig").splitlines()[0]
        cols = header.split(",")
        self.assertEqual(len(cols), 10, f"CSV 应为 10 列，实际 {len(cols)}：{cols}")
        self.assertEqual(cols[3], "命中键")

    def test_weight_output(self):
        """--pairs --weight：必须报出合并后的对数，并给出加权占比"""
        r = self.run_cli(str(EX / "pairs_matrix.txt"), "--pairs", "--weight")
        self.assertEqual(r.returncode, 0)
        out = r.stdout
        # 8 行 → 6 对（Tc1/PIF 与 Gypsy/Mutator 各合并一次）→ 17 条
        self.assertIn("共 8 行记录，合并镜像后 6 对", out)
        self.assertIn("涉及序列 17 条", out)
        self.assertIn("64.7%", out)          # 同目不同科加权后占比
        self.assertIn("35.3%", out)          # 跨类冲突加权后占比
        # 方向明细不能因为合并而丢失
        self.assertIn("Tc1_Mariner→PIF_Harbinger 7", out)

    def test_weight_requires_pairs(self):
        """--weight 单独用是误用，应当直接报错而不是静默忽略"""
        r = self.run_cli(str(EX / "test_names.txt"), "--weight")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--pairs", r.stderr)

    def test_pairs_csv_has_count_column(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        self.addCleanup(os.remove, path)
        r = self.run_cli(str(EX / "pairs_matrix.txt"), "--pairs", "--weight", "--out", path)
        self.assertEqual(r.returncode, 0)
        lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
        self.assertEqual(len(lines[0].split(",")), 11, "配对 CSV 应为 11 列（含次数）")
        self.assertTrue(lines[0].endswith("次数"))
        # 原始行不合并，审计用；权重必须落到列里
        self.assertEqual(len(lines) - 1, 8)
        self.assertTrue(any(line.endswith(",7") for line in lines[1:]))


# ---------------------------------------------------------------------------
# 第 5 层：显示工具（中文对齐 / 截断）
# ---------------------------------------------------------------------------
class TestDisplayHelpers(unittest.TestCase):

    def test_disp_width(self):
        self.assertEqual(N.disp_width("abc"), 3)
        self.assertEqual(N.disp_width("中文"), 4)      # 全角算 2 列
        self.assertEqual(N.disp_width("路径2 大小写"), 12)

    def test_pad(self):
        self.assertEqual(N.disp_width(N.pad("中文", 10)), 10)
        self.assertEqual(N.disp_width(N.pad("abc", 5)), 5)

    def test_clip(self):
        long = "DNA/ThisIsAVeryLongUnknownSuperfamilyNameHere"
        out = N.clip(long, 25)
        self.assertTrue(out.endswith("…"))
        self.assertLessEqual(N.disp_width(out), 25)
        # 短名字不动
        self.assertEqual(N.clip("DNA/hAT", 25), "DNA/hAT")


if __name__ == "__main__":
    unittest.main(verbosity=2)
