#!/usr/bin/env python3
"""
通用元件判读：给定 chr1 坐标，六框翻译 + REXdb v4.0 hmmsearch，输出域组成与判定。

用法：
    python classify_element.py 17528667 17536138            # 目标元件
    python classify_element.py 11148091 11152835 --label 对照  # 阳性对照

判据参考 TEsorter：命中模型名中的层级即为归属
    Class_I/LTR/Ty1_copia/<clade>:<domain>
    Class_I/LTR/Ty3_gypsy/<clade>:<domain>
"""
import sys
import re
import argparse
from collections import defaultdict

import pyhmmer

FA = "//wsl.localhost/Ubuntu-26.04/home/sean/edta_run/chr1.fa.mod"
HMM = ("//wsl.localhost/Ubuntu-26.04/home/sean/miniconda3/envs/edta/"
       "lib/python3.12/site-packages/TEsorter/database/"
       "REXdb_protein_database_viridiplantae_v4.0.hmm")
MIN_ORF_AA = 30

_B = "TCAG"
_A = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {}
_i = 0
for _x in _B:
    for _y in _B:
        for _z in _B:
            CODON[_x + _y + _z] = _A[_i]
            _i += 1
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CACHE = {}


def _load(seqid="1"):
    if seqid not in _CACHE:
        buf, want = [], False
        with open(FA, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    if want:
                        break
                    want = line[1:].split()[0] == seqid
                    continue
                if want:
                    buf.append(line.strip())
        _CACHE[seqid] = "".join(buf)
    return _CACHE[seqid]


def extract(start, end, seqid="1"):
    return _load(seqid)[start - 1:end]


def translate(seq, frame):
    return "".join(CODON.get(seq[i:i + 3].upper(), "X")
                   for i in range(frame, len(seq) - 2, 3))


def make_peptides(dna):
    """双链六框翻译，按终止切分，保留 >= MIN_ORF_AA 的片段。"""
    rc = dna.translate(COMP)[::-1]
    out, stats = [], {}
    for label, seq in (("plus", dna), ("minus", rc)):
        for f in range(3):
            prot = translate(seq, f)
            n_stop = prot.count("*")
            longest = max((len(x) for x in prot.split("*")), default=0)
            stats[f"{label}_f{f}"] = (len(prot), n_stop, longest)
            for k, chunk in enumerate(prot.split("*")):
                chunk = chunk.strip("X")
                if len(chunk) >= MIN_ORF_AA:
                    out.append((f"{label}_f{f}_s{k}", chunk))
    return out, stats


def search(peptides):
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=n.encode(), sequence=s).digitize(alphabet)
            for n, s in peptides]
    with pyhmmer.plan7.HMMFile(HMM) as hf:
        hmms = list(hf)
    rows = []
    for top in pyhmmer.hmmsearch(hmms, seqs):
        q = top.query
        qn = q.name.decode() if isinstance(q.name, bytes) else str(q.name)
        for hit in top:
            d = hit.best_domain
            ev = d.i_evalue if d is not None else 1.0
            if ev > 1e-3:
                continue
            rows.append((qn, hit.score, ev,
                         hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)))
    rows.sort(key=lambda x: -x[1])
    return rows


def report(label, dna, rows, stats):
    print(f"\n{'=' * 78}\n{label}  长度 {len(dna)} bp\n{'=' * 78}")
    print("终止密码子密度（判断是否还有编码能力）:")
    for k, (L, ns, lg) in stats.items():
        print(f"  {k}: {ns:3d}/{L:4d} = {ns / L * 100:4.1f}%  最长无终止段 {lg:3d} aa")
    if not rows:
        print("\n域比对结果: 无命中（E < 1e-3）")
        return
    print(f"\n域比对结果: {len(rows)} 条命中")
    dom = defaultdict(list)
    for name, score, ev, seq in rows:
        parts = name.split("/")
        branch = parts[2] if len(parts) > 2 else "?"
        clade = parts[3].split(":")[0] if len(parts) > 3 else "?"
        d = name.split(":")[-1]
        dom[branch].append((d, score, clade))
    for branch in sorted(dom, key=lambda b: -max(s for _, s, _ in dom[b])):
        items = sorted(dom[branch], key=lambda x: -x[1])
        print(f"\n  [{branch}] {len(items)} 个域")
        seen = set()
        for d, s, clade in items:
            k = (d, clade)
            if k in seen:
                continue
            seen.add(k)
            print(f"    {s:7.1f}  {d:<6} clade={clade}")
    # 判定
    agg = defaultdict(float)
    for b, items in dom.items():
        agg[b] = sum(s for _, s, _ in items)
    print("\n  分支累计得分:")
    for b, s in sorted(agg.items(), key=lambda x: -x[1]):
        print(f"    {b:<12} {s:8.1f}")
    if agg:
        win = max(agg, key=agg.get)
        print(f"\n  判定: {win}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("start", type=int)
    ap.add_argument("end", type=int)
    ap.add_argument("--label", default="目标元件")
    a = ap.parse_args()
    dna = extract(a.start, a.end)
    pep, stats = make_peptides(dna)
    print(f"六框翻译片段数: {len(pep)}（>= {MIN_ORF_AA} aa）")
    rows = search(pep)
    report(a.label, dna, rows, stats)


if __name__ == "__main__":
    main()
