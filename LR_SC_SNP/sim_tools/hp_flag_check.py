"""Share of TP / FP indel calls flagged by homopolymer rules. Usage: hp_flag_check.py EVAL_SITES"""
import sys
import pandas as pd
import pysam
fa = pysam.FastaFile("/abga/work/andreyp/data/reference/mouse/GRCm39.primary_assembly.genome.fa")
d = pd.read_csv(sys.argv[1], sep="\t")
d = d[d.type.isin(["INS", "DEL"]) & d.label.isin(["TP", "FP"])].copy()


def max_run(seq):
    best, cur = 1, 1
    for a, b in zip(seq, seq[1:]):
        cur = cur + 1 if a == b and a != "N" else 1
        best = max(best, cur)
    return best


def near(row, k):
    chrom, pos, ref, alt = row["#variant_id"].split(":")
    s0 = int(pos)  # 0-based first base after the anchor
    end = s0 + max(len(ref) - 1, 0)
    lo, hi = s0 - k, end + k
    # extend to whole runs touching the window
    seq = fa.fetch(chrom, max(0, lo - 20), hi + 20).upper()
    off = lo - max(0, lo - 20)
    # runs overlapping [lo, hi)
    best = 0
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j > off and i < off + (hi - lo):
            best = max(best, j - i)
        i = j
    return best


for k in (0, 1, 3):
    d["near%d" % k] = d.apply(lambda r: near(r, k), axis=1)
rows = []
for name, m in [("hp_len (caller) >=4", d.homopolymer_len >= 4), ("hp_len >=5", d.homopolymer_len >= 5),
                ("hp_len >=6", d.homopolymer_len >= 6),
                ("run>=4 within 1bp", d.near1 >= 4), ("run>=5 within 1bp", d.near1 >= 5),
                ("run>=4 within 3bp", d.near3 >= 4), ("run>=5 within 3bp", d.near3 >= 5),
                ("run>=3 within 1bp", d.near1 >= 3)]:
    row = {"rule": name}
    for (l, t), x in d.groupby(["label", "type"]):
        row["%s_%s" % (l, t)] = round(float(m[x.index].mean()), 3)
    rows.append(row)
print(pd.DataFrame(rows).set_index("rule").to_string())
print(d.groupby(["label", "type"]).size().to_dict())
