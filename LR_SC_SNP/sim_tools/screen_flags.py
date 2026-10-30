"""Candidate 'unreliable' flags: share of screened FP and TP calls each flag marks. Usage: screen_flags.py PREFIX"""
import sys
import numpy
import pandas as pd
pd.set_option("display.width", 250)
d = pd.read_csv(sys.argv[1] + ".sites.tsv", sep="\t", low_memory=False)
c = d[d.group.isin(["TP", "FP"])].copy()
# other calls within 50 bp (any type)
c["n_near50"] = 0
for chrom, x in c.groupby("chrom"):
    pos = numpy.sort(x.pos.values)
    lo = numpy.searchsorted(pos, x.pos.values - 50, "left")
    hi = numpy.searchsorted(pos, x.pos.values + 50, "right")
    c.loc[x.index, "n_near50"] = hi - lo - 1
s = c[c.screened == 1].copy()
flags = {
    "alt within 20bp of read end (>=50% alt reads)": s.alt_near_end20 >= 0.5,
    "alt within 20bp of read end (>=20% alt reads)": s.alt_near_end20 >= 0.2,
    "alt reads noisier (err >= ref err + 0.02)": s.alt_read_err >= s.ref_read_err + 0.02,
    ">=2 other calls within 50bp": s.n_near50 >= 2,
    "MAPQ<60 in >=20% alt reads": s.alt_mapq_lt60 >= 0.2,
    "several genes span (expressed)": s.genes_spanning > 1,
    "not exonic (expressed)": s.expr_tx_exonic == 0,
}
flags["any of first 4 (excl. 20% end)"] = flags["alt within 20bp of read end (>=50% alt reads)"] | \
    flags["alt reads noisier (err >= ref err + 0.02)"] | flags[">=2 other calls within 50bp"] | flags["MAPQ<60 in >=20% alt reads"]
rows = []
for name, m in flags.items():
    row = {"flag": name}
    for (g, t), x in s.groupby(["group", "type"]):
        row["%s_%s" % (g, t)] = round(float(m[x.index].mean()), 3)
    rows.append(row)
print("share of screened calls flagged")
print(pd.DataFrame(rows).set_index("flag").to_string())
print("screened:", s.groupby(["group", "type"]).size().to_dict())
anyf = flags["any of first 4 (excl. 20% end)"]
f = s[s.group == "FP"].copy()
f["vbin"] = pd.cut(f.vaf, [0, 0.02, 0.1, 0.5, 0.9, 1.01])
print("\nFP flagged (any of first 4) by VAF")
print(f.assign(flag=anyf[f.index]).groupby(["type", "vbin"], observed=True).flag.agg(["mean", "size"]).round(3).to_string())
print("\nunflagged FP: VAF, depth, homopolymer (medians) and share with vaf<0.02")
u = f[~anyf[f.index]]
print(u.groupby("type").agg(n=("site", "size"), vaf=("vaf", "median"), depth=("depth_umis", "median"), hp=("hp", "median"),
                           low_vaf=("vaf", lambda v: (v < 0.02).mean())).round(3).to_string())
