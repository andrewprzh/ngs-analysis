"""Where do alt reads of FP sites come from, relative to their source transcript? Usage: screen_mechanism.py PREFIX"""
import sys
from collections import Counter
import pandas as pd
import numpy
sys.argv = sys.argv[:2]
prefix = sys.argv[1]
GTF = "/abga/work/andreyp/data/reference/mouse/gencode.vM36.basic.annotation.gtf"
tx = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    v = line.split("\t")
    if v[2] != "exon":
        continue
    tid = v[8].split('transcript_id "')[1].split('"')[0]
    tx.setdefault(tid, [v[0], v[6], []])[2].append((int(v[3]) - 1, int(v[4])))
for t in tx.values():
    t[2].sort()
sites = pd.read_csv(prefix + ".sites.tsv", sep="\t", low_memory=False).set_index("site")
reads = pd.read_csv(prefix + ".alt_reads.tsv", sep="\t", low_memory=False)
reads = reads[reads.group == "FP"]
rows = []
for r in reads.itertuples():
    s = sites.loc[r.site]
    p = s.pos - 1 if s.type == "SNV" else s.pos
    tid = r.read_id.split("_")[2]
    chrom, strand, exons = tx[tid]
    if chrom != s.chrom or p < exons[0][0] - 10000 or p >= exons[-1][1] + 10000:
        where, d = "elsewhere", -1
    elif any(a <= p < b for a, b in exons):
        where, d = "exonic", 0
    else:
        d = min(min(abs(p - a), abs(p - (b - 1))) for a, b in exons)
        if p < exons[0][0] or p >= exons[-1][1]:
            # beyond the transcript: upstream / downstream in transcript orientation
            up = (p < exons[0][0]) == (strand == "+")
            where = "before_5end" if up else "after_3end"
        else:
            where = "intron"
    rows.append((r.site, s.type, s.vaf, where, d, r.qdist_end, r.clip))
m = pd.DataFrame(rows, columns=["site", "type", "vaf", "placement", "dist_src_exon", "qdist_end", "clip"])
print("alt reads at screened FP sites, by placement relative to the source transcript")
print(pd.crosstab(m.type, m.placement, normalize="index").round(3).to_string())
print(pd.crosstab(m.type, m.placement).to_string())
x = m[m["placement"].isin(["intron", "before_5end", "after_3end"])]
print("\nnon-exonic placements: distance to the nearest source exon (bp) and to the alignment end (query bp), quantiles")
print(x.groupby(["type", "placement"])[["dist_src_exon", "qdist_end"]].quantile([0.25, 0.5, 0.75, 0.9]).unstack().round(0).to_string())
# per site: dominant placement
dom = m.groupby("site")["placement"].agg(lambda w: Counter(w).most_common(1)[0][0])
st = sites.loc[dom.index]
print("\nFP sites by dominant alt read placement")
print(pd.crosstab(st.type, dom.values).to_string())
