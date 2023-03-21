#chr1_853391_854385_ENSG00000228794.10_+ 0       1       1       CROPseq Splicing_factor SRSF5   CROPseq_SF_sg203        AAAGGATAGGTGCGAT

import sys
from collections import defaultdict


inf = sys.argv[1]
tf = None
if len(sys.argv) > 2:
    tf = sys.argv[2]
tf_col = 6
REST = "other"

exon_dict = defaultdict(lambda: defaultdict(int))
for l in open(inf):
    v = l.strip().split()
    if tf is None:
        tf = v[tf_col]
    cur_tf = v[tf_col]
    exon_id = v[0]
    count = int(v[3])
    if cur_tf == tf:
        exon_dict[exon_id][cur_tf] += count
    else:
        exon_dict[exon_id][REST] += count


for exon_id in exon_dict:
    print("%d\t%d" % (exon_dict[exon_id][REST], exon_dict[exon_id][tf]))