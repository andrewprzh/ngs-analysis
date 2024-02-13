import sys
import pysam
import numpy
from collections import defaultdict


mapq_to_de_dict = defaultdict(list)
for a in pysam.AlignmentFile(sys.argv[1], "rb"):
    if a.is_secondary: continue
    try:    
        de = a.get_tag("de")
    except KeyError:
        continue
    mapq_to_de_dict[a.mapping_quality].append(de)

for k in sorted(mapq_to_de_dict.keys()):
    de_vals = mapq_to_de_dict[k]
    print("%d\t%.4f\t%.4f\t%.4f" % (k, numpy.quantile(de_vals, 0.1),
                                    numpy.mean(de_vals), numpy.quantile(de_vals, 0.9)))