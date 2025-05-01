import sys

count_dict = {}
for l in open(sys.argv[1]):
    if l.startswith("#") or l.startswith("transcript_id"):
        continue
    v = l.split("\t")
    count_dict[v[0]] = float(v[1])

scale_factor = sum(count_dict.values()) / 1000000.0
for tid in sorted(count_dict.keys()):
    print("%s\t%.2f\t%.6f" % (tid, count_dict[tid], count_dict[tid] / scale_factor))
