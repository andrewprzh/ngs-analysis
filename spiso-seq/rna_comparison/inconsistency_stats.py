import numpy
import sys
from collections import defaultdict

match_types = sorted(["extra_intron_novel", "fake_terminal_exon_3", "fake_terminal_exon_5", "extra_intron_5", "extra_intron_3"])

read_info_map = {}
for l in open(sys.argv[1]):
    t = l.strip().split()
    if t[2] != "inconsistent":
        continue
    read_id = t[0]
    intron_count = t[4].count(",")
    if intron_count not in read_info_map:
        read_info_map[intron_count] = defaultdict(int)

    for mt in match_types:
        if t[3].find(mt) != -1:
            read_info_map[intron_count][mt] += 1
    read_info_map[intron_count]["total_reads"] += 1

print("introns\t" + "\t".join(match_types))
for intron_count in sorted(read_info_map.keys()):
    print("%d\t%s" % (intron_count, "\t".join([str(read_info_map[intron_count][x]) for x in ["total_reads"] + match_types])))

print("introns\t" + "\t".join(match_types))
for intron_count in sorted(read_info_map.keys()):
    total_reads = read_info_map[intron_count]["total_reads"]
    fractions = [(read_info_map[intron_count][x] / total_reads) for x in match_types]
    print("%d\t%s" % (intron_count, "\t".join(["%.4f" % x for x in fractions])))



