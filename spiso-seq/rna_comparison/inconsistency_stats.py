import numpy
import sys
from collections import defaultdict


def print_stats(read_intron_dict, max_len = 5):
    # sum up everythin above max_len
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count <= max_len:
            continue
        for k in read_intron_dict[intron_count].keys():
            read_intron_dict[max_len][k] += read_intron_dict[intron_count][k]
            read_intron_dict[intron_count][k] = 0

    print("introns\ttotal\t" + "\t".join(match_types))
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count > max_len:
            continue
        print("%d\t%s" % (
        intron_count, "\t".join([str(read_intron_dict[intron_count][x]) for x in ["total_reads"] + match_types])))

    print("introns\t" + "\t".join(match_types))
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count > max_len:
            continue
        total_reads = read_intron_dict[intron_count]["total_reads"]
        fractions = [(100 * read_intron_dict[intron_count][x] / total_reads) for x in match_types]
        print("%d\t%s" % (intron_count, "\t".join(["%.2f" % x for x in fractions])))


match_types = (["extra_intron_novel", "fake_terminal_exon_3", "fake_terminal_exon_5", "extra_intron_5", "extra_intron_3", "alternative_structure_novel"])


supported_introns = set()
for l in open(sys.argv[2]):
    if l.strip("EN"):
        continue
    t = l.strip().split("-")
    supported_introns.add((int(t[0]), int(t[1])))

read_intron_dict = {}
confirmed_reads_introns = {}
for l in open(sys.argv[1]):
    t = l.strip().split()
    if t[2] != "inconsistent":
        continue
    read_id = t[0]
    intron_count = t[4].count(",")
    if intron_count not in read_intron_dict:
        read_intron_dict[intron_count] = defaultdict(int)
    if intron_count not in confirmed_reads_introns:
        confirmed_reads_introns[intron_count] = defaultdict(int)

    for mt in match_types:
        event_info = t[3]
        event_pos = event_info.find(mt)
        if event_pos != -1:
            read_intron_dict[intron_count][mt] += 1
            next_event = event_info.find(",", event_pos)
            if next_event == -1:
                intron_str = event_info[event_pos+len(mt)+1:]
            else:
                intron_str = event_info[event_pos + len(mt) + 1:next_event]
            coords = intron_str.split("-")
            intron = (int(coords[0]), int(coords[1]))
            if intron in supported_introns:
                confirmed_reads_introns[intron_count][mt] += 1

    read_intron_dict[intron_count]["total_reads"] += 1
    confirmed_reads_introns[intron_count]["total_reads"] += 1

print_stats(read_intron_dict)
print_stats(confirmed_reads_introns)




