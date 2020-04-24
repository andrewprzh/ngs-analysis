import os
import sys
import pysam
import numpy

def find_starts(bam_filename):
    bam = pysam.AlignmentFile(bam_filename, "rb")
    starts = []
    for alignment in bam.fetch():
        if alignment.reference_id == -1 or alignment.is_secondary:
            continue
        blocks = sorted(alignment.get_blocks())
        if len(blocks) == 0:
            continue

        start = blocks[0][0]
        starting_cigar = alignment.cigartuples[0]
        if starting_cigar[0] == 4:
            start -= starting_cigar[1]
        starts.append(start)
    return starts


def find_ends(bam_filename):
    bam = pysam.AlignmentFile(bam_filename, "rb")
    ends = []
    for alignment in bam.fetch():
        if alignment.reference_id == -1 or alignment.is_secondary:
            continue
        blocks = sorted(alignment.get_blocks())
        if len(blocks) == 0:
            continue

        end = blocks[-1][1]
        ending_cigar = alignment.cigartuples[-1]
        if ending_cigar[0] == 4:
            end += ending_cigar[1]
            ends.append(end)
    return ends


def print_hist(coords, flag):
    hist = {}
    for v in coords:
        if v not in hist:
            hist[v] = 0
        hist[v] += 1
    positions = sorted(hist.keys())
    if not flag.startswith("e"):
        interval = range(max(-50, positions[0]), min(101, positions[-1] + 1))
    else:
        interval = range(max(37000, positions[0]), min(37401, positions[-1] + 1))
    for k in interval:
        v = 0 if k not in hist else hist[k]
        print(str(k) + "\t" + str(v))


if sys.argv[2].startswith("e"):
    coords = find_ends(sys.argv[1])
else:
    coords = find_starts(sys.argv[1])

print_hist(coords,  sys.argv[2])
