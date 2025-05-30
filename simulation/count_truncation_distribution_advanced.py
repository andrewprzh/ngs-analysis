import math
import sys

import numpy
import pysam
from traceback import print_exc
import numpy as np
from collections import defaultdict

BIN_COUNT = 100
SMALL_BIN_COUNT = 10
HIST_STEP = 1.0 / BIN_COUNT


class TruncationStats:
    def __init__(self):
        self.right_truncations = [0 for _ in range(BIN_COUNT)]
        self.left_truncations = []
        for _ in range(BIN_COUNT):
            self.left_truncations.append([0 for _ in range(SMALL_BIN_COUNT)])

    def add(self, left_truncation, right_truncation):
        right_truncation_index = int(right_truncation * BIN_COUNT)
        if right_truncation_index == BIN_COUNT:
            right_truncation_index -= 1
        self.right_truncations[right_truncation_index] += 1

        left_truncation_index = int(left_truncation * SMALL_BIN_COUNT)
        if left_truncation_index == SMALL_BIN_COUNT:
            left_truncation_index -= 1
        self.left_truncations[right_truncation_index][left_truncation_index] += 1

    def normalize(self):
        s = sum(self.right_truncations)
        if s > 0:
            self.right_truncations = [float(x) / float(s) for x in self.right_truncations]
        for i in range(BIN_COUNT):
            s = sum(self.left_truncations[i])
            if s > 0:
                self.left_truncations[i] = [float(x) / float(s) for x in self.left_truncations[i]]


class TlenTruncationStats:
    TRANSCRIPT_LEN_BINS = list(sorted([500, 1000, 1500, 2000, 3000, 5000, 10000, math.inf]))
    def __init__(self):
        self.transcript_len_dict = {}
        for tlen in self.TRANSCRIPT_LEN_BINS:
            self.transcript_len_dict[tlen] = TruncationStats()
        self.read_lengths = []

    def add(self, transcript_len, reference_start, reference_end):
        left_truncation = float(reference_start) / float(transcript_len)
        right_truncation = float(transcript_len - reference_end) / float(transcript_len)
        for tlen in self.TRANSCRIPT_LEN_BINS:
            if transcript_len <= tlen:
                self.transcript_len_dict[tlen].add(left_truncation, right_truncation)
        self.read_lengths.append(reference_end - reference_start + 1)

    def print(self):
        print(numpy.mean(self.read_lengths))
        print(numpy.median(self.read_lengths))
        print()
        sys.stdout.write("{")
        for tlen in self.TRANSCRIPT_LEN_BINS:
            if tlen == math.inf:
                sys.stdout.write("math.inf: (")
            else:
                sys.stdout.write("%d: (" % tlen)
            self.transcript_len_dict[tlen].normalize()
            print(self.transcript_len_dict[tlen].right_truncations, ", ")
            print(self.transcript_len_dict[tlen].left_truncations, "),")
        print('}')


def main():
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
    truncation_stats = TlenTruncationStats()
    assigned_transcripts = defaultdict(int)
    t_lens = {}

    for alignment in bamfile:
        if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
              continue

        transcript_len = bamfile.get_reference_length(alignment.reference_name)
        truncation_stats.add(transcript_len, alignment.reference_start, alignment.reference_end)
        assigned_transcripts[alignment.reference_name] += 1
        t_lens[alignment.reference_name] = transcript_len

    truncation_stats.print()

    for t in t_lens.keys():
        if t_lens[t] > 8000:
            assigned_transcripts[t] = 0

    scale_factor = sum(assigned_transcripts.values()) / 1000000.0
    with open(sys.argv[2], "w") as outf:
        for tid in sorted(assigned_transcripts.keys()):
            outf.write("%s\t%.2f\t%.6f\n" % (tid, assigned_transcripts[tid], assigned_transcripts[tid] / scale_factor))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
