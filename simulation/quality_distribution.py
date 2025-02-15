############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import pysam
import numpy
from pyfaidx import Fasta
from traceback import print_exc
from enum import Enum, unique
from collections import defaultdict
from scipy.signal import savgol_filter


@unique
class BaseAlignmentType(Enum):
    match = 0
    mismatch = 1
    insertion = 2
    clip = 3
    deletion = 4


def process_read(alignment, chr_str, quality_dict, err_dict):
    read_pos = 0
    ref_pos = alignment.reference_start
    try:
        qualities = alignment.query_qualities
    except ValueError:
        return

    if not qualities: return

    for cigartuple in alignment.cigartuples:
        event = cigartuple[0]
        event_len = cigartuple[1]

        if event == pysam.CHARD_CLIP or event == pysam.CPAD:
            continue
        elif event == pysam.CREF_SKIP:
            ref_pos += event_len
            continue
        elif event == pysam.CDEL:
            err_dict[BaseAlignmentType.deletion] += event_len
            ref_pos += event_len
            continue
        for i in range(event_len):
            if event == pysam.CSOFT_CLIP:
                quality_dict[BaseAlignmentType.clip][qualities[read_pos]] += 1
                read_pos += 1
            elif event == pysam.CDIFF:
                quality_dict[BaseAlignmentType.mismatch][qualities[read_pos]] += 1
                err_dict[BaseAlignmentType.mismatch] += 1
                read_pos += 1
                ref_pos += 1
            elif event == pysam.CEQUAL:
                quality_dict[BaseAlignmentType.match][qualities[read_pos]] += 1
                err_dict[BaseAlignmentType.match] += 1
                read_pos += 1
                ref_pos += 1
            elif event == pysam.CMATCH:
                base_type = BaseAlignmentType.mismatch if alignment.query_sequence[read_pos] != chr_str[ref_pos] else BaseAlignmentType.match
                err_dict[base_type] += 1
                quality_dict[base_type][qualities[read_pos]] += 1
                read_pos += 1
                ref_pos += 1
            elif event == pysam.CINS:
                quality_dict[BaseAlignmentType.insertion][qualities[read_pos]] += 1
                err_dict[BaseAlignmentType.insertion] += 1
                read_pos += 1


def process_bam(in_bam, ref_dict):
    quality_dict = defaultdict(lambda: defaultdict(int))
    err_dict = defaultdict(int)
    current_chr_id = ""
    current_chr_str = ""

    for alignment in in_bam:
        if alignment.is_secondary or alignment.is_supplementary or alignment.reference_id == -1:
            continue

        chr_id = alignment.reference_name
        if chr_id != current_chr_id:
            current_chr_id = chr_id
            current_chr_str = ref_dict[current_chr_id]

        process_read(alignment, current_chr_str, quality_dict, err_dict)
    return quality_dict, err_dict


def process_counts(count_dict: defaultdict):
    total_count = sum(count_dict.values())
    max_val = max(count_dict.keys())

    normalized_probabilities = [count_dict[i] / total_count for i in range(max_val)]
    smoothed_probabilities = savgol_filter(normalized_probabilities, 6, 3)

    for i in range(len(smoothed_probabilities)):
        if smoothed_probabilities[i] < 0:
            smoothed_probabilities[i] = 0
    total_smoothed_count = sum(smoothed_probabilities)
    normalized_smooth_probabilities = [smoothed_probabilities[i] / total_smoothed_count for i in range(max_val)]

    return normalized_probabilities, normalized_smooth_probabilities


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="reference FASTA", type=str, required=True)
    parser.add_argument("--bam", "-b", help="BAM file to analyse", type=str, required=True)
    parser.add_argument("--output", "-o", help="output file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    print("Loading reference genome from %s" % args.fasta)
    ref_dict = Fasta(args.fasta)

    print("Processing %s" % args.bam)
    quality_dict, err_dict = process_bam(pysam.AlignmentFile(args.bam), ref_dict)
    #print(quality_dict.keys())

    out_stream = sys.stdout
    if args.output:
        out_stream = open(args.output, "w")

    for alignment_base_type in quality_dict.keys():
        print("Type: %s, total observed events: %d" % (alignment_base_type, sum(quality_dict[alignment_base_type].values())))
        normalized_probabilities, smoothed_probabilities = process_counts(quality_dict[alignment_base_type])
        out_stream.write("%s\n" % str(alignment_base_type))
        out_stream.write("\t".join(map(lambda x: "%.4f" % x, normalized_probabilities)) + "\n")
        out_stream.write("\t".join(map(lambda x: "%.4f" % x, smoothed_probabilities)) + "\n")

    total_read_bases = (err_dict[BaseAlignmentType.match] + err_dict[BaseAlignmentType.mismatch] +
                        err_dict[BaseAlignmentType.insertion]  + err_dict[BaseAlignmentType.deletion])
    total_errors = err_dict[BaseAlignmentType.mismatch] + err_dict[BaseAlignmentType.insertion]  + err_dict[BaseAlignmentType.deletion]
    out_stream.write("Total error rate: %.2f\n" % (100 * total_errors / total_read_bases))
    out_stream.write("Mismatch rate: %.2f\n" % (100 * err_dict[BaseAlignmentType.mismatch] / total_read_bases))
    out_stream.write("Insertion rate: %.2f\n" % (100 * err_dict[BaseAlignmentType.insertion] / total_read_bases))
    out_stream.write("Deletion rate: %.2f\n" % (100 * err_dict[BaseAlignmentType.deletion] / total_read_bases))

    out_stream.close()

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
