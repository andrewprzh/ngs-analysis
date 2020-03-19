############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import copy
import argparse
import numpy
import pysam
from common import *
from traceback import print_exc
from functools import partial

POS_DIFF = 10000


def read_info(info_file):
    barcode_map1 = {}
    barcode_map2 = {}
    read_pairs = []
    for l in open(info_file):
        tokens = l.strip().split()
        if len(tokens) != 4:
            continue
        barcode_map1[tokens[2]] = (tokens[0], tokens[1])
        barcode_map2[tokens[3]] = (tokens[0], tokens[1])
        read_pairs.append((tokens[2], tokens[3]))
    return barcode_map1, barcode_map2, read_pairs


def read_bam_file(infile, barcode_map):
    alignment_map = {}
    for alignment in pysam.AlignmentFile(infile, "rb"):
        if alignment.reference_id == -1:
            continue
        read_id = alignment.query_name
        if read_id in barcode_map:
            if read_id not in alignment_map:
                alignment_map[read_id] = []
            alignment_map[read_id].append(alignment)

    return alignment_map


def compare_alignment_sets(alignments1, alignments2):
    if len(alignments1) == 0:
        return "first_empty"
    if len(alignments2) == 0:
        return "second_empty"

    for a1 in alignments1:
        for a2 in alignments2:
            if a1.reference_id == a2.reference_id:
                if abs(a1.reference_start - a2.reference_start) < POS_DIFF:
                    if a1.is_secondary or a1.is_supplementary or a2.is_secondary or a2.is_supplementary:
                        return "close_secondary"
                    else:
                        return "close_primary"
                else:
                    return "same_chr"
    return "diff"


def compare_reads(read_pairs, alignment_map1, alignment_map2, barcode_map1, barcode_map2):
    stats = {}
    for read_pair in read_pairs:
        read_id1 = read_pair[0]
        read_id2 = read_pair[1]
#        print(read_id1 + " " + read_id2)
        if barcode_map1[read_id1] != barcode_map2[read_id2]:
            print("Unequal barcode/UMI")
            continue
        if read_id2 not in alignment_map2:
#            print(read_id2 + " was not found")
            stats[barcode_map1[read_id1]] = "second_abset"
            continue
        if read_id1 not in alignment_map1:
#            print(read_id1 + " was not found")
            stats[barcode_map1[read_id1]] = "first_abset"
            continue
        stats[barcode_map1[read_id1]] = compare_alignment_sets(alignment_map1[read_id1], alignment_map2[read_id2])

    return stats


def print_stats(stats, out_file):
    aggr_stats = {}
    outf = open(out_file, "w")
    for k in stats.keys():
        v = stats[k]
        if v not in aggr_stats:
            aggr_stats[v] = 0
        aggr_stats[v] += 1
        outf.write(k[0] + "\t" + k[1] + "\t" + v + "\n")

    outf.close()
    for k in aggr_stats.keys():
        print(k + "\t" + str(aggr_stats[k]))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bams", nargs=2, help="bam files", type=str)
    parser.add_argument("--read_info",  help="file read ids with the same UMI and barcode but mapped to different genes ", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    args = parser.parse_args()

    if args.bams is None or args.read_info is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()
    print("Reading read ids")
    barcode_map1, barcode_map2, read_pairs = read_info(args.read_info)
    print("Collecting infor from first BAM file")
    alignment_map1 = read_bam_file(args.bams[0], barcode_map1)
    print("Collecting infor from second BAM file")
    alignment_map2 = read_bam_file(args.bams[1], barcode_map2)

    print("Counting stats")
    stats = compare_reads(read_pairs, alignment_map1, alignment_map2, barcode_map1, barcode_map2)
    print_stats(stats, args.output_prefix)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
