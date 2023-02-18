#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import sys
import pysam
import argparse
from traceback import print_exc
from collections import defaultdict


COVERAGE_PERCENT = 0.95
MAX_OVERLAP = 0.2


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", "-b", type=str, help="SAM/BAM file to assess", required=True)
    args = parser.parse_args()
    return args


def count_coverage(samf):
    cov90 = set()
    match90 = set()
    unaligned = set()
    all_contigs = set()

    for a in samf:
        all_contigs.add(a.query_name)
        if a.is_unmapped:
            unaligned.add(a.query_name)
            continue

        ctg_len = a.query_length
        ctg_aligned = a.query_alignment_length
        ref_len = samf.get_reference_length(a.reference_name)
        ref_aligned = a.reference_length

        if ref_aligned >=  ref_len * COVERAGE_PERCENT:
            cov90.add(a.reference_name)
            if ctg_aligned >= ctg_len * COVERAGE_PERCENT:
                match90.add(a.reference_name)

    print("Total: ", len(all_contigs))
    print("Unmapped: ", len(unaligned))
    print("Genes covered: ", len(cov90))
    print("Genes covered no extra", len(match90))


def interval_len(interval):
    return max(0, interval[1] - interval[0] + 1)


def overlap_intervals(range1, range2):
    return (max(range1[0], range2[0]), min(range1[1], range2[1]))


def overlap_fraction(range1, range2):
    return interval_len(overlap_intervals(range1, range2)) / min(interval_len(range1), interval_len(range2))


def is_misassembled(alignment_list):
    primary_alignment = None
    for a in alignment_list:
        if a.is_secondary or a.is_supplementary:
            continue
        primary_alignment = a
        break

    if not primary_alignment:
        return False

    for a in alignment_list:
        if not a.is_supplementary:
            continue
        if primary_alignment.reference_name != a.reference_name and \
                overlap_fraction((primary_alignment.query_alignment_start, primary_alignment.query_alignment_end),
                                 (a.query_alignment_start, a.query_alignment_end)) <= MAX_OVERLAP:
            if min(abs(primary_alignment.query_alignment_start - primary_alignment.query_alignment_end), abs(a.query_alignment_start - a.query_alignment_end)) < 200:
                return False
            print("Misassemble: %s" % a.query_name)
            print("%s: %d - %d" % (primary_alignment.reference_name, primary_alignment.query_alignment_start, primary_alignment.query_alignment_end))
            print("%s: %d - %d" % (a.reference_name, a.query_alignment_start, a.query_alignment_end))
            return True

    return False


def count_misassemblies(samf):
    alignment_dict = defaultdict(list)

    for a in samf:
        if a.is_unmapped:
            continue
        alignment_dict[a.query_name].append(a)

    candidates = 0
    misassemblies = 0
    for contig_name in alignment_dict.keys():
        if len(alignment_dict[contig_name]) == 1:
            continue
        candidates += 1
#
        if is_misassembled(alignment_dict[contig_name]):
            misassemblies += 1
            # print("Misassembly: %s: %s" % (contig_name, str([a.reference_name for a in alignment_dict[contig_name]])))

    print("Candidates: %d, misassemblies %d" % (candidates, misassemblies))

def main():
    args = parse_args()
    samf = pysam.AlignmentFile(args.bam, "r")
    count_coverage(samf)
    samf = pysam.AlignmentFile(args.bam, "r")
    count_misassemblies(samf)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


