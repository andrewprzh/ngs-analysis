#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2022 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import gzip
import numpy
import pysam
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name", default="quality_stats.tsv")
    parser.add_argument("--input", "-i", type=str, help="input reads, fastq/fastq.gz/sam/bam", required=True)
    parser.add_argument("--mapped_only", action="store_true", help="use only mapped reads (effective for sam/bam only)", default=False)

    args = parser.parse_args()
    return args


def read_from_bam(inf, mapped_only=False):
    quality_dict = {}
    for a in pysam.AlignmentFile(inf, "r"):
        if a.is_unmapped and mapped_only:
            continue
        read_id = a.query_name
        if read_id in quality_dict or a.query_qualities is None:
            continue
        mean_q = float(sum(a.query_qualities)) / float(len(a.query_qualities))
        quality_dict[read_id] = mean_q

    return quality_dict


def read_from_fastq(inf):
    quality_dict = {}
    ref_name, outer_ext = os.path.splitext(os.path.basename(inf))
    low_ext = outer_ext.lower()
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(inf, "rt")
    else:
        handle = open(inf, "r")
    for r in SeqIO.parse(handle, "fastq"):
        read_id = r.id
        qual_array = r.letter_annotations['phred_quality']

        mean_q = float(sum(qual_array)) / float(len(qual_array))
        quality_dict[read_id] = mean_q

    return quality_dict


def count_stats(outf, quality_dict):
    q_values = list(quality_dict.values())
    mean_q = float(sum(q_values)) / float(len(q_values))
    max_q = int(max(q_values))
    bins = [i for i in range(max_q + 3)]
    q_hist, q_bins = numpy.histogram(q_values, bins)
    print("Average quality %.6f" % mean_q)
    print("Phred\tCount")
    for i in range(len(q_hist)):
        print("%d\t%d" % (q_bins[i], q_hist[i]))

    with open(outf, "w") as handle:
        for k in quality_dict.keys():
            handle.write("%s\t%.4f\n" % (k, quality_dict[k]))


def main():
    args = parse_args()
    ref_name, outer_ext = os.path.splitext(os.path.basename(args.input))
    low_ext = outer_ext.lower()
    if low_ext in ['.gz', '.gzip', '.fastq', '.fq']:
        quality_dict = read_from_fastq(args.input)
    elif low_ext in ['.sam', '.bam']:
        quality_dict = read_from_bam(args.input, args.mapped_only)
    else:
        print("Unknown format")
        exit(-1)

    count_stats(args.output, quality_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



