#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gffutils
from Bio import SeqIO
import gzip
from traceback import print_exc
import numpy


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output FASTQ file", default="ont_truncated.fastq")
    parser.add_argument("--all_info", "-a", type=str, help="all info file", required=True)
    parser.add_argument("--ont", "-i", type=str, help="input file with ONT sequences", required=True)
    parser.add_argument("--read_length", type=int, help="illumina read length [76]", default=76)
    parser.add_argument("--mean_frag", type=int, help="mean fragment length [300]", default=300)
    parser.add_argument("--stdev_frag", type=int, help="standard deviation of the fragment length [30]", default=30)

    args = parser.parse_args()
    return args


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def read_all_info(all_info_file):
    read_dict = {}
    for l in open(all_info_file):
        vals = l.strip().split()
        if len(vals) < 7:
            continue
        read_id = vals[0]
        if read_id[0] == '@':
            read_id = read_id[1:]

        if vals[2] != 'none':
            read_dict[read_id] = 1
        elif vals[5] != 'none':
            read_dict[read_id] = -1
    return read_dict


def truncate_reads(infastq, read_dict, outfastq_stream):
    pass



def main():
    args = parse_args()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
