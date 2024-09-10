#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
import os
import sys
import pysam
from collections import defaultdict


def read_called_barcodes(inf):
    barcode_dict = {}
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) >= 2 and v[1] == "*":
            continue
        barcode_dict[v[0]] = v[1]
    print("Loaded %d read ids" % len(barcode_dict))
    return barcode_dict


def has_intron(cigartules):
    return any(x[0] == 3 for x in cigartules)


def count_stats(in_file_name, barcode_dict=None, intron_chr_set=None):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    total_reads = 0
    stat_dict = defaultdict(int)

    for read in inf:
        total_reads += 1

        if read.reference_id == -1:
            mapped = "unmapped"
        elif not read.is_secondary and not read.is_supplementary:
            mapped = "primary"
        else:
            mapped = "non-primary"

        if barcode_dict is not None and read.query_name in barcode_dict:
            barcoded = "barcoded"
        else:
            barcoded = "no barcode"

        if has_intron(read.cigartuples):
            spliced = "spliced"
        else:
            spliced = "unspliced"

        stat_dict[(mapped, barcoded, spliced)] += 1

        if total_reads % 100000 == 0:
            sys.stdout.write("Processed " + str(total_reads) + " reads\r")

    print("Processed %d reads" % total_reads)
    inf.close()
    for k in sorted(stat_dict.keys()):
        print("%s\t%d" % (k, stat_dict[k]))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--output", "-o", type=str, help="output folder (same as input file if not set), "
    #                                                      "or file name")
    parser.add_argument("--barcodes", type=str, help="read barcodes")
    parser.add_argument("--bam", "-b", nargs="+", type=str, help="BAM file ", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    barcode_dict = None
    if args.barcodes:
        barcode_dict = read_called_barcodes(args.barcodes)

    for in_bam in args.bam:
        print("Processing %s" % in_bam)

        count_stats(in_bam, barcode_dict)
        print("End processing %s" % in_bam)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
