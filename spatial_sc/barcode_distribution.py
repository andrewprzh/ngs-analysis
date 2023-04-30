#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
from numpy import histogram


def barcode_distribution(inf):
    barcode_counts = defaultdict(int)
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) < 10 or v[6] == "*":
            continue
        barcode_counts[v[6]] += 1
    return barcode_counts


def create_whitelist(barcode_counts, min_count, max_count, out_filename):
    with open(out_filename, "w") as outf:
        for b in barcode_counts:
            if min_count > barcode_counts[b]:
                continue
            if max_count > 0 and barcode_counts[b] > max_count:
                continue
            outf.write(b + "\n")


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name")
    parser.add_argument("--min_count", type=int, help="minimal barcode count", default=0)
    parser.add_argument("--max_count", type=int, help="maximal barcode count", default=0)
    parser.add_argument("--input", "-i", type=str, help="input barcodes", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not args.output:
        file_full_name = os.path.basename(args.input)
        fname, outer_ext = os.path.splitext(file_full_name)
        args.output = fname + ".filtered_barcodes" + outer_ext
    barcode_counts = barcode_distribution(args.input)
    max_count = max(barcode_counts.values())
    bin_size = 100
    bins = [i * 10 for i in range(21)] +  [i * bin_size for i in range(3, int(max_count / bin_size) + 2)]
    hist, b = histogram(list(barcode_counts.values()), bins)
    for i, v in enumerate(hist):
        print("%d\t%d" % (b[i], hist[i]))

    create_whitelist(barcode_counts, args.min_count, args.max_count, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
