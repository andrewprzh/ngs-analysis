#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
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

import numpy


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--guide_info", "-g", type=str, help="Table with guides / factors", required=True)
    parser.add_argument("--guide_columns", type=str, help="comma-separated guide columns, first one is the guide column", default="0,3,1")
    parser.add_argument("--barcode2guide", "-b", type=str, help="Table with barcodes and guides", required=True)
    parser.add_argument("--barcode_umi", "-u", type=str, help="barcode-UMI table", required=True)
    args = parser.parse_args()
    return args


def read_barcodes(inf_name):
    umi_dict = defaultdict(set)
    for l in open(inf_name):
        v = l.strip().split('\t')
        umi = v[3]
        bc = v[2]
        umi_dict[bc].add(umi)

    return umi_dict


def parse_barcode_table(guides, barcodes, columns):
    cols = list(map(int, columns.split(',')))
    guide_column = cols[0]
    print(cols)
    min_cols = max(cols)
    cols = cols[1:]
    guide_dict = {}
    for l in open(guides):
        if l.startswith("#"):
            continue
        v = l.strip().split('\t')
        if len(v) <= min_cols:
            continue
        guide_dict[v[guide_column]] = [v[col].replace(" ", "_") for col in cols]

    barcode_dict = {}
    for l in open(barcodes):
        v = l.strip().split('\t')
        barcode = v[0]
        guide = "_".join(v[1].split("_")[:3])
        group = [guide_dict[guide][1], guide, barcode]
        barcode_dict[barcode] = group

    return barcode_dict


def compute_counts(barcode_dict, umi_dict, col):
    umi_counts = defaultdict(int)
    bc_counts = defaultdict(int)
    for bc in umi_dict.keys():
        if bc not in barcode_dict:
            continue
        umi_counts[barcode_dict[bc][col]] += len(umi_dict[bc])
        bc_counts[barcode_dict[bc][col]] += 1
    return umi_counts, bc_counts


def print_dict(barcode_dict, output_file):
    outf = open(output_file, "w")
    for bc in sorted(barcode_dict.keys()):
        outf.write("%s\t%s\n" % (bc, barcode_dict[bc]))
    outf.close()


def save_hist(count_dict, bin_len, bin_count, out_fname):
    umi_counts = list(count_dict.values())
    bins = [i * bin_len for i in range(bin_count)]
    hist = numpy.histogram(umi_counts, bins=bins)
    with open(out_fname, "w") as outf:
        for i, count in enumerate(hist[0]):
            outf.write("%d\t%d\n" % ((hist[1][i] + hist[1][i + 1]) / 2, count))


def main():
    #set_logger(logger)
    args = parse_args()
    barcode_dict = parse_barcode_table(args.guide_info, args.barcode2guide, args.guide_columns)
    umi_dict = read_barcodes(args.barcode_umi)
    print(len(umi_dict), len(barcode_dict))
    bin_sizes = [1000, 1000, 100]
    for i in range(3):
        umi_counts, bc_counts = compute_counts(barcode_dict, umi_dict, i)
        save_hist(umi_counts, bin_sizes[i], 100, args.output + ".UMI_hist.%d" % i + ".tsv")
        save_hist(bc_counts, 2, 30, args.output + ".BC_hist.%d" % i + ".tsv")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



