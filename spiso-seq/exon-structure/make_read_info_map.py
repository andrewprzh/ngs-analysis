############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import copy
import argparse
from traceback import print_exc
from functools import partial


def read_barcode_clusters(args):
    CONTRADICTORY_CLUSTERS = 0
    barcode_map = {}
    for l in open(args.barcode_to_property):
        tokens = l.strip().split()
        barcode = tokens[args.barcode_column]
        property = tokens[args.cluster_column]
        if barcode in barcode_map and barcode_map[barcode] != property:
            CONTRADICTORY_CLUSTERS += 1

        barcode_map[barcode] = property
    print("Total barcodes read " + str(len(barcode_map)))
    print("Contradictory barcode assignments " + str(CONTRADICTORY_CLUSTERS))
    return barcode_map


def process_reads(args, barcode_map):
    MISSED_BARCODES = 0
    outf = open(args.output, "w")
    for l in open(args.read_to_barcode):
        tokens = l.strip().split()
        read_id = tokens[args.readid_column]
        barcodes = []
        for i in args.read_barcode_columns:
            if tokens[i] != "none":
                barcodes.append((tokens[i], int(tokens[i + 1])))

        if len(barcodes) != 1:
            continue
        barcode = barcodes[0][0]
        kmers = barcodes[0][1]
        if barcode.find(',') != -1 or kmers < 3:
            continue

        if barcode not in barcode_map:
            MISSED_BARCODES += 1
        else:
            outf.write(read_id + "\t" + barcode + "\t" + barcode_map[barcode] + "\n")
    outf.close()
    print("Reads with non-informative barcodes " + str(MISSED_BARCODES))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--barcode_to_property", help="file with barcodes and cluster info", type=str)
    parser.add_argument("--barcode_column", help="0-based barcode column index", type=int, default=0)
    parser.add_argument("--cluster_column", help="0-based cluster property column index", type=int, default=1)

    parser.add_argument("--read_to_barcode", help="file with reads and assigned barcodes ", type=str)
    parser.add_argument("--readid_column", help="0-based read id column index", type=int, default=0)
    parser.add_argument("--read_barcode_columns", nargs='+', help="0-based barcode is column indices "
                                                                  "(will select first which is not none", type=int, default=1)

    parser.add_argument("--output", "-o", help="output prefix", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print("Reading clusters from " + args.barcode_to_property)
    barcode_map = read_barcode_clusters(args)
    print("Reading reads from " + args.read_to_barcode)
    process_reads(args, barcode_map)
    print("Read info map is written to " + args.output)


main()
