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

MISSED_BARCODES = 0
CONTRADICTORY_CLUSTERS = 0

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
    return barcode_map


def process_reads(args, barcode_map):
    MISSED_BARCODES = 0
    outf = open(args.output, "w")
    for l in open(args.read_to_barcode):
        tokens = l.strip().split()
        read_id = tokens[args.readid_column]
        barcode = "none"
        for i in args.read_barcode_columns:
            barcode = tokens[i]
            if barcode != "none":
                break

        if barcode not in barcode_map:
            MISSED_BARCODES += 1
        else:
            outf.write(read_id + "\t" + barcode + "\t" + barcode_map[barcode] + "\n")
    outf.close()


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
    print("Contradictory barcode assignments " + str(CONTRADICTORY_CLUSTERS))
    print("Reads with non-informative barcodes " + str(MISSED_BARCODES))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)