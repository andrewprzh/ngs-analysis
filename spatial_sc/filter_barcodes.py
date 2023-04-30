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


def filter_barcodes(inf, outf, filter_umis=False, min_score=None, whitelist=None):
    with open(outf, "w") as filtered_barcodes:
        for l in open(inf):
            v = l.strip().split('\t')
            if len(v) < 10 or v[6] == "*":
                continue
            if filter_umis and v[9] != "True":
                continue
            if min_score is not None and int(v[8]) < min_score:
                continue
            if whitelist and v[6] not in whitelist:
                continue
            filtered_barcodes.write(l)


def load_barcode_set(inf):
    barcodes = set()
    for l in open(inf):
        barcodes.add(l.strip())
    return barcodes


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name")
    parser.add_argument("--trusted_umis", action="store_true", help="keep only trusted UMIs", default=False)
    parser.add_argument("--min_score", "-m", type=int, help="score to filter (default: do not filter)")
    parser.add_argument("--input", "-i", type=str, help="input barcodes", required=True)
    parser.add_argument("--white_list", type=str, help="barcode whitelist")


    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not args.output:
        file_full_name = os.path.basename(args.input)
        fname, outer_ext = os.path.splitext(file_full_name)
        args.output = fname + ".filtered" + outer_ext
    whitelist = set() if not args.white_list else load_barcode_set(args.white_list)
    filter_barcodes(args.input, args.output, args.trusted_umis, args.min_score, whitelist=whitelist)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
