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


def filter_assignments(inf, outf, whitelist):
    with open(outf, "w") as filtered_assignments:
        for l in open(inf):
            v = l.strip().split('\t')
            if v[0] not in whitelist:
                continue
            filtered_assignments.write(l)


def load_readid_set(inf):
    barcodes = set()
    for l in open(inf):
        barcodes.add(l.strip())
    return barcodes


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name")
    parser.add_argument("--input", "-i", type=str, help="input assignmetns / all info", required=True)
    parser.add_argument("--white_list", type=str, help="read whitelist")


    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not args.output:
        file_full_name = os.path.basename(args.input)
        fname, outer_ext = os.path.splitext(file_full_name)
        args.output = fname + ".filtered" + outer_ext
    whitelist = set() if not args.white_list else load_readid_set(args.white_list)
    if whitelist:
        print("Loaded %d read ids" % len(whitelist))
    filter_assignments(args.input, args.output, whitelist)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
