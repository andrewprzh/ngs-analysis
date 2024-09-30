#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import argparse
from traceback import print_exc
from collections import defaultdict

# 82718f26-b81f-4561-ac75-59b987941431    ACGCGTTTAAGACG  ATCGGCTAG       14      True    -       82      40      49      66
def get_umi_dict(inf, trusted_only=False):
    umi_dict = defaultdict(set)
    for l in inf:
        if l.startswith("#"): continue
        v = l.strip().split("\t")
        if v[1] == '*' or v[2] == '*': continue
        if trusted_only and v[4] != "True": continue
        umi_dict[v[1]].add(v[2])

    return umi_dict


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with UMI count tables", required=True)
    parser.add_argument("--input", "-i", nargs="+", type=str, help="input TSV(s) with barcoded reads", required=True)
    parser.add_argument("--only_trusted_umi", default=False, action="store_true", help="keep only spliced reads")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    dicts = []
    for tsv in args.input:
        print("Loading %s" % tsv)
        dicts.append(get_umi_dict(tsv, args.only_trusted_umi))

    print("Merging barcodes")
    all_barcodes = set()
    for d in dicts:
        all_barcodes.update(d.keys())

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Barcode\t" + "\t".join(args.input) + "\n")
        for b in all_barcodes:
            counts = []
            for d in dicts:
                counts.append(len(d[b]))
            outf.write(b + "\t" + "\t".join(map(str, counts)) + "\n")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
