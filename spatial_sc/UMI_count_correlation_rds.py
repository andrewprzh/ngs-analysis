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
# LH00376:28:22GVKLLT3:4:1101:2655:1016   41      -1      8       25      +       CAACTGGCCGGGTA  TACCGTCGT       13      False


def get_umi_dict(inf, trusted_only=False):
    umi_dict = defaultdict(set)
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split('\t')
        if v[-1] in ['True', 'False']:
            # old format
            if v[6] == '*' or v[7] == '*': continue
            if trusted_only and v[-1] != "True": continue
            umi_dict[v[6]].add(v[7])
            continue

        if v[1] == '*' or v[2] == '*': continue
        if trusted_only and v[4] != "True": continue
        umi_dict[v[1]].add(v[2])

    count_dict = defaultdict(int)
    for k in umi_dict.keys():
        count_dict[k] = len(umi_dict[k])
    return count_dict


def load_counts(inf):
    barcode_counts = defaultdict(int)

    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 2: continue
        barcode_counts[v[0]] = int(v[1])
    return barcode_counts


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with UMI count tables", required=True)
    parser.add_argument("--lr_barcodes", type=str, help="input TSV with barcoded reads", required=True)
    parser.add_argument("--sr_counts", type=str, help="short read counts from RDS", required=True)
    parser.add_argument("--only_trusted_umi", default=False, action="store_true", help="keep only spliced reads")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    count_dicts = []

    print("Loading %s" % args.lr_barcodes)
    count_dicts.append(get_umi_dict(args.lr_barcodes, args.only_trusted_umi))

    print("Loading %s" % args.sr_counts)
    count_dicts.append(load_counts(args.sr_counts))

    print("Merging barcodes")
    all_barcodes = set()
    for d in count_dicts:
        all_barcodes.update(d.keys())

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Barcode\t" + "\t".join([args.lr_barcodes, args.sr_counts]) + "\n")
        for b in all_barcodes:
            counts = []
            for d in count_dicts:
                counts.append(d[b])
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
