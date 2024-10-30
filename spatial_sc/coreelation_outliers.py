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
import glob


def load_counts(inf):
    counts = defaultdict(int)

    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 2: continue
        counts[v[0]] = int(v[1])
    return counts


def fold_change(count_dicts):
    fold_changes = defaultdict(list)
    for d in count_dicts:
        for k in d.keys():
            fold_changes.


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with stats", required=True)
    parser.add_argument("--counts", nargs='+', type=str, help="input TSV with gene/barcodes counts LR vs SR", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    count_dicts = []

    for f in args.counts:
        for inf in glob.glob(f):
            print("Loading %s" % inf)
            count_dicts.append(load_counts(inf))

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
