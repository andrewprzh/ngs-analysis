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

import numpy


def sort_dict(d):
    res = []
    for k in d.keys():
        res.append((k, d[k]))
    res = sorted(res, key=lambda x: x[1])
    return res


def dump_pairs(outf, lp):
    for k, v in lp:
        outf.write("%s\t%.6f\n" % (k, v))


def load_counts(inf):
    counts = defaultdict(tuple)

    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 3: continue
        counts[v[0]] = (int(v[1]), int(v[2]))
    return counts


def count_fold_change(count_dicts):
    fold_changes = defaultdict(list)
    for d in count_dicts:
        for k in d.keys():
            ont_count = d[k][0]
            sr_count = d[k][1]
            if sr_count == 0:
                fold_changes[k].append(-1)
            else:
                fold_changes[k].append(float(ont_count / sr_count))

    sr_dominated = defaultdict(int)
    lr_dominated = defaultdict(int)
    mean_fold_changes = defaultdict(float)
    for k in fold_changes.keys():
        sr_dominated[k] = fold_changes[k].count(0)
        lr_dominated[k] = fold_changes[k].count(-1)
        fold_change_vals = filter(lambda x: x > 0, fold_changes[k])
        if len(fold_change_vals) > 3:
            mean_fold_changes[k] = numpy.mean(fold_change_vals)

    return mean_fold_changes, lr_dominated, sr_dominated


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

    print("Merging counts")
    mean_fold_changes, lr_dominated, sr_dominated = count_fold_change(count_dicts)

    fold_change_val = sort_dict(mean_fold_changes)
    lr_dominated_val = sort_dict(lr_dominated)
    sr_dominated_val = sort_dict(sr_dominated)

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Fold change values\n")
        dump_pairs(outf, fold_change_val)
        outf.write("\nMost represented in SR\n")
        dump_pairs(outf, sr_dominated_val)
        outf.write("\nMost represented in LR\n")
        dump_pairs(outf, lr_dominated_val)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)