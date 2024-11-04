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
        if isinstance(d[k], tuple):
            res.append((k, d[k][0], d[k][1]))
        else:
            res.append((k, d[k]))
    res = sorted(res, key=lambda x: x[1])
    return res


def sort_dict2(d):
    res = []
    for k in d.keys():
        if isinstance(d[k], tuple):
            res.append((k, d[k][0], d[k][1]))
        else:
            res.append((k, d[k]))
    res = sorted(res, key=lambda x: (x[1],x[2]), reverse=True)
    return res


def dump_pairs(outf, lp):
    for k, v in lp:
        outf.write("%s\t%.6f\n" % (k, v))

def dump_triples(outf, lp):
    for k, v1, v2 in lp:
        outf.write("%s\t%d\t%.6f\n" % (k, v1, v2))


def load_counts(inf):
    counts = defaultdict(tuple)

    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 3: continue
        try:
            counts[v[0]] = (int(v[1]), int(v[2]))
        except ValueError:
            pass
    return counts


def count_fold_change(count_dicts):
    fold_changes = defaultdict(list)
    sr_counts = defaultdict(list)
    lr_counts = defaultdict(list)
    for d in count_dicts:
        for k in d.keys():
            ont_count = d[k][0]
            sr_count = d[k][1]
            sr_counts[k].append(sr_count)
            lr_counts[k].append(ont_count)
            if sr_count == 0:
                fold_changes[k].append(-1)
            else:
                fold_changes[k].append(float(ont_count / sr_count))

    sr_dominated = {}
    lr_dominated = {}
    mean_fold_changes = defaultdict(float)
    for k in fold_changes.keys():
        sr_present = fold_changes[k].count(0)
        if sr_present >= 3:
            sr_dominated[k] = (sr_present, numpy.mean(sr_counts[k]))
        lr_present = fold_changes[k].count(-1)
        if lr_present >= 3:
            lr_dominated[k] = (lr_present, numpy.mean(lr_counts[k]))
        fold_change_vals = list(filter(lambda x: x > 0, fold_changes[k]))
        if len(fold_change_vals) > 3:
            fc = numpy.mean(fold_change_vals)
            if fc > 16 or fc < (1/16):
                mean_fold_changes[k] = fc

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
    lr_dominated_val = sort_dict2(lr_dominated)
    sr_dominated_val = sort_dict2(sr_dominated)

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Fold change values\t%d\t%d\n" % (len(list(filter(lambda x: x[1] < 1, fold_change_val))), len(list(filter(lambda x: x[1] > 1, fold_change_val)))))
        dump_pairs(outf, fold_change_val)
        outf.write("\nMost represented in SR\t%d\n" % len(sr_dominated_val))
        dump_triples(outf, sr_dominated_val)
        outf.write("\nMost represented in LR\t%d\n" % len(lr_dominated_val))
        dump_triples(outf, lr_dominated_val)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
