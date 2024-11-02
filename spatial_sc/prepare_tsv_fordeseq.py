#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import os.path
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import glob

import numpy


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


def extract_counts(count_dicts):
    sr_counts = defaultdict(list)
    lr_counts = defaultdict(list)
    all_genes = set()
    for d in count_dicts:
        for k in d.keys():
            all_genes.add(k)
    for d in count_dicts:
        for k in all_genes:
            ont_count = d[k][0]
            sr_count = d[k][1]
            sr_counts[k].append(sr_count)
            lr_counts[k].append(ont_count)

    return all_genes, sr_counts, lr_counts


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with large TSV", required=True)
    parser.add_argument("--counts", "-c", nargs='+', type=str, help="input TSV with gene/barcodes counts LR vs SR", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    count_dicts = []

    samples_labels = []
    for f in args.counts:
        for inf in glob.glob(f):
            print("Loading %s" % inf)
            count_dicts.append(load_counts(inf))
            samples_labels.append(os.path.splitext(os.path.basename(inf))[0])

    print("Merging counts")
    all_genes, sr_counts, lr_counts = extract_counts(count_dicts)

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("gene_id\t" + "\t".join(map(lambda x: "SR_" + x, samples_labels)) +
                   "\t".join(map(lambda x: "LR_" + x, samples_labels)) + "\n")
        for g in all_genes:
            outf.write("\t".join(map(str, sr_counts[g])))
            outf.write("\n")
            outf.write("\t".join(map(str, lr_counts[g])))
            outf.write("\n")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
