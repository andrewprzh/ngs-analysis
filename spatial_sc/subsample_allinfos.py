#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
import random
import time
from collections import defaultdict
import numpy


def load_allinfo(inf):
    print("Loading allinfo from %s" % inf)
    genes_buckets = defaultdict(list)
    count = 0
    for l in open(inf):
        v = l.split("\t")
        genes_buckets[v[1]].append(l)
        count += 1
    print("Loaded %d reads and %d genes" % (count, len(genes_buckets)))
    return genes_buckets


def select_high_count_genes(gene_buckets, min_read_count):
    selected_genes = set()
    for gene_id in gene_buckets.keys():
        if len(gene_buckets[gene_id]) < min_read_count:
            continue
        selected_genes.add(gene_id)
    print("Selected %d genes with >= %d reads" % (len(selected_genes), min_read_count))
    return selected_genes


def subsample(gene_buckets1, gene_buckets2, common_genes, subsample_genes_count, subsample_read_count):
    selected_genes = numpy.random.choice(list(common_genes), subsample_genes_count)
    subsampled_buckets1 = {}
    subsampled_buckets2 = {}
    for gene_id in selected_genes:
        subsampled_buckets1[gene_id] = list(numpy.random.choice(gene_buckets1[gene_id], subsample_read_count))
        subsampled_buckets2[gene_id] = list(numpy.random.choice(gene_buckets2[gene_id], subsample_read_count))
    return subsampled_buckets1, subsampled_buckets2


def dump_allinfo(gene_buckets, iteration, output_prefix, output_suffix):
    outf = open(os.path.join(output_prefix, "%d_%s.allinfo" % (iteration, output_suffix)), "w")
    for gene_id in gene_buckets.keys():
        for l in gene_buckets[gene_id]:
            outf.write(l)
    outf.close()


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder name", required=True)
    parser.add_argument("--input", "-i", nargs=2, type=str, help="2 input allinfo files", required=True)
    parser.add_argument("--iterations", type=int, help="numer of iterations", default=100)
    parser.add_argument("--subsample_genes_count", "-g", type=int, help="number of genes to select at each iteration", default=100)
    parser.add_argument("--min_read_count", "-m", type=int, help="minimap number of reads per gene", default=40)
    parser.add_argument("--subsample_read_count", "-r", type=int, help="number of reads per gene to select at each iteration", default=20)
    parser.add_argument("--seed", "-s", type=int, help="default randomness seed (-1 for time seed)", default=11)

    args = parser.parse_args(argv)
    return args


def main(argv):
    args = parse_args(argv)
    if args.seed == -1:
        random.seed(time.time())
    else:
        random.seed(args.seed)
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    gene_buckets1 = load_allinfo(args.input[0])
    genes1 = select_high_count_genes(gene_buckets1, args.min_read_count)
    gene_buckets2 = load_allinfo(args.input[1])
    genes2 = select_high_count_genes(gene_buckets2, args.min_read_count)

    selected_genes = genes1.intersection(genes2)
    print("Total number of common genes: %d" % len(selected_genes))
    subsample_genes_count = args.subsample_genes_count
    if len(selected_genes) < args.subsample_genes_count * 2:
        subsample_genes_count = len(selected_genes) // 2
        print("The number will be reduced to %d" % subsample_genes_count)

    for i in range(args.iterations):
        subsampled_buckets1, subsampled_buckets2 = subsample(gene_buckets1, gene_buckets2, selected_genes,
                                                             subsample_genes_count, args.subsample_read_count)
        dump_allinfo(subsampled_buckets1, i, args.output, "1")
        dump_allinfo(subsampled_buckets2, i, args.output, "2")
    print("Done. Results are in %s" % args.output)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
