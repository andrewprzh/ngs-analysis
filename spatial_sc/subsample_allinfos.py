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
    # layer -> gene -> (sample, read)
    layer_backets = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    count = 0
    for l in open(inf):
        v = l.strip().split("\t")
        layer = v[2]
        gene_id = v[1]
        sample = v[11]
        layer_backets[layer][gene_id][sample].append("\t".join(v[:-1]))
        count += 1
    print("Loaded %d reads and %d layer-age-conditions" % (count, len(layer_backets)))
    return layer_backets


def select_high_count_genes(gene_buckets, min_read_count, max_from_single_sample, subsample_read_count):
    selected_genes = set()
    for gene_id in gene_buckets.keys():
        total_reads_per_gene = sum([len(gene_buckets[gene_id][sample]) for sample in gene_buckets[gene_id].keys()])
        if total_reads_per_gene < min_read_count:
            continue
        total_reads_per_gene_wlim = sum([min(max_from_single_sample, len(gene_buckets[gene_id][sample]))
                                         for sample in gene_buckets[gene_id].keys()])
        if total_reads_per_gene_wlim < subsample_read_count:
            continue
        selected_genes.add(gene_id)
    print("Selected %d genes with >= %d reads" % (len(selected_genes), min_read_count))
    return selected_genes


def pre_sample_reads(gene_bucket, max_from_single_sample):
    pre_samples_reads = []
    for sample in gene_bucket.keys():
        reads = gene_bucket[sample]
        if len(reads) > max_from_single_sample:
            pre_samples_reads += list(numpy.random.choice(reads, max_from_single_sample))
        else:
            pre_samples_reads += reads
    return pre_samples_reads


def subsample(gene_buckets1, gene_buckets2, common_genes, subsample_genes_count,
              subsample_read_count, max_from_single_sample):
    selected_genes = numpy.random.choice(list(common_genes), subsample_genes_count, replace=False)
    subsampled_buckets1 = {}
    subsampled_buckets2 = {}
    for gene_id in selected_genes:
        pre_sampled_reads1 = pre_sample_reads(gene_buckets1[gene_id], max_from_single_sample)
        pre_sampled_reads2 = pre_sample_reads(gene_buckets2[gene_id], max_from_single_sample)
        subsampled_buckets1[gene_id] = list(numpy.random.choice(pre_sampled_reads1, subsample_read_count, replace=False))
        subsampled_buckets2[gene_id] = list(numpy.random.choice(pre_sampled_reads2, subsample_read_count, replace=False))
    return subsampled_buckets1, subsampled_buckets2


def dump_allinfo(gene_buckets, out_fname):
    outf = open(out_fname, "a")
    for gene_id in gene_buckets.keys():
        for l in gene_buckets[gene_id]:
            outf.write(l + "\n")
    outf.close()


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder name", required=True)
    parser.add_argument("--input", "-i", type=str, help= "input allinfo file", required=True)
    parser.add_argument("--iterations", type=int, help="numer of iterations", default=100)
    parser.add_argument("--subsample_genes_count", "-g", type=int, help="number of genes to select at each iteration", default=200)
    parser.add_argument("--min_read_count", type=int, help="minimap number of reads per gene", default=75)
    parser.add_argument("--max_read_sample_count", type=int, help="minimap number of reads per gene per sample", default=17)
    parser.add_argument("--subsample_read_count", "-r", type=int, help="number of reads per gene to select at each iteration", default=50)
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

    layer_buckets = load_allinfo(args.input)
    print("Per-layer-age gene counts:")
    for layer in layer_buckets.keys():
        print("  %s: %d" % (layer, len(layer_buckets[layer])))

    high_count_layer_buckets = {}
    for layer in layer_buckets.keys():
        selected_genes = select_high_count_genes(layer_buckets[layer],
                                                 args.min_read_count,
                                                 args.max_read_sample_count,
                                                 args.subsample_read_count)
        high_count_layer_buckets[layer] = {}
        for g in selected_genes:
            high_count_layer_buckets[layer][g] = layer_buckets[layer][g]

    print("Per-layer-age gene counts after filtering:")
    for layer in high_count_layer_buckets.keys():
        print("  %s: %d" % (layer, len(high_count_layer_buckets[layer])))

    layer_no_age_set = set()
    for layer in high_count_layer_buckets.keys():
        if layer.startswith("Young"):
            layer_no_age_set.add(layer[5:])
        elif layer.startswith("Old"):
            layer_no_age_set.add(layer[3:])

    overlapping_genes = defaultdict(set)
    print("Total number of common genes:")
    for lna in layer_no_age_set:
        young_gene_set = set(high_count_layer_buckets["Young" + lna].keys())
        old_gene_set = set(high_count_layer_buckets["Old" + lna].keys())
        overlapping_genes[lna] = young_gene_set.intersection(old_gene_set)
        print("  %s: %d" % (lna, len(overlapping_genes[lna])))

    subsample_genes_count = args.subsample_genes_count
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    for i in range(args.iterations):
        print("Iteration %d" % (i + 1))
        outdir = os.path.join(args.output, "Run%d" % i)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outf = os.path.join(outdir, "sampled_AllInfo.tsv")
        open(outf, "w").close()
        for lna in layer_no_age_set:
            gene_buckets1 = high_count_layer_buckets["Young" + lna]
            gene_buckets2 = high_count_layer_buckets["Old" + lna]
            subsampled_buckets1, subsampled_buckets2 = subsample(gene_buckets1, gene_buckets2, overlapping_genes[lna],
                                                                 subsample_genes_count, args.subsample_read_count,
                                                                 args.max_read_sample_count)
            dump_allinfo(subsampled_buckets1, outf)
            dump_allinfo(subsampled_buckets2, outf)
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
