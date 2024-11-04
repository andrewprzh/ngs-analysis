#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import gffutils
import numpy


# Intron length distribution (internal v. not interal)
# terminal exon lengths
# introns w GTAG, etc.
# DS with 100% data, 90%, 80%, and so on


def load_genes(inf):
    gene_list = []
    for l in open(inf):
        v = l.strip().split()
        gene_name = v[0]
        if v[0].find(".") != -1:
            continue
        gene_list.append(gene_name)

    return gene_list


def load_genedb(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    gene_dict = {}
    print("Loading genedb from %s" % genedb)
    for g in gffutils_db.features_of_type('gene'):
        if "gene_name" in g.attributes:
            gene_name = g.attributes["gene_name"][0]
        else:
            continue
        if "gene_type" in g.attributes:
            gene_type = g.attributes["gene_type"][0]
        else:
            continue
        exon_counts = []
        tlens = []
        for t in gffutils_db.children(g, featuretype=("transcript", "mRNA")):
            ex = 0
            tlen = 0
            for e in gffutils_db.children(t, featuretype="exon"):
                ex += 1
                tlen += e.end - e.start + 1
            exon_counts.append(ex)
            tlens.append(tlen)
        gene_dict[gene_name] = (gene_type, len(exon_counts), numpy.mean(exon_counts), max(tlens))
    return gene_dict


def count_stats(gene_list, gene_dict):
    total = 0
    unspliced = 0
    gene_lengths = []
    spliced_gene_lengths = []
    pcg = []

    gene_types = defaultdict(int)

    for g in gene_list:
        if g not in gene_dict: continue
        total += 1
        gene_info = gene_dict[g]
        gene_types[gene_info[0]] += 1
        if gene_info[0] == 'protein_coding':
            pcg.append(g)
        gene_lengths.append(gene_info[3])
        if gene_info[2] == 1:
            unspliced += 1
        else:
            spliced_gene_lengths.append(gene_info[3])

    print("Total genes found: %d" % total)
    print("Of them unspliced: %d" % unspliced)
    print(gene_types)
    print("Mean gene len %.2f" % numpy.mean(gene_lengths))
    print("Mean spliced gene len %.2f" % numpy.mean(spliced_gene_lengths))
    for g in sorted(pcg):
        print(g)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--de_genes", "-d", nargs='+', type=str, help="2 list of genes (up/down regulated)", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    gene_dict = load_genedb(args.genedb)
    for f in args.de_genes:
        gene_list = load_genes(f)
        count_stats(gene_list, gene_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
