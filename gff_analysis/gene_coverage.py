#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import subprocess
import sys
import argparse
from collections import defaultdict
import pysam
import gffutils
from common import *


def count_coverage(chr_id, start, end, bam):
    gene_len = end - start + 1
    coverage = [0 for _ in range(gene_len)]
    for a in bam.fetch(chr_id, start - 1, end):
        covered_start = max(a.reference_start + 1, start)
        covered_end = min(a.reference_end + 1, end)
        for pos in range(covered_start, covered_end + 1):
            coverage[pos - start] += 1
    return 1 - (coverage.count(0) / gene_len)


def get_gene_stats(gene_db, bam):
    gene_cov_dict = {}
    for g in gene_db.features_of_type('CDS', order_by=('seqid', 'start')):
        gene_name = g.id
        gene_cov = count_coverage(g.seqid, g.start, g.end, bam)
        gene_cov_dict[gene_name] = gene_cov
    print("Genes processed %d" % len(gene_cov_dict))
    return gene_cov_dict


def count_higher_than(container, value):
    count = 0
    for c in container:
        if c > value:
            count += 1
    return count


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", default="gtf_stats.tsv")
    parser.add_argument("--genedb", "-d", type=str, help="gene db", required=True)
    parser.add_argument("--bam", "-b", type=str, help="bam to count coverage", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    bam = pysam.AlignmentFile(args.bam, "rb")
    gene_cov_dict = get_gene_stats(gffutils_db, bam)
    cov_values = gene_cov_dict.values()
    for v in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]:
        print("> %.2f: %d" % (v, count_higher_than(cov_values, v)))

    with open(args.output, 'w') as outf:
        for k in sorted(gene_cov_dict.keys()):
            outf.write("%s\t%.3f\n" % (k, gene_cov_dict[k]))
    

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

