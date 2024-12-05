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
import pyfaidx


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


def load_genes(inf):
    gene_dict = {}
    for l in open(inf):
        v = l.split()
        gene_id = v[0]
        pa1 = v[6]
        pa2 = v[7]
        if pa2 == "NA": continue

        gene_dict[gene_id] = {pa1, pa2}
    return gene_dict


# 00000cbd-aa5d-4bf5-b193-3cead774ea12    ENSG00000240972.2       OldL5-6 GTGGGTCAGACTTT  ATGGAAGCA       ;%;chr22_23894583_23894771_+;%;chr22_23894945_23895039_+        chr22_23894402_23894402_+       chr22_23895224_23895224_+       ;%;chr22_23894402_23894582_+;%;chr22_23894772_23894944_+;%;chr22_23895040_23895223_+       known   2       ENST00000215754.8       protein_coding  chr22_23895103_23895103_+       121     ;%;chr22_23895104_23895224_+
def process_allinfo_with_utr(allinfo_utr, gene_polya_dict):
    gene_utr_counts = defaultdict(int)
    for l in open(allinfo_utr):
        v = l.split()
        gene_id = v[1]
        if gene_id not in gene_polya_dict: continue
        polya = v[7]
        if polya not in gene_polya_dict[gene_id]: continue
        group = "Old" if v[2].startswith("Old") else "Young"
        utr = v[-1]
        gene_utr_counts[(gene_id, utr, group)] += 1

    return gene_utr_counts


def get_utr_classification(gene_utr_counts, chr_dict):
    gene_dict = defaultdict(list)
    for gene_id, utr, group in gene_utr_counts.keys():
        gene_dict[gene_id].append((utr, group))


    for gene_id in gene_dict.keys():
        print("Processing gene %s" % gene_id)
        utrs = set([x[0] for x in gene_dict[gene_id]])
        print(utrs)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file prefix", required=True)
    parser.add_argument("--allinfo", "-a", type=str, help="input ALLINFO with UTRs", required=True)
    parser.add_argument("--gene_list", "-l", type=str, help="gene list with polyAs and FDRs", required=True)
    parser.add_argument("--genome", "-g", type=str, help="genome in FASTA", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gene_polya_dict = load_genes(args.gene_list)
    gene_utr_counts = process_allinfo_with_utr(args.allinfo, gene_polya_dict)
    get_utr_classification(gene_utr_counts, None)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
