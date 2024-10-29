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
import gffutils

# 82718f26-b81f-4561-ac75-59b987941431    ACGCGTTTAAGACG  ATCGGCTAG       14      True    -       82      40      49      66
# LH00376:28:22GVKLLT3:4:1101:2655:1016   41      -1      8       25      +       CAACTGGCCGGGTA  TACCGTCGT       13      False


def load_allinfo(inf, gene_name_dict=None):
    print("Loading allinfo from %s" % inf)
    # read_id -> (chr, strand, gene, isoform)
    gene_count_dict = defaultdict(int)

    for l in open(inf):
        v = l.strip().split("\t")
        gene_id = v[1]
        if not gene_name_dict:
            gene_count_dict[gene_id] += 1
        else:
            if gene_id in gene_name_dict:
                gene_name = gene_name_dict[gene_id]
                gene_count_dict[gene_name] += 1

    print("Loaded %d genes" % len(gene_count_dict))
    return gene_count_dict


def load_counts(inf):
    barcode_counts = defaultdict(int)
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 2: continue
        barcode_counts[v[0]] = int(v[1])
    return barcode_counts


def load_gene_names(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    gene_names = {}
    print("Loading genedb from %s" % genedb)
    for g in gffutils_db.features_of_type(('gene')):
        if "gene_name" in g.attributes:
            gene_name = g.attributes["gene_name"][0]
        else:
            continue
        gene_names[g.id] = gene_name
    return gene_names


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with UMI count tables", required=True)
    parser.add_argument("--lr_allinfo", type=str, help="input TSV with ALLINFO", required=True)
    parser.add_argument("--sr_counts", type=str, help="short read counts from RDS", required=True)
    parser.add_argument("--genedb", type=str, help="genedb for gene names", required=True)
    parser.add_argument("--only_trusted_umi", default=False, action="store_true", help="keep only spliced reads")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    count_dicts = []

    gene_names = load_gene_names(args.genedb)

    print("Loading %s" % args.lr_allinfo)
    count_dicts.append(load_allinfo(args.lr_allinfo, gene_names))

    print("Loading %s" % args.sr_counts)
    count_dicts.append(load_counts(args.sr_counts))

    print("Merging barcodes")
    all_barcodes = set()
    for d in count_dicts:
        all_barcodes.update(d.keys())

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Barcode\t" + "\t".join([args.lr_allinfo, args.sr_counts]) + "\n")
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
