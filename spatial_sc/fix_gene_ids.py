#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
import gffutils


def load_transcripts(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    gene_ids = {}
    print("Loading genedb from %s" % genedb)
    for g in gffutils_db.features_of_type(('gene')):
        short_id = g.id.split(".")[0]
        gene_ids[short_id] = g.id
    return gene_ids


def fix_ids(gene_ids, gene_list, outf):
    with open(outf, "w") as output_file:
        for l in open(gene_list):
            gene_id = l.strip().split()[0]
            if gene_id not in gene_ids:
                continue
            output_file.write("%s\n" % gene_ids[gene_id])


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    parser.add_argument("--gene_list", nargs="+", type=str, help="gene list(s) to analyse", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gene_ids = load_transcripts(args.genedb)

    for gene_list in args.gene_lists:
        fix_ids(gene_ids, gene_list, gene_list + ".fixed.tsv")



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
