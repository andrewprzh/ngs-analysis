#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import gffutils
import argparse
from traceback import print_exc
from collections import defaultdict


def transcript2gene(db):
    transcript2gene_dict = {}
    gene_db = gffutils.FeatureDB(db)
    for g in gene_db.features_of_type('gene', order_by=('seqid', 'start')):
        for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):
            transcript2gene_dict[t.id] = g.id
    return transcript2gene_dict


def convert_counts(transcript_counts_file_name, transcript2gene, gene_counts_file_name):
    with open(gene_counts_file_name, "w") as outf:
        gene_counts = defaultdict(float)
        for l in open(transcript_counts_file_name):
            if l.startswith("#"):
                outf.write(l)
                continue
            v = l.strip().split("\t")
            gene_counts[transcript2gene[v[0]]] += float(v[1])

        scale_factor = sum(gene_counts.values()) / 1000000.0
        print("Saving counts to %s" % gene_counts_file_name)
        total_tpm = 0
        for gid in sorted(gene_counts.keys()):
            tpm = gene_counts[gid] / scale_factor
            total_tpm += tpm
            outf.write("%s\t%.2f\t%.6f\n" % (gid, gene_counts[gid], tpm))
        print("Done. Total counts: %d; total TPM: %.2f" % (sum(gene_counts.values()), total_tpm))




def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with gene counts", required=True)
    parser.add_argument("--input", "-i", type=str, help="input file with transcript counts", required=True)
    parser.add_argument('--genedb', "-g", type=str, help="gffutils gene db", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    transcript2gene_dict = transcript2gene(args.genedb)
    convert_counts(args.input, transcript2gene_dict, args.output)
    

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
