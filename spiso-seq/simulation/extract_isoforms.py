############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from Bio import SeqIO
import numpy
from traceback import print_exc
import gffutils


def extract_isoforms(args, isoforms):
    isoform_records = []

    for record in SeqIO.parse(args.fasta, 'fasta'):
        isoform_id = record.id.strip()
        if isoform_id not in isoforms:
            continue
        isoform_records.append(record)

    SeqIO.write(isoform_records, args.output, 'fasta')


def get_isoform_set(args):
    isoforms = set()
    if args.genedb is None or args.gene_id is None:
        return isoforms

    db = gffutils.FeatureDB(args.genedb)
    gene = db[args.gene_id]

    for t in db.children(gene, featuretype='transcript'):
        isoforms.add(t.id)
    return isoforms


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="initial file with sequences", type=str)
    parser.add_argument("--output", "-o", help="output folder", type=str)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--gene_id", help="gene id to take isoforms from", type=str)


    args = parser.parse_args()

    if args.fasta is None or args.genedb is None or args.gene_id is None or args.output is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()
    isoforms = get_isoform_set(args)
    extract_isoforms(args, isoforms)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)