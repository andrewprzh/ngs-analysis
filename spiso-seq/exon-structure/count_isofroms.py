############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import operator
import gffutils
import argparse
from traceback import print_exc

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", help="gene database in gffutils db format", type=str)
    parser.add_argument("--geneid", help="gene id", type=str)
    parser.add_argument("--infile", help="isoform assignments", type=str)

    args = parser.parse_args()
    return args


def parse_assignments(inf, gene_name, isoform_counts):
    gene_found = False
    for l in open(inf):
        if l.startswith('E'):
            if gene_found:
                gene_found = False
                continue
            elif l.find(gene_name) != -1:
                gene_found = True
        elif not gene_found:
            continue

        tokens = l.strip().split()
        if len(tokens) != 2:
            continue

        isoform_id = tokens[1]
        if isoform_id not in isoform_counts:
            continue
        read_id = tokens[0]
        cell_group = read_id.split(':')[0]
        if cell_group not in isoform_counts[isoform_id]:
            isoform_counts[isoform_id][cell_group] = 0
        isoform_counts[isoform_id][cell_group] += 1


def print_counts(isoform_counts, outf):
    f = open(outf, 'w')
    for i in isoform_counts:
        for cg in isoform_counts[i]:
            f.write(i + '\t' + cg + '\t' + str(isoform_counts[i][cg]) + '\n')
    f.close()


def main():
    args = parse_args()

    if not os.path.isfile(args.genedb):
        raise Exception("Gene database " + args.genedb + " does not exist")
    db = gffutils.FeatureDB(args.genedb, keep_order=True)

    gene_name = args.geneid
    gene_db = db[gene_name]
    isoform_counts = {}
    for t in db.children(gene_db, featuretype='transcript'):
        if t.id not in isoform_counts:
            isoform_counts[t.id] = {}

    parse_assignments(args.infile, gene_name, isoform_counts)
    print_counts(isoform_counts, gene_name + '.iso_counts.tsv')


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
