############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import pysam
from Bio import SeqIO
import gffutils


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="initial file with sequences", type=str)
    parser.add_argument("--bam", "-b", help="mapped sequences", type=str)
    parser.add_argument("--tsv", "-t", help="assigned isoforms", type=str)
    parser.add_argument("--gene_db", "-g", help="gene database", type=str)

    args = parser.parse_args()
    return args


def count_subset_stats(seq_set, mapped_seqs):
    mapped_count = 0
    for seq in seq_set:
        if seq in mapped_seqs:
            mapped_count += 1
    return mapped_count


def count_mapping_stats(mapped_isoforms, gene_db):
    pass
    #TODO: check whether sequence mapped to the same gene

def count_assignment_stats(assigned_isoforms):
    #TODO: check whether correct isoforms are assigned
    correct_assignments = 0

    for seq in assigned_isoforms.keys():
        if seq.split('_') == assigned_isoforms[seq]


#TODO: print ambiguous genes

def main():
    args = parse_args()

    seq_set = set()
    for record in SeqIO.parse(args.fasta, "fasta"):
        seq_set.add(record.id)

    #TODO: add mapping coords
    mapped_seqs = {}
    bamfile_in = pysam.AlignmentFile(args.bam, "rb")
    for alignment in bamfile_in:
        if alignment.reference_id == -1 or alignment.is_secondary:
            continue
        mapped_seqs.add(alignment.query_name)

    assigned_isoforms = {}
    for l in open(args.tsv):
        tokens = l.strip().split()
        if len(tokens) != 2:
            print("Error processing " + l)
            continue
        seq_id = tokens[0]
        if seq_id in seq_set:
            assigned_isoforms[seq_id] = tokens[1]

    gene_db = gffutils.FeatureDB(args.gene_db, keep_order=True)
    for g in self.db.features_of_type('gene', order_by=('seqid', 'start')):
        if current_chromosome != g.seqid:
            current_chromosome = g.seqid
            print("Processing chromosome " + current_chromosome)
        gene_name = g.id
        gene_db = self.db[gene_name]



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

