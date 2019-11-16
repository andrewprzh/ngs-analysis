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
from common import *


class StatCounter:

    def __init__(self, args):
        self.args = args
        self.seq_set = set()
        for record in SeqIO.parse(args.fasta, "fasta"):
            self.seq_set.add(record.id)

        self.mapped_seqs = {}
        self.secondary_mapped_seqs = {}
        bamfile_in = pysam.AlignmentFile(args.bam, "rb")
        for alignment in bamfile_in:
            if alignment.reference_id == -1:
                continue
            blocks = sorted(alignment.get_blocks())
            if len(blocks) == 0:
                continue
            if alignment.is_secondary:
                self.secondary_mapped_seqs.add[alignment.query_name] = (blocks[0][0], blocks[-1][1])
            else:
                self.mapped_seqs.add[alignment.query_name] = (blocks[0][0], blocks[-1][1])

        self.assigned_isoforms = {}
        for l in open(args.tsv):
            tokens = l.strip().split()
            if len(tokens) != 2:
                print("Error processing " + l)
                continue
            seq_id = tokens[0]
            if seq_id in seq_set:
                self.assigned_isoforms[seq_id] = tokens[1]

        self.db = gffutils.FeatureDB(args.gene_db, keep_order=True)
        self.isoform_to_gene_map = {}
        self.gene_to_isoforms_map = {}
        for g in self.db.features_of_type('gene'):
            gene_name = g.id
            self.gene_to_isoforms_map[gene_name] = set()
            gene_db = self.db[gene_name]
            for t in self.db.children(gene_db, featuretype='transcript'):
                self.isoform_to_gene_map[t.id] = gene_name
                self.gene_to_isoforms_map[gene_name].add(t.id)

    def seqid_to_isoform(self, seq):
        return seq.split('_')[0]

    def count_mapping_stats(sefl, mapped_seqs):
        mapped_count = 0
        correctly_mapped = 0
        for seq in self.seq_set:
            if seq in mapped_seqs:
                mapped_count += 1
                gene_name = self.isoform_to_gene_map[sefl.seqid_to_isoform(seq)]
                gene_db = self.db[gene_name]
                if overlaps(self.mapped_seqs[seq], (gene_db.start, gene_db.end)):
                    correctly_mapped += 1
                #TODO add stats on mismapped genes
        return mapped_count, correctly_mapped

    def count_assignment_stats(sefl):
        #TODO: check whether correct isoforms are assigned
        correct_assignments = 0
        incorrect_same_gene = 0
        incorrect_other_gene = 0

        for seq in assigned_isoforms.keys():
            isoform_id = self.seqid_to_isoform(seq)
            if isoform_id == assigned_isoforms[seq]:
                correct_assignments += 1
            else:
                gene_name = self.isoform_to_gene_map[isoform_id]
                if isoform_id in self.gene_to_isoforms_map[gene_name]:
                    incorrect_same_gene += 1
                else:
                    incorrect_other_gene += 1
        return correct_assignments, incorrect_same_gene, incorrect_other_gene


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="initial file with sequences", type=str)
    parser.add_argument("--bam", "-b", help="mapped sequences", type=str)
    parser.add_argument("--tsv", "-t", help="assigned isoforms", type=str)
    parser.add_argument("--gene_db", "-g", help="gene database", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    stat_counter = StatCounter(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

