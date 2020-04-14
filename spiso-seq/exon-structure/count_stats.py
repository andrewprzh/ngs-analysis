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
        self.seqid_to_isoform = dict()
        for record in SeqIO.parse(args.fasta, "fasta"):
            if '_' in record.description:
                fs = record.description.split('_')
                seq_id, isoform_id = record.description, fs[1]
            else:  # SQANTI2 requires sequence id as >PB.X.Y <isoform_id>
                fs = record.description.split()
                seq_id, isoform_id = fs[0], fs[1]
            self.seqid_to_isoform[seq_id] = isoform_id
            self.seq_set.add(seq_id)

        self.mapped_seqs = {}
        self.secondary_mapped_seqs = {}
        alignment_file_in = pysam.AlignmentFile(args.bam or args.sam, "rb")
        for alignment in alignment_file_in:
            if alignment.reference_id == -1:
                continue
            blocks = sorted(alignment.get_blocks())
            if len(blocks) == 0:
                continue
            if alignment.is_secondary or alignment.is_supplementary:
                self.secondary_mapped_seqs[alignment.query_name] = (blocks[0][0], blocks[-1][1])
            else:
                self.mapped_seqs[alignment.query_name] = (blocks[0][0], blocks[-1][1])

        self.assigned_isoforms = {}
        for l in open(args.tsv):
            tokens = l.strip().split()
            seq_id = tokens[0]
            if seq_id in self.seq_set:
                if len(tokens) == 2:
                    self.assigned_isoforms[seq_id] = tokens[1]
                elif len(tokens) > 10:  ## SQANTI2
                    self.assigned_isoforms[seq_id] = tokens[7]
                else:
                    print("Error processing " + l)
                    continue

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

    # def seqid_to_isoform(self, seq):
    #    return seq.split('_')[0]

    def count_mapping_stats(self):
        mapped_count = 0
        correctly_mapped = 0
        mismapped_seqs = set()
        for seq_id in self.seq_set:
            if seq_id in self.mapped_seqs:
                mapped_count += 1
                gene_name = self.isoform_to_gene_map[self.seqid_to_isoform[seq_id]]
                gene_db = self.db[gene_name]
                if overlaps(self.mapped_seqs[seq_id], (gene_db.start, gene_db.end)):
                    correctly_mapped += 1
                else:
                    mismapped_seqs.add(seq_id)
                #TODO add stats on mismapped genes
        return mapped_count, correctly_mapped, mismapped_seqs

    def count_assignment_stats(self, mismapped_seqs=None):
        correct_assignments = 0
        incorrect_same_gene = 0
        incorrect_other_gene = 0
        incorrect_novel = 0  # SQANTI2 output

        '''with open("unassigned.txt", "w") as f:
            for seq_id in self.seq_set:
                if seq_id not in self.assigned_isoforms:
                    f.write(seq_id + "\n")'''
        for seq_id in self.assigned_isoforms.keys():
            real_isoform_id = self.seqid_to_isoform[seq_id]
            assigned_isoform_id = self.assigned_isoforms[seq_id]
            if mismapped_seqs and seq_id in mismapped_seqs:
                continue
            if assigned_isoform_id == real_isoform_id:
                correct_assignments += 1
            elif assigned_isoform_id == 'novel':
                incorrect_novel += 1
            else:
                gene_name = self.isoform_to_gene_map[real_isoform_id]
                if assigned_isoform_id in self.gene_to_isoforms_map[gene_name]:
                    incorrect_same_gene += 1
                else:
                    incorrect_other_gene += 1
        return correct_assignments, incorrect_same_gene, incorrect_other_gene, incorrect_novel


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="initial file with sequences", type=str)
    parser.add_argument("--bam", "-b", help="mapped sequences", type=str)
    parser.add_argument("--sam", "-s", help="mapped sequences", type=str) ## SQANTI2 output
    parser.add_argument("--tsv", "-t", help="assigned isoforms", type=str)
    parser.add_argument("--gene_db", "-g", help="gene database", type=str)

    args = parser.parse_args()
    return args


def calc_precision(tp, fp):
    return tp * 100.0 / (tp + fp)


def calc_recall(tp, fn):
    return tp * 100.0 / (tp + fn)


def main():
    args = parse_args()
    stat_counter = StatCounter(args)
    mapped_count, correctly_mapped, mismapped_seqs = stat_counter.count_mapping_stats()
    total_reads = len(stat_counter.seq_set)
    mismapped_reads = mapped_count - correctly_mapped
    unmapped_reads = total_reads - mapped_count
    print("Correct mapping %d, incorrect mapping %d, unmapped %d" %
          (correctly_mapped, mismapped_reads, unmapped_reads))
    print("Mapping precision: %f, recall %f" %
          (calc_precision(correctly_mapped, mismapped_reads),
           calc_recall(correctly_mapped, unmapped_reads)))

    # use only correctly mapped reads
    correct_assignments, incorrect_same_gene, incorrect_other_gene, incorrect_novel = \
        stat_counter.count_assignment_stats(mismapped_seqs)
    assigned_reads = correct_assignments + incorrect_same_gene + incorrect_other_gene
    unassigned_mapped_reads = correctly_mapped - assigned_reads
    print("Correct isoform %d, incorrect assignment %d, unassigned %d" %
          (correct_assignments, assigned_reads - correct_assignments, unassigned_mapped_reads))
    print("Isoform assignment precision: %f, recall %f" %
          (calc_precision(correct_assignments, incorrect_same_gene+incorrect_other_gene),
           calc_recall(correct_assignments, unassigned_mapped_reads)))

    # use all reads
    correct_assignments, incorrect_same_gene, incorrect_other_gene, incorrect_novel = \
        stat_counter.count_assignment_stats()
    assigned_reads = correct_assignments + incorrect_same_gene + incorrect_other_gene
    wrong_assignments = assigned_reads - correct_assignments
    unassigned_total_reads = total_reads - assigned_reads
    print("ALL READS Correct isoform %d, incorrect assignment %d, unassigned %d" %
          (correct_assignments, wrong_assignments, unassigned_total_reads))
    print("ALL READS Isoform assignment precision: %f, recall: %f" %
          (calc_precision(correct_assignments, wrong_assignments),
           calc_recall(correct_assignments, unassigned_total_reads)))

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

