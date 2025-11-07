############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Region  Position        Reference       Strand  Coverage-q30    MeanQ   BaseCount[A,C,G,T]      AllSubs Frequency       gCoverage-q30   gMeanQ  gBaseCount[A,C,G,T]     gAllSubs        gFrequency
# GL000219.1      54843   T       2       49      46.53   [0, 0, 0, 49]   -       0.00    38      39.24   [0, 0, 0, 38]   -       0.00

import os
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import gffutils
import glob


class PositionFilter:
    def __init__(self, args):
        # chr + strand -> position
        self.position_dict = defaultdict(set)
        self.position_rna_cov = defaultdict(list)
        self.position_rna_sub_cov = defaultdict(list)
        self.args = args

    def update(self, reditools_output):
#        self.position_dict = defaultdict(set)
        for l in open(reditools_output):
            if l.startswith("Region"): continue
            v = l.strip().split("\t")
            rna_sub = v[7]
            try:
                rna_freq = float(v[8])
            except ValueError:
                rna_freq = 0.0
            try:
                rna_cov = int(v[4])
            except ValueError:
                rna_cov = 0

            dna_sub = v[12]
            try:
                dna_freq = float(v[13])
            except ValueError:
                dna_freq = 0.0
            try:
                dna_cov = int(v[9])
            except ValueError:
                dna_cov = 0
            position = int(v[1])
            chr_id = v[0]
            strand = v[3]

            if dna_cov < self.args.min_dna_cov:
                continue

            self.position_rna_cov[(chr_id, position)].append(rna_cov)
            self.position_rna_sub_cov[(chr_id, position)].append(rna_cov * rna_freq)

            if rna_sub != dna_sub or rna_freq > dna_freq:
                if rna_sub == "AG":
                     self.position_dict[chr_id + strand].add(position)
                elif rna_sub == "TC":
                    self.position_dict[chr_id + "-"].add(position)

    def filter_coverage(self):
        new_position_dict = defaultdict(set)
        for chr_strand in self.position_dict:
            for pos in self.position_dict[chr_strand]:
                chr_id = chr_strand[:-1]
                if sum(self.position_rna_sub_cov[(chr_id, pos)]) < len(self.position_rna_sub_cov[(chr_id, pos)]) / 2:
                    continue
                if self.args.min_cov_strategy == 1:
                    if all(x >= self.args.min_rna_cov for x in self.position_rna_cov[(chr_id, pos)]):
                        new_position_dict[chr_strand].add(pos)
                elif self.args.min_cov_strategy == 2:
                    if any(x >= self.args.min_rna_cov for x in self.position_rna_cov[(chr_id, pos)]):
                        new_position_dict[chr_strand].add(pos)
                elif self.args.min_cov_strategy == 3:
                    cov_vals = self.position_rna_cov[(chr_id, pos)]
                    if sum(cov_vals) >= self.args.min_rna_cov * len(cov_vals) :
                        new_position_dict[chr_strand].add(pos)
                else:
                    print("Unsupported srtatgy %d" % self.args.min_cov_strategy)
                    exit(-1)
        self.position_dict = new_position_dict

    def filter_table(self, reditools_output, filtered_output_name):
        with open(filtered_output_name, "w") as outf:
            outf.write("Region\tPosition\tReference\tStrand\tCoverage-q30\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q30\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\n")
            for l in open(reditools_output):
                if l.startswith("Region"):
                    continue
                v = l.strip().split("\t")
                # if v[9] == v[10] == '-':
                #    continue
                position = int(v[1])
                chr_id = v[0]
                if position in self.position_dict[chr_id + "+"] or position in self.position_dict[chr_id + "-"]:
                    outf.write(l)

    def filter_genes(self, genedb_filename):
        gffutils_db = gffutils.FeatureDB(genedb_filename)
        gene_dict = defaultdict(list)
        for g in gffutils_db.features_of_type('gene', order_by=('seqid', 'start')):
            gene_dict[g.seqid + g.strand].append((g.start, g.end, g.id))

        print("Total genes loaded %d" % sum([len(s) for s in gene_dict.values()]))

        new_dict = defaultdict(set)
        pos_gene_dict = {}
        for chr_info in self.position_dict.keys():
            if chr_info not in gene_dict:
                continue
            gene_list = gene_dict[chr_info]
            pos_list = list(sorted(self.position_dict[chr_info]))
            pos_index = 0
            gene_index = 0

            while pos_index < len(pos_list) and gene_index < len(gene_list):
                position = pos_list[pos_index]
                gene = gene_list[gene_index]
                if position < gene[0]:
                    pos_index += 1
                elif position > gene[1]:
                    gene_index += 1
                else:
                    pos_index += 1
                    new_dict[chr_info].add(position)
                    pos_gene_dict[(chr_info[:-1], chr_info[-1], position)] = gene[2]

        self.position_dict = new_dict
        return pos_gene_dict


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tables", "-t", nargs="+", help="REDItools output tables", required=True, type=str)
    parser.add_argument("--genedb", "-g", help="gffutils genedb (can be obtained using gtf2db.py);"
                                               "no gene filtring will be performed if not given", type=str)
    parser.add_argument("--output", "-o", help="output folder", type=str, default="./")
    parser.add_argument("--min_rna_cov", help="minimal # RNA reads", type=int, default=10)
    parser.add_argument("--min_cov_strategy", help="minimal # RNA reads is required in 1: all samples, 2: average, 3: at least one (default: all)", type=int, default=1)
    parser.add_argument("--min_dna_cov", help="minimal # DNA reads", type=int, default=1)


    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    position_filter = PositionFilter(args)
    for f in args.tables:
        for inf in glob.glob(f):
            print("Processing %s" % inf)
            position_filter.update(inf)

    print("Total positions collected %d" % sum([len(s) for s in position_filter.position_dict.values()]))

    position_filter.filter_coverage()
    print("Total positions after filtering by coverage %d" % sum([len(s) for s in position_filter.position_dict.values()]))

    if args.genedb:
        print("Filtering using gene annotation %s" % args.genedb)
        position_dict = position_filter.filter_genes(args.genedb)

        print("Total positions after gene filtering %d" % sum([len(s) for s in position_filter.position_dict.values()]))

        with open("position_info.tsv", "w") as outf:
            for pos in position_dict.keys():
                outf.write("%s\t%s\t%d\t%s\n" % (pos[0], pos[1], pos[2], position_dict[pos]))

    for inf in args.tables:
        filtered_table = os.path.join(args.output, os.path.basename(inf) + ".filtered.tsv")
        print("Filtering %s to %s" % (inf, filtered_table))
        position_filter.filter_table(inf, filtered_table)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
