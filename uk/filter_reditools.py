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


class PositionFilter:
    def __init__(self):
        # chr + strand -> position
        self.position_dict = defaultdict(set)

    def update(self, reditools_output):
        for l in open(reditools_output):
            v = l.strip().split("\t")
            rna_sub = v[7]
            try:
                rna_freq = float(v[8])
            except ValueError:
                rna_freq = 0.0
            dna_sub = v[12]
            try:
                dna_freq = float(v[13])
            except ValueError:
                dna_freq = 0.0
            position = int(v[1])
            chr_id = v[0]

            if rna_sub != dna_sub or rna_freq > dna_freq:
                if rna_sub == "AG":
                    self.position_dict[chr_id + "+"].add(position)
                elif rna_sub == "TC":
                    self.position_dict[chr_id + "-"].add(position)

    def filter_table(self, reditools_output, filtered_output_name):
        with open(filtered_output_name, "w") as outf:
            for l in open(reditools_output):
                v = l.strip().split("\t")
                position = int(v[1])
                chr_id = v[0]
                if position in self.position_dict[chr_id + "+"] or position in self.position_dict[chr_id + "-"]:
                    outf.write(l)

    def filter_genes(self, genedb_filename):
        gffutils_db = gffutils.FeatureDB(genedb_filename)
        gene_dict = defaultdict(list)
        for g in gffutils_db.features_of_type('gene', order_by=('seqid', 'start')):
            gene_dict[g.seqid + g.strand].append((g.start, g.end, g.id))

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

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    position_filter = PositionFilter()
    for inf in args.tables:
        print("Processing %s" % inf)
        position_filter.update(inf)

    if args.genedb:
        print("Filtering using gene annotation %s" % args.genedb)
        position_dict = position_filter.filter_genes(args.genedb)

        with open("position_info.tsv", "w") as outf:
            for pos in position_dict.keys():
                outf.write("%s\t%s\t%d\t%s\n" % (pos[0], pos[1], pos[2], position_dict[pos]))

    for inf in args.tables:
        filtered_table = inf + "filtered.tsv"
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
