############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import argparse
from assign_isoforms_to_barcodes import *
from traceback import print_exc
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def hamming_rel_dist(l1, l2):
    if len(l1) != len(l2):
        return 1
    return float(hamming(l1, l2)) / float(len(l1))


class IsoformPhylogeneticTreeConstructor:

    def __init__(self, db, gene_list):
        self.db = db
        self.chr_id = ""

        gene_db_list = []
        for gene_name in gene_list:
            gene_db = self.db[gene_name]
            if self.chr_id == "":
                self.chr_id = gene_db.seqid
            elif self.chr_id != gene_db.seqid:
                print("Gene " + gene_name + " is from different chromosome")
                continue

            gene_db_list.append(gene_db)

        self.gene_info = GeneInfo(gene_db_list, db)

    def construct_distance_matrix(self):
        exon_profiles = self.gene_info.all_rna_profiles.exon_profiles
        self.isoform_list = sorted(exon_profiles.keys())
        dm = []
        for i1 in range(len(self.isoform_list)):
            matrix_row = []
            for i2 in range(i1):
                isoform1 = self.isoform_list[i1]
                isoform2 = self.isoform_list[i2]
                matrix_row.append(hamming_rel_dist(exon_profiles[isoform1][1:-1], exon_profiles[isoform2][1:-1]))
            matrix_row.append(0.0)
            dm.append(matrix_row)

        self.distance_matrix = DistanceMatrix(self.isoform_list, dm)
        print(self.distance_matrix)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(self.distance_matrix)
        print(tree)
        tree2 = constructor.upgma(self.distance_matrix)
        print(tree2)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gene_ids', metavar='BAM_FILE', nargs='+', type=str,  help='gene ids for tree construction')
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if not os.path.isfile(args.genedb):
        raise Exception("Gene database " + args.genedb + " does not exist")
    db = gffutils.FeatureDB(args.genedb, keep_order=True)

    isoform_tree_constructor = IsoformPhylogeneticTreeConstructor(db, args.gene_ids)
    isoform_tree_constructor.construct_distance_matrix()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)





