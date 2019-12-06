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
from assign_isoforms_to_barcodes import *
from traceback import print_exc
import Bio.Phylo as Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import networkx
import matplotlib.pyplot as plt


MIN_NODE_SIZE = 20
MAX_NODE_SIZE = 1000
EDGE_LEN_SCALE = 20
DEFAULT_NODE_SIZE = 200

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

    def make_tree(self, isoform_counts = {}):
        exon_profiles = self.gene_info.all_rna_profiles.exon_profiles
        self.isoform_counts = isoform_counts
        if len(self.isoform_counts) == 0:
            self.isoform_list = sorted(exon_profiles.keys())
        else:
            self.isoform_list = sorted(self.isoform_counts.keys())

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
        constructor = DistanceTreeConstructor()
        self.tree = constructor.upgma(self.distance_matrix)

    def draw_tree(self):
        max_isoform = ""
        if len(self.isoform_counts) > 0:
            max_isoform = max(self.isoform_counts.iteritems(), key=operator.itemgetter(1))[0]
        if max_isoform != "":
            self.tree.root_with_outgroup(max_isoform)
        self.tree.clade.name = "InnerX"
        net = Phylo.to_networkx(self.tree)

        nice_tree = networkx.DiGraph()
        for node in net.nodes:
            nice_tree.add_node(node.name, size='100')
        for e in net.edges:
            edge_len = EDGE_LEN_SCALE * e[1].branch_length
            if edge_len == 0:
                edge_len = 0.1
            nice_tree.add_edge(e[0].name, e[1].name, len=edge_len)

        node_sizes = []
        node_colors = []
        node_labels = {}
        for node in nice_tree.nodes:
            if node.startswith("Inner"):
                node_sizes.append(1)
                node_colors.append('black')
                node_labels[node] = ' '
            else:
                node_labels[node] = node
                if len(self.isoform_counts) > 0:
                    count_div = max(self.isoform_counts.values()) - min(self.isoform_counts.values())
                    if count_div == 0:
                        size_scale = MAX_NODE_SIZE / (2 * max(self.isoform_counts.values()))
                    else:
                        size_scale = (MAX_NODE_SIZE - MIN_NODE_SIZE) / count_div
                    size = MIN_NODE_SIZE + (self.isoform_counts[node] - min(self.isoform_counts.values())) * size_scale
                    node_sizes.append(size)
                    node_colors.append('red' if node == max_isoform else 'green')
                else:
                    node_sizes.append(DEFAULT_NODE_SIZE)
                    node_colors.append('green')
                    node_labels[node] = node

        pos = networkx.drawing.nx_pydot.pydot_layout(nice_tree)
        networkx.drawing.nx_pydot.write_dot(nice_tree, "1.dot")
        networkx.draw(nice_tree, pos, with_labels = True, labels=node_labels, node_color=node_colors, node_size=node_sizes)
        plt.show()


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
    isoform_counts = {'ENSMUST00000217664.1' : 61, 'ENSMUST00000217810.1' : 215, 'ENSMUST00000218469.1' : 29, 'ENSMUST00000218480.1' : 81}
    isoform_tree_constructor.make_tree(isoform_counts)
    isoform_tree_constructor.draw_tree()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)





