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


MIN_NODE_SIZE = 40
MAX_NODE_SIZE = 600
EDGE_LEN_SCALE = 10
DEFAULT_NODE_SIZE = 200

CELL_GROUPS = ['P7Hipp_ExcitNeuron', 'P7Hipp_InhibNeuron', 'P7Hipp_GranuleNB', 'P7PFC_ExcitNeuron', 'P7PFC_InhibNeuron']

def hamming_rel_dist(l1, l2):
    if len(l1) != len(l2):
        return 1
    return float(hamming(l1, l2)) / float(len(l1))


class IsoformPhylogeneticTreeConstructor:

    # cell_groups_isoform_counts: cell_group -> {isoform -> counts}
    def __init__(self, db, gene_list, cell_groups_isoform_counts = {}):
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

        self.exon_profiles = self.gene_info.all_rna_profiles.exon_profiles
        self.cell_groups_isoform_counts = cell_groups_isoform_counts
        if len(self.cell_groups_isoform_counts) == 0:
            self.all_isoform_list = sorted(exon_profiles.keys())
        else:
            self.all_isoform_list = set()
            self.max_count = 0
            self.min_count = 100000000
            isoform_totals = {}

            self.max_isoform = ""
            for cg in self.cell_groups_isoform_counts.keys():
                self.max_count = max(self.max_count, max(self.cell_groups_isoform_counts[cg].values()))
                self.min_count = min(self.min_count, min(self.cell_groups_isoform_counts[cg].values()))
                for isoform in self.cell_groups_isoform_counts[cg].keys():
                    self.all_isoform_list.add(isoform)
                    if isoform not in isoform_totals:
                        isoform_totals[isoform] = 0
                    isoform_totals[isoform] += self.cell_groups_isoform_counts[cg][isoform]
            self.all_isoform_list = sorted(list(self.all_isoform_list))
            self.max_isoform = max(isoform_totals.iteritems(), key=operator.itemgetter(1))[0]

            count_div = self.max_count - self.min_count
            if count_div == 0:
                self.size_scale = MAX_NODE_SIZE / (2 * self.max_count)
            else:
                self.size_scale = (MAX_NODE_SIZE - MIN_NODE_SIZE) / count_div

    def make_tree(self, isoform_list):
        dm = []
        for i1 in range(len(isoform_list)):
            matrix_row = []
            for i2 in range(i1):
                isoform1 = isoform_list[i1]
                isoform2 = isoform_list[i2]
                matrix_row.append(hamming_rel_dist(self.exon_profiles[isoform1][1:-1], self.exon_profiles[isoform2][1:-1]))
            matrix_row.append(0.0)
            dm.append(matrix_row)

        distance_matrix = DistanceMatrix(isoform_list, dm)
        constructor = DistanceTreeConstructor()
        return constructor.upgma(distance_matrix)

    def draw_tree(self, tree, isoform_counts = {}, xshift = 0, yshift = 0, cell_group_name = ''):
        tree.root_with_outgroup(self.max_isoform)
        tree.clade.name = "InnerX"
        net = Phylo.to_networkx(tree)

        nice_tree = networkx.DiGraph()
        for node in net.nodes:
            nice_tree.add_node(node.name)
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
                if len(isoform_counts) > 0:
                    if node not in isoform_counts:
                        size = MIN_NODE_SIZE
                    else:
                        size = MIN_NODE_SIZE + (isoform_counts[node] - self.min_count) * self.size_scale
                    node_sizes.append(size)
                    node_colors.append('red' if node == self.max_isoform else 'green')
                else:
                    node_sizes.append(DEFAULT_NODE_SIZE)
                    node_colors.append('green')
                    node_labels[node] = node


        pos = networkx.drawing.nx_pydot.pydot_layout(nice_tree)
        new_pos = {}
        for k in pos:
            new_pos[k] = (pos[k][0] + xshift, pos[k][1] + yshift)
        # TODO: manually draw labels
        networkx.draw(nice_tree, new_pos, with_labels=True, labels=node_labels, node_color=node_colors,
                      node_size=node_sizes, font_size = '4')
#        plt.show()

    def draw_trees(self, outf_name):
        hipp_count = 0
        pfc_count = 0
        for cg in CELL_GROUPS:
            yshift = 0
            xshift = 0
            if cg.startswith("P7Hipp"):
                xshift = 600 * hipp_count
                hipp_count += 1
            else:
                yshift = 400
                xshift = 600 * pfc_count
                pfc_count += 1

            t = self.make_tree(self.all_isoform_list)
            self.draw_tree(t, self.cell_groups_isoform_counts[cg], xshift, yshift)
        plt.savefig(outf_name, format='svg')


def read_counts(inf):
    isoform_counts = {}
    for l in open(inf):
        tokens = l.strip().split()
        if len(tokens) != 3:
            continue
        isoform = tokens[0]
        cg = tokens[1]
        for c in CELL_GROUPS:
            if cg.startswith(c):
                cg = c
        count = int(tokens[2])
        if cg not in isoform_counts:
            isoform_counts[cg] = {}
        if isoform not in isoform_counts[cg]:
            isoform_counts[cg][isoform] = 0
        isoform_counts[cg][isoform] += count
    return isoform_counts


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gene_id', type=str,  help='gene ids for tree construction')
    parser.add_argument('--counts', type=str, help='counts')
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if not os.path.isfile(args.genedb):
        raise Exception("Gene database " + args.genedb + " does not exist")
    db = gffutils.FeatureDB(args.genedb, keep_order=True)

    isoform_counts = read_counts(args.counts)
    isoform_tree_constructor = IsoformPhylogeneticTreeConstructor(db, [args.gene_id], isoform_counts)
    isoform_tree_constructor.draw_trees(args.gene_id + ".svg")

#    isoform_counts = {'ENSMUST00000217664.1' : 61, 'ENSMUST00000217810.1' : 215, 'ENSMUST00000218469.1' : 29, 'ENSMUST00000218480.1' : 81}
#    isoform_tree_constructor.make_tree(isoform_counts)
#    isoform_tree_constructor.draw_tree()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)





