#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gffutils
from traceback import print_exc
import glob


def equals_any_feature(exon, feature_list):
    for f in feature_list:
        if exon == f:
            return True
    return False


def overlaps_any_feature(exon, feature_list):
    ovlp = 0.0
    best_feature = None
    for f in feature_list:
        feature_overlap = (min(exon[1], f[1]) - max(exon[0], f[0]) + 1) / (exon[1] - exon[0] + 1)
        if feature_overlap > ovlp:
            ovlp = feature_overlap
            best_feature = f
    return best_feature


def process_gene(gene_db, gene):
    exons = set()
    if isinstance(gene, str):
        g = gene_db[gene]
    else:
        g = gene
    for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):
        for c in gene_db.children(t, featuretype='exon', order_by='start'):
            exons.add((c.start, c.end))

    return exons


def load_exon_info(genedb_path):
    gene_db = gffutils.FeatureDB(genedb_path)
    gene_dicts = []
    for g in gene_db.features_of_type('gene', order_by=('seqid', 'start')):
        gene_dicts[g.id] = process_gene(gene_db, g)
    return gene_dicts


def load_exon_pairs(infile):
    exon_pairs = []
    for l in open(infile):
        if l.startswith("#"): continue
        v = l.strip().split()
        exon_pairs.append((v[0], v[1]))
    return exon_pairs


def process_exons(exon_gene_pairs, gene_db, gene_dicts):
    exon_dict = {}
    for exon_pair in exon_gene_pairs:
        gene_id = exon_pair[1]
        if gene_id not in gene_dicts:
            gene_dicts[gene_id] = process_gene(gene_db, gene_id)

        exon_id = exon_pair[0]
        exon_data = exon_id.split("_")
        exon = (int(exon_data[1]), int(exon_data[2]))
        exons = gene_dicts[gene_id]

        if equals_any_feature(exon, exons):
            exon_dict[exon_id] = (gene_id, exon_id)
            continue

        if exon_id in exon_dict and exon_dict[exon_id][1] == exon_dict:
            continue

        overlapping_exon = overlaps_any_feature(exon, exons)
        if overlapping_exon is None:
            print("No matches found yet %d-%d" % (exon[0], exon[1]))
            exon_dict[exon_id] = (gene_id, "None")

        if overlapping_exon[0] != exon[0] and overlapping_exon[1] != exon[1]:
            print("Unequal splice junctions for %d-%d and %d-%d" % (exon[0], exon[1], overlapping_exon[0], overlapping_exon[1]))
        exon_dict[exon_id] = (gene_id, "_".join([exon_data[0], str(overlapping_exon[0]),
                                                 str(overlapping_exon[1]), exon_data[3]]))
    return exon_dict


def dump_dict(exon_dict, out_fname):
    with open(out_fname, "w") as outf:
        for e in sorted(exon_dict.keys()):
            exon_pair = exon_dict[e]
            outf.write("%s\t%s\t%s\n" % (e, exon_pair[0], exon_pair[1]))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output dir", required=True)
    parser.add_argument("--genedb", type=str, help="gffutils gene DB", required=True)
    parser.add_argument("--input", "-i", type=str, nargs='+',
                        help="one or more files/dirs with input exon lists", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gene_db = gffutils.FeatureDB(args.genedb)
    gene_dicts = {}


    for inf in args.input:
        if os.path.isdir(inf):
            for f in glob.glob(os.path.join(inf, "*")):
                if os.path.isfile(f):
                    name = os.path.basename(f)
                    exon_dict = process_exons(load_exon_pairs(f), gene_db, gene_dicts)
                    dump_dict(exon_dict, os.path.join(args.output, name + ".normalized.tsv"))
        elif os.path.isfile(inf):
            name = os.path.basename(inf)
            exon_dict = process_exons(load_exon_pairs(inf), gene_db, gene_dicts)
            dump_dict(exon_dict, os.path.join(args.output, name + ".normalized.tsv"))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
