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


def contains_any_feature(exon, feature_list, delta=0):
    for f in feature_list:
        if exon[0] - delta <= f[0] and f[1] + delta <= exon[1]:
            return True
    return False


def equals_any_feature(exon, feature_list):
    for f in feature_list:
        if exon == f:
            return True
    return False


def overlaps_any_feature(exon, feature_list):
    ovlp = 0.0
    for f in feature_list:
        ovlp = max(ovlp, (min(exon[1], f[1]) - max(exon[0], f[0]) + 1) / (exon[1] - exon[0] + 1))
    return ovlp


class ExonInfo:
    def __init__(self, exon, gene_id, is_cds, contains_start, contains_stop, cds_overlap=0.0):
        self.exon = exon
        self.gene_ids = {gene_id}
        self.is_cds = is_cds
        self.contains_start = contains_start
        self.contains_stop = contains_stop
        if self.is_cds:
            self.cds_overlap = 1.0
        else:
            self.cds_overlap = cds_overlap
        self.whole_codon_count = ((self.exon[1] - self.exon[0] + 1) % 3) == 0

    def merge(self, other):
        self.gene_ids.update(other.gene_ids)
        self.is_cds = self.is_cds or other.is_cds
        self.contains_start = self.contains_start or other.contains_start
        self.contains_stop = self.contains_stop or other.contains_stop
        self.cds_overlap = max(self.cds_overlap, other.cds_overlap)


def process_gene(gene_db, gene):
    cds = set()
    start_codons = set()
    stop_codons = set()
    if isinstance(gene, str):
        g = gene_db[gene]
    else:
        g = gene
    for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):

        for s in gene_db.children(t, featuretype='start_codon', order_by='start'):
            start_codons.add((s.start, s.end))
        for s in gene_db.children(t, featuretype='stop_codon', order_by='start'):
            stop_codons.add((s.start, s.end))
        for c in gene_db.children(t, featuretype='CDS', order_by='start'):
            cds.add((c.start, c.end))

    return cds, start_codons, stop_codons


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
        cds, start_codons, stop_codons = gene_dicts[gene_id]

        is_cds = equals_any_feature(exon, cds)
        has_start = contains_any_feature(exon, start_codons)
        has_stop = contains_any_feature(exon, stop_codons)
        cds_overlap = overlaps_any_feature(exon, cds) if not is_cds else 1.0

        if exon_id not in exon_dict:
            exon_dict[exon_id] = ExonInfo(exon, gene_id, is_cds, has_start, has_stop, cds_overlap)
        else:
            exon_dict[exon_id].merge(ExonInfo(exon, gene_id, is_cds, has_start, has_stop, cds_overlap))

    return exon_dict


def count_stats(exon_info_dict):
    exon_count = 0
    cds_count = 0
    overlaps_cds = 0
    whole_codon_count = 0
    whole_codon_cds_count = 0
    start_count = 0
    stop_count = 0
    avg_cds_overlap = 0.0

    for exon_id in exon_info_dict.keys():
        exon_info = exon_info_dict[exon_id]
        exon_count += 1
        if exon_info.is_cds:
            cds_count += 1
        if exon_info.whole_codon_count:
            whole_codon_count += 1
            if exon_info.is_cds:
                whole_codon_cds_count += 1
        if exon_info.contains_start:
            start_count += 1
        if exon_info.contains_stop:
            stop_count += 1
        if exon_info.cds_overlap > 0:
            overlaps_cds += 1
            if exon_info.cds_overlap == 1 and not exon_info.is_cds:
                cds_count += 1
        avg_cds_overlap += exon_info.cds_overlap

    avg_cds_overlap /= exon_count
    return exon_count, cds_count, overlaps_cds, whole_codon_count, whole_codon_cds_count, start_count, stop_count, avg_cds_overlap

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name", required=True)
    parser.add_argument("--genedb", type=str, help="gffutils gene DB", required=True)
    parser.add_argument("--input", "-i", type=str, nargs='+',
                        help="one or more files/dirs with input exon lists", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gene_db = gffutils.FeatureDB(args.genedb)
    gene_dicts = {}

    with open(args.output, "w") as outf:
        outf.write("#name\texons\tcds\toverlaps_cds\t3Xbp\t3Xbp_CDS\tstart_codons\tstop_codons\tavg_cds_overlap\n")
        for inf in args.input:
            if os.path.isdir(inf):
                for f in glob.glob(os.path.join(inf, "*")):
                    if os.path.isfile(f):
                        name = os.path.basename(f)
                        exon_info_dict = process_exons(load_exon_pairs(f), gene_db, gene_dicts)
                        outf.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n" % (name, *count_stats(exon_info_dict)))
            elif os.path.isfile(inf):
                name = os.path.basename(inf)
                exon_info_dict = process_exons(load_exon_pairs(inf), gene_db, gene_dicts)
                outf.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n" % (name, *count_stats(exon_info_dict)))




if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
