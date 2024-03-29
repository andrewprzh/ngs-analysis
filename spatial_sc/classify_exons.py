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
from collections import OrderedDict


def contains_any_feature(exon, ordered_feature_list, delta=0):
    for f in ordered_feature_list:
        if exon[0] - delta <= f[0] and f[1] + delta <= exon[1]:
            return True
    return False


def equals_any_feature(exon, ordered_feature_list):
    for f in ordered_feature_list:
        if exon == f:
            return True
    return False


def overlaps_any_feature(exon, ordered_feature_list):
    ovlp = 0.0
    for f in ordered_feature_list:
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


def process_gene(gene_db, g):
    exon_dict = {}

    for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):
        cds = []
        start_codons = []
        stop_codons = []
        for s in gene_db.children(t, featuretype='start_codon', order_by='start'):
            start_codons.append((s.start, s.end))
        for s in gene_db.children(t, featuretype='stop_codon', order_by='start'):
            stop_codons.append((s.start, s.end))
        for c in gene_db.children(t, featuretype='CDS', order_by='start'):
            cds.append((c.start, c.end))

        for e in gene_db.children(t, featuretype='exon', order_by='start'):
            exon = (e.start, e.end)
            is_cds = equals_any_feature(exon, cds)
            has_start = contains_any_feature(exon, start_codons)
            has_stop = contains_any_feature(exon, stop_codons)
            cds_overlap = overlaps_any_feature(exon, cds) if not is_cds else 1.0

            exon_id = "%s_%d_%d_%s" % (e.seqid, e.start, e.end, e.strand)
            if exon_id not in exon_dict:
                exon_dict[exon_id] = ExonInfo(exon, g.id, is_cds, has_start, has_stop, cds_overlap)
            else:
                exon_dict[exon_id].merge(ExonInfo(exon, g.id, is_cds, has_start, has_stop, cds_overlap))

    return exon_dict


def load_exon_info(genedb_path):
    gene_db = gffutils.FeatureDB(genedb_path)
    chr_dicts = []
    exon_dict = OrderedDict()
    current_seqid = ""
    for g in gene_db.features_of_type('gene', order_by=('seqid', 'start')):
        if current_seqid != g.seqid:
            chr_dicts.append(exon_dict)
            exon_dict = OrderedDict()
            current_seqid = g.seqid

        gene_dict = process_gene(gene_db, g)
        for e in sorted(gene_dict.keys()):
            if e in exon_dict:
                exon_dict[e].merge(gene_dict[e])
            else:
                exon_dict[e] = gene_dict[e]
    return chr_dicts


def print_exon_info(dict_list, out_file_name):
    with open(out_file_name, "w") as outf:
        outf.write("#exon\tCDS\tstart\tstop\tCDS_overlap\twhole_codon_count\n")
        for d in dict_list:
            for e, info in d.items():
                outf.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%d\n" %
                           (e, list(info.gene_ids)[0], info.is_cds, info.contains_start, info.contains_stop,
                            info.cds_overlap, info.whole_codon_count))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name", required=True)
    parser.add_argument("--genedb", type=str, help="gffutils gene DB", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    dict_list = load_exon_info(args.genedb)
    print_exon_info(dict_list, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
