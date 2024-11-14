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
from Bio import SeqIO
from pybedtools.cbedtools import defaultdict
import numpy


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


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


class ExonInfo:
    def __init__(self, exon, gene_id, is_cds, contains_start, contains_stop, preceding_introns, subsequent_introns, cds_overlap=0.0):
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
        self.preceding_introns = set(preceding_introns)
        self.subsequent_introns = set(subsequent_introns)

    def merge(self, other):
        self.gene_ids.update(other.gene_ids)
        self.is_cds = self.is_cds or other.is_cds
        self.contains_start = self.contains_start or other.contains_start
        self.contains_stop = self.contains_stop or other.contains_stop
        self.cds_overlap = max(self.cds_overlap, other.cds_overlap)
        self.preceding_introns.update(other.preceding_introns)
        self.subsequent_introns.update(other.subsequent_introns)


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

        exon_list = []
        for e in gene_db.children(t, featuretype='exon', order_by='start'):
            exon_list.append((e.start, e.end))

        seqid = t.seqid
        strand = t.strand
        for i in range(len(exon_list)):
            exon = exon_list[i]
            is_cds = equals_any_feature(exon, cds)
            has_start = contains_any_feature(exon, start_codons)
            has_stop = contains_any_feature(exon, stop_codons)
            cds_overlap = overlaps_any_feature(exon, cds) if not is_cds else 1.0
            if strand == "+":
                preceding_intron = [] if i == 0 else junctions_from_blocks([exon_list[i-1], exon])
                subsequent_intron = [] if i == len(exon_list)-1 else junctions_from_blocks([exon, exon_list[i+1]])
            else:
                subsequent_intron = [] if i == 0 else junctions_from_blocks([exon_list[i-1], exon])
                preceding_intron = [] if i == len(exon_list)-1 else junctions_from_blocks([exon, exon_list[i+1]])

            exon_id = "%s_%d_%d_%s" % (seqid, exon[0], exon[1], strand)
            if exon_id not in exon_dict:
                exon_dict[exon_id] = ExonInfo(exon, g.id, is_cds, has_start, has_stop, preceding_intron, subsequent_intron, cds_overlap)
            else:
                exon_dict[exon_id].merge(ExonInfo(exon, g.id, is_cds, has_start, has_stop, preceding_intron, subsequent_intron, cds_overlap))

    return exon_dict


def load_exon_info(genedb_path):
    gene_db = gffutils.FeatureDB(genedb_path)
    chr_dicts = {}
    for g in gene_db.features_of_type('gene'):
        if g.seqid not in chr_dicts:
            chr_dicts[g.seqid] = {}
        exon_dict = chr_dicts[g.seqid]
        gene_dict = process_gene(gene_db, g)
        for e in sorted(gene_dict.keys()):
            if e in exon_dict:
                exon_dict[e].merge(gene_dict[e])
            else:
                exon_dict[e] = gene_dict[e]
    return chr_dicts


def gc_content(seq):
    gc = 0
    for c in seq:
        if c.upper() in ['C', 'G']:
            gc += 1
    return float(gc) / float(len(seq))


#>chr8_67163732_67163798_+_ENSG00000104218.15_downstream
def process_fasta(in_fasta, chr_dicts):
    total_exons = 0
    cds_count = 0
    introns_lengths = defaultdict(list)
    cds_fractions = []
    gc_content = []

    for r in SeqIO.parse(in_fasta, "fasta"):
        seq_name = r.id
        v = seq_name.split("_")
        chr_id = v[0]
        exon_id = "%s_%s_%s_%s" % (chr_id, v[1], v[2], v[3])
        exon_info = chr_dicts[chr_id][exon_id]
        total_exons += 1
        cds_fractions.append(exon_info.cds_overlap)
        if exon_info.is_cds: cds_count += 1
        gc_content.append(gc_content(str(r.seq)))

        preceding_intron_lengths = list(map(interval_len, exon_info.preceding_introns))
        subsequent_intron_lengths = list(map(interval_len, exon_info.subsequent_introns))
        introns_lengths["upstream_min"].append(min(preceding_intron_lengths))
        introns_lengths["upstream_max"].append(max(preceding_intron_lengths))
        introns_lengths["upstream_all"] += preceding_intron_lengths
        introns_lengths["downstream_min"].append(min(subsequent_intron_lengths))
        introns_lengths["downstream_max"].append(max(subsequent_intron_lengths))
        introns_lengths["downstream_all"] += subsequent_intron_lengths
    return total_exons, cds_count, cds_fractions, gc_content, introns_lengths


def print_hist(numpy_hist, stream):
    for i in range(len(numpy_hist[0])):
        stream.write("%.2f\t%d\n" % ((numpy_hist[1][i] + numpy_hist[1][i+1]) / 2, numpy_hist[0][i]))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils gene DB", required=True)
    parser.add_argument("--fasta", "-f", nargs='+', type=str, help="FASTA files with exons", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print("Loading %s" % args.genedb)
    chr_dicts = load_exon_info(args.genedb)

    for f in args.fasta:
        print("Processing %s" % f)
        total_exons, cds_count, cds_fractions, gc_content, introns_lengths = process_fasta(f, chr_dicts)
        name = os.path.splitext(os.path.basename(f))[0]
        with open(os.path.join(args.output, name) + ".tsv") as outf:
            outf.write("CDS count\t%d / %d\n" % (cds_count, total_exons))
            cds_cov_hist = numpy.histogram(cds_fractions, [0.1 * i for i in range(11)])
            outf.write("\nCDS coverage hist\t%.2f\t%.2f\n" % (numpy.mean(cds_fractions), numpy.median(cds_fractions)))
            print_hist(cds_cov_hist, outf)

            gc_hist = numpy.histogram(gc_content, [5 * i for i in range(21)])
            outf.write("\nGC content\t%.2f\t%.2f\n" % (numpy.mean(gc_content), numpy.median(gc_content)))
            print_hist(gc_hist, outf)

            for k in sorted(introns_lengths.keys()):
                intron_lens = introns_lengths[k]
                hist = numpy.histogram(intron_lens, [100 * i for i in range(11)] + [1000 * i for i in range(2, 11)] + [100000, 1000000])
                outf.write("\nIntron lengths %s\t%.2f\t%.2f\n" % (k, numpy.mean(intron_lens), numpy.median(intron_lens)))
                print_hist(hist, outf)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
