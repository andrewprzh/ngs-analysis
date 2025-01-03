#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import gffutils
import numpy
from pyfaidx import Fasta


# Intron length distribution (internal v. not interal)
# terminal exon lengths
# introns w GTAG, etc.
# DS with 100% data, 90%, 80%, and so on


CANONICAL_FWD_SITES = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
CANONICAL_REV_SITES = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


class TranscriptCollector:
    def __init__(self, genedb, reference):
        gffutils_db = gffutils.FeatureDB(genedb)
        self.transcript_dict = {}
        print("Loading genedb from %s" % genedb)
        for t in gffutils_db.features_of_type(('transcript', 'mRNA')):
            exons = []
            for e in gffutils_db.children(t, order_by='start', featuretype="exon"):
                exons.append((e.start, e.end))
            self.transcript_dict[t.id] = (t.seqid, t.strand, exons, junctions_from_blocks(exons))
        print("Loading reference from %s" % reference)
        self.reference_record_dict = Fasta(reference)
        self.canonical_info = {}

    def get_lengths(self, transcript_id):
        assert transcript_id in self.transcript_dict
        internal_exons = []
        terminal_exons = []
        intron_lengths = []

        exons = self.transcript_dict[transcript_id][2]
        exon_count = len(exons)
        if exon_count == 1:
            return internal_exons, terminal_exons, intron_lengths

        for i in range(exon_count):
            exon = exons[i]
            length = exon[1] - exon[0] + 1
            if i == exon_count - 1 or i == 0:
                terminal_exons.append(length)
            else:
                internal_exons.append(length)

        for intron in self.transcript_dict[transcript_id][3]:
            length = intron[1] - intron[0] + 1
            intron_lengths.append(length)

        return internal_exons, terminal_exons, intron_lengths

    def check_sites_are_canonical(self, chr_id, introns, strand):
        fasta_record = self.reference_record_dict[chr_id]
        canonical_count = 0
        non_canonical_count = 0
        for intron in introns:
            if intron not in self.canonical_info:
                intron_left_pos = intron[0]
                intron_right_pos = intron[1]
                left_site = str(fasta_record[intron_left_pos:intron_left_pos+2])
                right_site = str(fasta_record[intron_right_pos - 1:intron_right_pos + 1])
                if strand == '+':
                    self.canonical_info[intron] = (left_site, right_site) in CANONICAL_FWD_SITES
                else:
                    self.canonical_info[intron] = (left_site, right_site) in CANONICAL_REV_SITES

            if self.canonical_info[intron]:
                canonical_count += 1
            else:
                non_canonical_count += 1

        return canonical_count, non_canonical_count

    def get_canonical(self, transcript_id):
        chr_id = self.transcript_dict[transcript_id][0]
        introns = self.transcript_dict[transcript_id][3]
        strand  = self.transcript_dict[transcript_id][1]
        return self.check_sites_are_canonical(chr_id, introns, strand)


def load_allinfo(inf):
    print("Loading allinfo from %s" % inf)
    # read_id -> (chr, strand, gene, isoform)
    read_dict = {}

    for l in open(inf):
        v = l.strip().split("\t")
        gene_id = v[1]
        transcript_id = v[11]
        exons = v[8].split(";%;")
        exon1 = exons[1].split("_")
        chr_id = exon1[0]
        strand = exon1[3]
        read_dict[v[0]] = (chr_id, strand, gene_id, transcript_id)

    print("Loaded %d reads" % len(read_dict))
    return read_dict


def process_allinfo(read_dict: dict, transcript_collector: TranscriptCollector):
    internal_exons_lengths = []
    terminal_exons_lengths = []
    introns_lengths = []
    non_canonical_dict = defaultdict(int)
    for readid in read_dict:
        transcript_id = read_dict[readid][3]
        internal_exons, terminal_exons, introns = transcript_collector.get_lengths(transcript_id)
        internal_exons_lengths += internal_exons
        terminal_exons_lengths += terminal_exons
        introns_lengths += introns
        canonical_count, non_canonical_count = transcript_collector.get_canonical(transcript_id)
        non_canonical_dict[non_canonical_count] += 1

    return internal_exons_lengths, terminal_exons_lengths, introns_lengths, non_canonical_dict


def count_genes(read_dict: dict, percentage: float):
    assert 0.0 <= percentage <= 1.0
    gene_set = defaultdict(int)
    for readid in read_dict:
        r = random.random()
        if r > percentage:
            continue
        gene_set[read_dict[readid][2]] += 1
    return sum(v >= 10 for v in gene_set.values())


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--allinfo", "-a", nargs='+', type=str, help="allinfo files", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    parser.add_argument("--fasta", "-f", type=str, help="reference genome in fasta", required=True)
    parser.add_argument("--seed", type=int, help="randomizer seed", default=11 )
    args = parser.parse_args()
    return args


MIN_FRAC = 0.1
ITERATIONS = 20
EXON_BINS = [20 * i for i in range(51)] + [10000]
LONG_EXON_BINS = [50 * i for i in range(51)] + [10000]
SHORT_INTRON_BINS = [5 * i for i in range(51)] + [1000000]
INTRON_BINS = [2000 * i for i in range(51)] + [1000000]


def print_hist(bins, val_lists, name):
    with open(name + ".tsv", "w") as outf:
        for i in range(len(bins) - 1):
            outf.write("%d\t%s\n" % (bins[i], "\t".join(map(str, [val[i] for val in val_lists]))))

    with open(name + ".norm.tsv", "w") as outf:
        total_counts = [sum(v) for v in val_lists]
        for i in range(len(bins) - 1):
            normalized_values = []
            for j in range(len(val_lists)):
                normalized_values.append(val_lists[j][i] / total_counts[j])
            outf.write("%d\t%s\n" % (bins[i], "\t".join(map(lambda x: ("%.4f" % x), normalized_values))))


def main():
    args = parse_args()
    random.seed(args.seed)

    transcript_collector = TranscriptCollector(args.genedb, args.fasta)
    for allinfo in args.allinfo:
        name = os.path.splitext(os.path.basename(allinfo))[0]
        read_dict = load_allinfo(allinfo)
        internal_exons_lengths, terminal_exons_lengths, introns_lengths, non_canonical_dict = process_allinfo(read_dict, transcript_collector)

        gene_count_dict = defaultdict(list)
        percentage = MIN_FRAC
        while percentage < 1.0:
            for i in range(ITERATIONS):
                gene_count_dict[percentage].append(count_genes(read_dict, percentage))
            percentage += MIN_FRAC


        print(numpy.median(internal_exons_lengths))
        print(numpy.median(terminal_exons_lengths))
        print(numpy.median(introns_lengths))
        print(numpy.mean(internal_exons_lengths))
        print(numpy.mean(terminal_exons_lengths))
        print(numpy.mean(introns_lengths))

        internal_hist = numpy.histogram(internal_exons_lengths, bins=EXON_BINS)
        terminal_hist = numpy.histogram(terminal_exons_lengths, bins=EXON_BINS)
        print_hist(terminal_hist[1], [internal_hist[0], terminal_hist[0]], args.output + name + ".exon_lengths")

        internal_hist = numpy.histogram(internal_exons_lengths, bins=LONG_EXON_BINS)
        terminal_hist = numpy.histogram(terminal_exons_lengths, bins=LONG_EXON_BINS)
        print_hist(terminal_hist[1], [internal_hist[0], terminal_hist[0]], args.output + name + ".long_exon_lengths")

        intron_hist = numpy.histogram(introns_lengths, bins=INTRON_BINS)
        print_hist(intron_hist[1], [intron_hist[0]], args.output + name + ".intron_lengths")

        intron_hist = numpy.histogram(introns_lengths, bins=SHORT_INTRON_BINS)
        print_hist(intron_hist[1], [intron_hist[0]], args.output + name + ".short-intron_lengths")

        with open(args.output + name + ".canonical_stats.tsv", "w") as outf:
            for k in sorted(non_canonical_dict.keys()):
                outf.write("%d\t%d\n" % (k, non_canonical_dict[k]))

        with open(args.output + name + ".gene_cov.tsv", "w") as outf:
            for p in sorted(gene_count_dict.keys()):
                quantiles = numpy.quantile(gene_count_dict[p], [0.5, 0.25, 0.75])
                outf.write("%.2f\t%s\n" % (p * 100, "\t".join(map(str, quantiles))))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
