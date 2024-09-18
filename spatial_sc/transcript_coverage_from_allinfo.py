#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
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


# % Transcript covered


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


def left_of(range1, range2):
    return range1[1] < range2[0]


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def transcript_coverage_fraction(read_range_list, isoform_range_list):
    intersection = 0
    pos1 = 0
    pos2 = 0

    while pos1 < len(read_range_list) and pos2 < len(isoform_range_list):
        block1 = read_range_list[pos1]
        block2 = isoform_range_list[pos2]
        if overlaps(block1, block2):
            intersection += min(block1[1], block2[1]) - max(block1[0], block2[0]) + 1
            if block2[1] < block1[1]:
                pos2 += 1
            else:
                pos1 += 1
        elif left_of(block2, block1):
            pos2 += 1
        else:
            pos1 += 1

    transcript_length = intervals_total_length(isoform_range_list)
    return float(intersection) / float(transcript_length)


def load_transcripts(genedb):
    gffutils_db = gffutils.FeatureDB(genedb)
    transcript_dict = {}
    print("Loading genedb from %s" % genedb)
    for t in gffutils_db.features_of_type(('transcript', 'mRNA')):
        exons = []
        for e in gffutils_db.children(t, order_by='start', featuretype="exon"):
            exons.append((e.start, e.end))
        transcript_dict[t.id] = (t.seqid, t.strand, exons)
    return transcript_dict


# 3f9af4d7-aea6-4be6-9a98-a9e924e7ae65    ENSG00000228463.10      None    ATTCGCCGTGGCGA  GCGGCAACG       ;%;     chr1_259027_259027_-    chr1_258545_258545_-    ;%;chr1_258545_259027_- novel   0       ENST00000450734.1       transcribed_processed_pseudogene
def load_allinfo(inf):
    print("Loading allinfo from %s" % inf)
    # read_id -> (chr, strand, gene, isoform, exons)
    read_dict = {}

    for l in open(inf):
        v = l.strip().split("\t")
        gene_id = v[1]
        transcript_id = v[11]
        exons = v[8].split(";%;")
        exon1 = exons[1].split("_")
        chr_id = exon1[0]
        strand = exon1[3]
        exon_coords = []
        for e in exons:
            if not e: continue
            ex = e.split("_")
            exon_coords.append((int(ex[1]), int(ex[2])))
        read_dict[v[0]] = (chr_id, strand, gene_id, transcript_id, exon_coords)

    print("Loaded %d reads" % len(read_dict))
    return read_dict


def load_genes(inf):
    print("Loading genes from %s" % inf)
    gene_set = set()
    for l in open(inf):
        v = l.strip().split("\t")
        if not v[0]: continue
        gene_set.add(v[0])
    return gene_set


def process_allinfo(read_dict, transcript_dict, gene_set=None):
    transcript_cov_fractions = []
    for readid in read_dict:
        read_info = read_dict[readid]
        gene_id = read_info[2]
        if gene_set and gene_id not in gene_set: continue
        transcript_id = read_info[3]
        isoform_exons = transcript_dict[transcript_id][2]
        read_exons = read_info[4]
        transcript_cov_fractions.append(transcript_coverage_fraction(read_exons, isoform_exons))
    return transcript_cov_fractions


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--allinfo", "-a", nargs='+', type=str, help="allinfo files", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    parser.add_argument("--gene_list", nargs="+", type=str, help="gene list(s) to analyse")
    args = parser.parse_args()
    return args


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
    transcript_dict = load_transcripts(args.genedb)

    for allinfo in args.allinfo:
        read_dict = load_allinfo(allinfo)
        transcript_cov_fractions = process_allinfo(read_dict, transcript_dict)
        print("All reads for %s" % allinfo)
        print("Mean\t%.2f" % numpy.mean(transcript_cov_fractions))
        print("Med\t%.2f" % numpy.median(transcript_cov_fractions))
        print("Q25\t%.2f" % numpy.quantile(transcript_cov_fractions, 0.25))
        print("Q75\t%.2f" % numpy.quantile(transcript_cov_fractions, 0.75))
        for gene_list in args.gene_lists:
            gene_set = load_genes(gene_list)
            transcript_cov_fractions = process_allinfo(read_dict, transcript_dict, gene_set)
            print("Stats for %s" % gene_list)
            print("Mean\t%.2f" % numpy.mean(transcript_cov_fractions))
            print("Med\t%.2f" % numpy.median(transcript_cov_fractions))
            print("Q25\t%.2f" % numpy.quantile(transcript_cov_fractions, 0.25))
            print("Q75\t%.2f" % numpy.quantile(transcript_cov_fractions, 0.75))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
