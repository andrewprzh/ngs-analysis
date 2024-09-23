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
import gzip


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


# 9ba552b4-fb77-4f4b-8795-fdceaf13cf00    chr1    -       ENST00000488147.1       ENSG00000227232.5       inconsistent    alternative_structure_novel:14830-16139,tes_match_precise:4     14400-14829,16140-16224 gene_assignment=inconsistent; PolyA=False; Classification=novel_not_in_catalog;
def load_read_assignments(read_assignments):
    # read_id -> (chr, strand, gene, isoform, exons)
    read_dict = {}

    fname, outer_ext = os.path.splitext(os.path.basename(read_assignments))
    low_ext = outer_ext.lower()

    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(read_assignments, "rt")
    else:
        handle = open(read_assignments)

    for l in handle:
        if l.startswith("#"): continue
        v = l.split("\t")
        assigniment_type = v[5]
        if not assigniment_type.startswith("unique"): continue
        gene_id = v[4]
        transcript_id = v[3]
        chr_id = v[1]
        strand = v[2]
        exons = v[7].split(",")
        exon_coords = []
        for e in exons:
            if not e: continue
            ex = e.split("-")
            exon_coords.append((int(ex[0]), int(ex[1])))
        read_dict[v[0]] = (chr_id, strand, gene_id, transcript_id, exon_coords)
    return read_dict


def load_genes(inf):
    print("Loading genes from %s" % inf)
    gene_set = set()
    for l in open(inf):
        v = l.strip().split("\t")
        if not v[0]: continue
        gene_set.add(v[0])
    return gene_set


def process_allinfo(read_dict, transcript_dict, gene_set=None, spliced_only=False):
    transcript_cov_fractions = []
    transcript_exon_count = []
    gene_ids = set()
    for readid in read_dict:
        read_info = read_dict[readid]
        gene_id = read_info[2]
        if gene_set and gene_id not in gene_set: continue
        transcript_id = read_info[3]
        read_exons = read_info[4]
        if spliced_only and len(read_exons) == 1:
            continue
        gene_ids.add(gene_id)
        isoform_exons = transcript_dict[transcript_id][2]
        transcript_cov_fractions.append(transcript_coverage_fraction(read_exons, isoform_exons))
        transcript_exon_count.append(len(read_exons))
    return transcript_cov_fractions, transcript_exon_count, gene_ids


def common_unique_genes(read_dict1, read_dict2):
    gene_counts1 = defaultdict(int)
    gene_counts2 = defaultdict(int)
    for read_id in read_dict1:
        gene_counts1[read_dict1[read_id][2]] += 1
    for read_id in read_dict2:
        gene_counts2[read_dict2[read_id][2]] += 1

    gene_set1 = set(gene_counts1.keys())
    gene_set2 = set(gene_counts2.keys())
    unique1 = gene_set1.difference(gene_set2)
    unique2 = gene_set2.difference(gene_set1)
    return (("common", gene_set1.intersection(gene_set2)),
            ("unique1", unique1),
            ("unique1_10", set(filter(lambda x: gene_counts1[x] >= 10, unique1))),
            ("unique2", unique2),
            ("unique2_10", set(filter(lambda x: gene_counts2[x] >= 10, unique2))))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", nargs='+', type=str, help="allinfo/read assignments files", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    parser.add_argument("--gene_list", nargs="+", type=str, help="gene list(s) to analyse", default=[])
    parser.add_argument("--data_type", '-d', choices=['auto', 'allinfo', 'ra'], help="input data type", default="auto")
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


def print_stats(header, read_dict, transcript_dict, gene_set=None, spliced_only=False):
    transcript_cov_fractions, transcript_exon_count, gene_set = process_allinfo(read_dict, transcript_dict, gene_set, spliced_only)
    if not transcript_cov_fractions:
        print("No data for %s" % header)
        return
    print(header)
    print("Reads\t%d" % len(transcript_cov_fractions))
    print("Genes\t%d" % len(gene_set))
    print("Mean \t%.5f" % numpy.mean(transcript_cov_fractions))
    print("Med  \t%.5f" % numpy.median(transcript_cov_fractions))
    print("E.avg \t%.5f" % numpy.mean(transcript_exon_count))

    # print("Q25\t%.5f" % numpy.quantile(transcript_cov_fractions, 0.25))
    # print("Q75\t%.5f" % numpy.quantile(transcript_cov_fractions, 0.75))
    print()


def main():
    args = parse_args()
    transcript_dict = load_transcripts(args.genedb)
    read_dicts = []

    for inf in args.input:
        if args.data_type == "auto":
            fname, outer_ext = os.path.splitext(os.path.basename(inf))
            low_ext = outer_ext.lower()
            if low_ext == ".allinfo":
                read_dict = load_allinfo(inf)
            else:
                read_dict = load_read_assignments(inf)
        elif args.data_type == "allinfo":
            read_dict = load_allinfo(inf)
        else:
            read_dict = load_read_assignments(inf)

        read_dicts.append(read_dict)
        print_stats("All reads for %s" % inf, read_dict, transcript_dict)
        for gene_list in args.gene_list:
            gene_set = load_genes(gene_list)
            print_stats("Stats for %s" % os.path.basename(gene_list), read_dict, transcript_dict, gene_set)

        print_stats("SPLICED reads for %s" % inf, read_dict, transcript_dict, spliced_only=True)
        for gene_list in args.gene_list:
            gene_set = load_genes(gene_list)
            print_stats("SPLICED Stats for %s" % os.path.basename(gene_list), read_dict, transcript_dict, gene_set, spliced_only=True)

    if len(read_dicts) == 2:
        for header, gene_set in common_unique_genes(read_dicts[0], read_dicts[1]):
            print("Using gene set %s of size %d" % (header, len(gene_set)))
            print_stats("Stats for allinfo %s" % (os.path.basename(args.allinfo[0])), read_dicts[0], transcript_dict, gene_set)
            print_stats("Stats for allinfo %s" % (os.path.basename(args.allinfo[1])), read_dicts[1], transcript_dict, gene_set)

            print("Using gene set %s of size %d, SPLICED ONLY" % (header, len(gene_set)))
            print_stats("Stats for allinfo %s" % (os.path.basename(args.allinfo[0])), read_dicts[0], transcript_dict, gene_set, spliced_only=True)
            print_stats("Stats for allinfo %s" % (os.path.basename(args.allinfo[1])), read_dicts[1], transcript_dict, gene_set, spliced_only=True)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
