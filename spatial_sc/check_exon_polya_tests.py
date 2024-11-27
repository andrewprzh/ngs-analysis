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
from collections import defaultdict
import numpy
from pyfaidx import Fasta

DELTA = 0


base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverse_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


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


def equal_ranges(range1, range2, delta=0):
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta


def equal_pos(pos1, pos2, delta=0):
    return abs(pos1 - pos2) <= delta


def get_geneid_map(genedb):
    gene_id2name = {}
    for g in genedb.features_of_type('gene'):
        g_id = g.id.split('.')[0]
        gene_id2name[g_id] = g.id
    return gene_id2name


def get_transcripts(genedb, gene):
    transcripts = defaultdict(list)
    if isinstance(gene, str):
        g = genedb[gene]
    else:
        g = gene

    for t in genedb.children(g, featuretype=('transcript', 'mRNA')):
        for c in genedb.children(t, featuretype='exon', order_by='start'):
            transcripts[t.id].append((c.start, c.end))

    return transcripts


def process_test(polya1, polya2, exon, transcripts, strand):
    combinations = set()
    # print(polya1, polya2, exon, strand)
    diff = abs(polya1 - polya2)
    p1_a = False
    p2_a = False
    for t in transcripts:
        exons = transcripts[t]
        # print(exons, strand)
        polya_site = exons[-1][1] if strand == '+' else exons[0][0]
        has_polya1 = equal_pos(polya1, polya_site, 1)
        has_polya2 = equal_pos(polya2, polya_site, 1)
        if polya1 == polya_site: p1_a = True
        if polya2 == polya_site: p2_a = True
        if not has_polya1 and not has_polya2: continue

        contains_exon = any(equal_ranges(e, exon, DELTA) for e in exons)
        if has_polya1:
            combinations.add((polya1, contains_exon))
        if has_polya2:
            combinations.add((polya2, contains_exon))

#    print("Combinations %d %s" % (len(combinations), str(combinations)))
#    if len(combinations) < 2 or (len(combinations) == 2 and len(set([x[0] for x in combinations])) == 1):
#        print("+++ ODD POLYA! +++")
    return len(combinations), diff, p1_a, p2_a


# chr2_152694389_152694389_-,chr2_152672551_152672551_-,chr2_152714550_152714629_ENSG00000196504.19_-
def process_csv(inf, genedb, gene_id_dict, chr_dict, args):
    combination_stats = defaultdict(int)
    close_polya = defaultdict(int)
    diffs = defaultdict(int)
    genomic_strs = {}
    for l in open(inf):
        if l.startswith("PolyA"): continue
        v = l.strip().split(',')
        assert len(v) == 3
        polya1 = int(v[0].split('_')[1])
        polya2 = int(v[1].split('_')[1])
        exon_values = v[2].split('_')
        exon = (int(exon_values[1]), int(exon_values[2]))
        chr_id = exon_values[0]
        gene_id = gene_id_dict[exon_values[3].split('.')[0]]
        strand = exon_values[4]

        transcripts = get_transcripts(genedb, gene_id)
        comb, diff, p1_a, p2_a = process_test(polya1, polya2, exon, transcripts, strand)

        if comb == 4 and diff >= 10:
            #print(l.strip())
            diffs[diff] += 1
        for c in [1, 3, 6, 10, 20, 50]:
            if diff <= c:
                close_polya[(c, p1_a and p2_a)] += 1
        combination_stats[comb] += 1

        polya_min = min(polya1, polya2)
        polya_max = max(polya1, polya2)
        polya_region = (chr_id, polya_min, polya_max)

        polya_percent = 0
        if chr_dict and polya_region not in genomic_strs and diff < 50:
            seq = str(chr_dict[chr_id][polya_min:polya_max+1])
            if strand == '-':
                seq = reverse_complement(seq)
            genomic_strs[polya_region] = seq
            polya_percent = seq.count('A') / len(seq)

        # ALL CONDITIONS ARE HERE
        # comb - number of possible combinations
        # p1_a, p2_a - whether polya1 and polya2 are annotated
        # diff - distance
        # polya_percent - percent of A bases inbetween these polyas
        if comb >= args.min_combinations and diff >= args.min_distance and polya_percent <= args.max_a_fraction and (not args.annotated or (p1_a and p2_a)):
            print(l.strip())

    return combination_stats, genomic_strs


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--output", "-o", type=str, help="output dir", required=True)
    parser.add_argument("--genedb", type=str, help="gffutils gene DB", required=True)
    parser.add_argument("--reference", type=str, help="genome in FASTA format with .fai index (samtools faidx genome.fa)")
    parser.add_argument("--csv", type=str, help="CSV file with polyA1/polyA2/exon (Lieke format)", required=True)
    parser.add_argument("--delta", type=int, help="coordinate comparison delta", default=0)
    parser.add_argument("--max_a_fraction", type=float, help="maximum allowed fraction of A bases between "
                                                             "polyA sites closer than 50bp (requires reference genome) "
                                                             "[0.5]", default=0.5)
    parser.add_argument("--min_combinations", type=int, help="minimal number of possible polyA-exon combinations [4]", default=4)
    parser.add_argument("--min_distance", type=int, help="minimal distance between polyA sites [10]", default=10)
    parser.add_argument("--annotated", type=int, help="coordinate comparison delta", default=0)


    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    global DELTA
    DELTA = args.delta
    genedb = gffutils.FeatureDB(args.genedb)
    reference_record_dict = Fasta(args.reference, indexname=args.reference + ".fai") if args.reference else None
    gene_ids = get_geneid_map(genedb)
    combination_stats, genomic_strs = process_csv(args.csv, genedb, gene_ids, reference_record_dict, args)

    if genomic_strs:
        a_percent = []
        lens = []
        for k in genomic_strs:
            s = genomic_strs[k]
            #print(len(s), "\t", s.count("A") / len(s))
            if s.count("A") * 2 <= len(s):
                a_percent.append(s.count("A") / len(s))

            lens.append(len(s))

#        print(numpy.histogram(a_percent, bins=[0.1 * i for i in range(11)])[0])
#        print(numpy.histogram(lens, bins=[ i for i in range(5, 51)])[0])
#        print(len(a_percent))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
