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

DELTA = 0


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
        has_polya1 = equal_pos(polya1, polya_site, 50)
        has_polya2 = equal_pos(polya2, polya_site, 50)
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
def process_csv(inf, genedb, gene_id_dict):
    combination_stats = defaultdict(int)
    close_polya = defaultdict(int)
    diffs = defaultdict(int)
    for l in open(inf):
        if l.startswith("PolyA"): continue
        v = l.strip().split(',')
        assert len(v) == 3
        polya1 = int(v[0].split('_')[1])
        polya2 = int(v[1].split('_')[1])
        exon_values = v[2].split('_')
        exon = (int(exon_values[1]), int(exon_values[2]))
        gene_id = gene_id_dict[exon_values[3].split('.')[0]]
        strand = exon_values[4]

        transcripts = get_transcripts(genedb, gene_id)
        #print("=== %s ===" % gene_id)
        comb, diff, p1_a, p2_a = process_test(polya1, polya2, exon, transcripts, strand)
        if comb == 4 and diff >= 10:
            print(l.strip())
            diffs[diff] += 1
        for c in [1, 3, 6, 10, 20, 50]:
            if diff <= c:
                close_polya[(c, p1_a and p2_a)] += 1
        combination_stats[comb] += 1
        #print("<<< %s >>>" % gene_id)

    #print(close_polya)

#    print(sum(diffs.values()))
#    for k in sorted(diffs.keys()):
#        print(k, diffs[k])
    return combination_stats


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--output", "-o", type=str, help="output dir", required=True)
    parser.add_argument("--genedb", type=str, help="gffutils gene DB", required=True)
    parser.add_argument("--csv", type=str, help="CSV file with polyA1/polyA2/exon", required=True)
    parser.add_argument("--delta", type=int, help="coordinate comparison delta", default=0)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    global DELTA
    DELTA = args.delta
    genedb = gffutils.FeatureDB(args.genedb)
    gene_ids = get_geneid_map(genedb)
    process_csv(args.csv, genedb, gene_ids)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
