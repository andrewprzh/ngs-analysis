############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from collections import defaultdict, namedtuple
from traceback import print_exc
import gffutils
import random


TSSRecord = namedtuple("TSSRecord", ("chr", "coordinate", "strand", "isoform_id"))


def read_expression(expr_f):
    expr_dict = {}
    for l in open(expr_f):
        if l.startswith("#"):
            continue
        tokens = l.strip().split()
        if len(tokens) < 2:
            continue
        expr_dict[tokens[0]] = (int(tokens[1]), float(tokens[2]))
    return expr_dict


def generate_tss(genedb):
    db = gffutils.FeatureDB(genedb, keep_order=True)
    tss_list = []
    for t in db.features_of_type('transcript'):
        exons = list(db.children(t, featuretype='exon', order_by='start'))
        chr = t.seqid
        strand = t.strand
        if strand == '+':
            coordinate = exons[0][0]
        else:
            coordinate = exons[-1][1]
        tss_list.append(TSSRecord(chr, coordinate, strand, t.id))
    return tss_list


def generate_cage_intervals(args):
    expr_dict = read_expression(args.expr)
    tss_list = generate_tss(args.genedb)
    cage_peaks = {}

    i = 0
    prev_c
    while i < len(tss_list):

        t_id = tss_record
        if t_id not in expr_dict or expr_dict[t_id][0] < args.min_count or expr_dict[t_id][1] < args.min_tpm:
            continue

        if random.random() < args.skip_probability / expr_dict[t_id][1]:
            continue

        interval_len = random.randint(args.min_interval_len, args.max_interval_len)
        left_delta = random.randint(1, interval_len - 1)
        tss = tss_dict[t_id][2]
        cage_interval = (tss - left_delta, tss - left_delta + interval_len)
        cage_peaks



def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expr", "-e", help="isoform expression", type=str)
    parser.add_argument("--output", "-o", help="output BED file", type=str)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--min_tpm", help="minimal TPM value", type=int, default=0)
    parser.add_argument("--min_count", help="minimal count", type=int, default=0)
    parser.add_argument("--skip_probability", help="probability coefficient of skipping TSS (divided by TPM)",
                        type=float, default=0.0)
    # parser.add_argument("--decoy_probability", help="probability of random peak to be inserted",
    #                     type=float, default=0.0)
    parser.add_argument("--min_interval_len", help="minimal CAGE interval length", type=int, default=5)
    parser.add_argument("--max_interval_len", help="maximal CAGE interval length", type=int, default=50)

    args = parser.parse_args()

    if args.genedb is None or args.expr is None or args.output is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)