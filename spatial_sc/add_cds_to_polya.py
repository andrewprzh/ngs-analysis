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


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


def find_cds(gffutils_db, t_id, polya_pos):
    t = gffutils_db[t_id]
    cds = []
    exons = []
    for e in gffutils_db.children(t, featuretype="CDS", order_by="start"):
        cds.append((e.start, e.end))
    for e in gffutils_db.children(t, featuretype="exons", order_by="start"):
        cds.append((e.start, e.end))

    if not cds:
        return "NC"

    utr = []
    if t.strand == "+":
        cds_end_pos = cds[-1][1]
        for i in range(len(exons)):
            e = exons[i]
            if e[0] <= cds_end_pos <= e[1]:
                utr.append((cds_end_pos + 1, e[1]))
            elif cds_end_pos <= e[0]:
                utr.append(e)

        utr[-1] = (utr[-1][0], polya_pos)
    else:
        cds_end_pos = cds[0][0]
        for i in range(len(exons)):
            e = exons[i]
            if e[0] <= cds_end_pos <= e[1]:
                utr.append((e[0], cds_end_pos - 1))
                break
            elif e[1] <= cds_end_pos:
                utr.append(e)

        utr[0] = (polya_pos, utr[0][1])
        
    utr = sorted(utr)
    utr_str = ";%;" + ";%;".join("%s_%d_%d_%s" % (t.seqid, e[0], e[1], t.strand) for e in utr)

    return "%s_%d_%d_%s" % (t.seqid, cds_end_pos, cds_end_pos, t.strand), str(intervals_total_length(utr)), utr_str


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", required=True)
    parser.add_argument("--input", "-i", type=str, help="input TSV", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    random.seed(args.seed)
    gffutils_db = gffutils.FeatureDB(args.genedb)
    cds_dict = {}
    with open(args.output, "w") as outf:
        for l in open(args.input):
            if l.startswith("#"):
                outf.write(l)
                continue
            l = l.strip()
            v = l.strip().split('\t')
            t_id = v[1]
            polya_pos = int(v[4].split("_")[1])
            if t_id in cds_dict:
                cds, utr_len, utrs = cds_dict[t_id]
            else:
                cds, utr_len, utrs = find_cds(gffutils_db, t_id, polya_pos)
                cds_dict[t_id] = (cds, utr_len, utrs)
            outf.write(l + "\t%s\t%s\t%s\n" % (cds, utr_len, utrs))



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
