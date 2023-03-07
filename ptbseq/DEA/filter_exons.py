#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", type=str, help="input counts", required=True)
    parser.add_argument("--output", "-o", type=str, help="output file name")
    parser.add_argument("--cutoff", type=int, help="cutoff for total counts", default=20)
    parser.add_argument("--keep_constitutive", action="store_true", help="do not remove constitutive exons", default=False)
    args = parser.parse_args()
    return args


# == range operations ==
def filter_counts(inf, outf, args):
    exon_inc_counts = defaultdict(int)
    exon_exc_counts = defaultdict(int)
    for l in open(inf, "r"):
        v = l.strip().split()
        if len(v) < 4:
            continue
        exon_id = v[0]
        exon_inc_counts[exon_id] += int(v[1])
        exon_exc_counts[exon_id] += int(v[2])

    with open(outf, "w") as filtered_file:
        for l in open(inf, "r"):
            v = l.strip().split()
            if len(v) < 4:
                continue
            exon_id = v[0]
            if exon_inc_counts[exon_id] + exon_exc_counts[exon_id] < args.cutoff:
                continue
            if not args.keep_constitutive and min(exon_inc_counts[exon_id], exon_exc_counts[exon_id]) == 0:
                continue
            filtered_file.write(l)


def main():
    args = parse_args()
    basename, ext = os.path.splitext(args.input)
    if not args.output:
        args.output = basename + ".filtered_%d.tsv" % args.cutoff
    filter_counts(args.input, args.output, args)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



