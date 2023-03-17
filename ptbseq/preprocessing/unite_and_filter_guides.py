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


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output name", required=True)
    parser.add_argument("--guide_tsvs", "-g", nargs='+', type=str, help="barcode to guides TSVs", required=True)
    parser.add_argument("--method", "-m", type=str, help="how to process tsvs", choices=['union', 'intersection', 'merge'], default='intersection')
    parser.add_argument("--keep_ambiguous", action='store_true', help="keep ambiguous barcodes", default=False)

    args = parser.parse_args()
    return args


def load_barcodes(inf):
    barcode2guide = defaultdict(set)
    for l in open(inf):
        v = l.strip().split()
        barcode2guide[v[0]].add(v[1])
    print("Loaded %d barcodes from %s" % (len(barcode2guide), inf))
    return barcode2guide


def union_barcodes(barcode2guide1, barcode2guide2):
    result = defaultdict(set)
    for bc in barcode2guide1:
        result[bc].update(barcode2guide1[bc])

    for bc in barcode2guide2:
        result[bc].update(barcode2guide2[bc])
    return result


def intersect_barcodes(barcode2guide1, barcode2guide2):
    result = defaultdict(set)
    for bc in barcode2guide1:
        if bc in barcode2guide2:
            result[bc].update(barcode2guide1[bc].intersection(barcode2guide2[bc]))

    return result


def merge_barcodes(barcode2guide1, barcode2guide2):
    result = defaultdict(set)
    for bc in barcode2guide1:
        if bc in barcode2guide2:
            result[bc].update(barcode2guide1[bc].intersection(barcode2guide2[bc]))
        else:
            result[bc].update(barcode2guide1[bc])
    for bc in barcode2guide2:
        if bc not in barcode2guide1:
            result[bc].update(barcode2guide2[bc])
    return result


def filter_barcodes(barcode2guide):
    result = {}
    for bc in barcode2guide:
        if len(barcode2guide[bc]) == 1:
            result[bc] = barcode2guide[bc]
    return result


func_dict = {'union': union_barcodes, 'intersection': intersect_barcodes, 'merge': merge_barcodes}


def main():
    #set_logger(logger)
    args = parse_args()
    resulting_dict = defaultdict(set)
    input_tsvs = []
    for inf in args.guide_tsvs:
        input_tsvs.append(load_barcodes(inf))

    func = func_dict[args.method]
    res = func(input_tsvs[0], input_tsvs[1])
    for i in range(2, len(input_tsvs)):
        res = func(res, input_tsvs[i])

    print("%s yields %d barcodes" % (args.method, len(res)))
    if not args.keep_ambiguous:
        res = filter_barcodes(res)

    print("Barcodes after filtering %d" % len(res))
    with open(args.output, "w") as outf:
        for bc in res:
            for guide in res[bc]:
                outf.write("%s\t%s\n" % (bc, guide))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



