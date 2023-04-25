#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc

from ssw import AlignmentMgr
import editdistance
from collections import defaultdict
from numpy import histogram
from kmer_indexer import *


class BarcodeStatCounter:
    def __init__(self, barcodes):
        self.barcodes = barcodes
        self.barcode_stats = []
        self.barcode_indexer = KmerIndexer(barcodes, kmer_size=4)

    def _align(self, bc1, bc2):
        align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
        align_mgr.set_read(bc1)
        align_mgr.set_reference(bc2)
        alignment = align_mgr.align(gap_open=1, gap_extension=1)
        return alignment.optimal_score

    def count_edit_distances(self):
        for i in range(len(self.barcodes)):
            max_score = 0
            best_ed = 100
            best_pair = None
            if i % 100 == 0:
                sys.stdout.write("Processed %d barcodes\r" % i)

            bc1 = self.barcodes[i]
            matching_barcodes = self.barcode_indexer.get_occurrences(bc1, hits_delta=7, ignore_equal=True)
            for bc2 in matching_barcodes.keys():
                if bc1 == bc2:
                    continue
                score = self._align(bc1, bc2)
                if score >= max_score:
                    ed = editdistance.eval(bc1, bc2)
                    if ed < best_ed:
                        best_pair = bc2
                        max_score = score
                        best_ed = ed

            self.barcode_stats.append(best_ed)

    def print_hist(self):
        bins = [i for i in range(20)]
        h, b = histogram(self.barcode_stats, bins)
        print(h)


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--barcodes", "-b", type=str, help="barcode list", required=True)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    barcodes = load_barcodes(args.barcodes)
    stat_counter = BarcodeStatCounter(barcodes)
    stat_counter.count_edit_distances()
    stat_counter.print_hist()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
