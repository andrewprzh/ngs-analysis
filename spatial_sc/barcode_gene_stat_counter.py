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
import logging



logger = logging.getLogger('GeneBarcodeStats')


class GeneBarcodeStats:
    def __init__(self, args):
        self.args = args
        self.max_edit_distance = args.min_distance
        self.barcode_dict = self.load_barcodes(args.barcodes)
        self.read_to_gene_dict = self.load_assignments(args.read_assignmetns)
        self.good_barcode_score = 13

    def load_barcodes(self, in_file):
        # afaf0413-4dd7-491a-928b-39da40d68fb3    99      56      66      83      +       GCCGATACGCCAAT  CGACTGAAG       13      True
        logger.info("Loading barcodes from " + in_file)
        barcode_dict = {}
        for l in open(in_file):
            v = l.strip().split("\t")
            if len(v) < 10: continue
            barcode = v[6]
            if barcode == "*":
                continue
            umi = v[7]
            score = int(v[8])
            trusted_umi = v[9]
            barcode_dict[v[0]] = (barcode, umi, score, trusted_umi)
        logger.info("Loaded %d barcodes " % len(barcode_dict))
        return barcode_dict

    def load_assignments(self, assignment_file):
        logger.info("Loading barcodes from " + assignment_file)
        read_to_gene_dict = defaultdict(list)
        assignment_count = 0
        # 1251f521-d4c2-4dcc-96d7-85070cc44e12    chr1    -       ENST00000477740.5       ENSG00000238009.6       ambiguous       ism_internal    112699-112804,120721-120759     PolyA=False
        for l in open(assignment_file):
            if l.startswith("#"): continue
            v = l.strip().split("\t")
            read_id = v[0]
            gene_id = v[4]
            assignment_type = v[5]
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))
            assignment_count += 1
            read_to_gene_dict[read_id].append((assignment_type, gene_id, len(exon_blocks)))
        logger.info("Loaded %d barcodes " % assignment_count)
        return read_to_gene_dict

    def count_stats(self):
        total_reads = len(self.read_to_gene_dict.keys())
        spliced = 0
        assigned = 0
        assigned_barcoded = 0
        assigned_good_barcoded = 0
        unique = 0
        unique_barcoded = 0
        unique_good_barcoded = 0
        unique_spliced = 0
        unique_spliced_barcoded = 0
        unique_spliced_good_barcoded = 0

        for read_id in self.read_to_gene_dict.keys():
            read_info = self.read_to_gene_dict[read_id]
            assignment_type = read_info[0][0]
            if len(read_info[0][2]) > 1:
                spliced += 1
            if assignment_type.startswith("noninformative"):
                continue
            barcode_info = self.barcode_dict[read_id] if read_id in self.barcode_dict else None
            assigned += 1
            if barcode_info:
                assigned_barcoded += 1
                if barcode_info[2] >= self.good_barcode_score:
                    assigned_good_barcoded += 1

            if len(read_info) > 1:
                continue
            unique += 1
            if barcode_info:
                unique_barcoded += 1
                if barcode_info[2] >= self.good_barcode_score:
                    unique_good_barcoded += 1
            if len(read_info[0][2]) <= 1:
                continue
            unique_spliced += 1
            if barcode_info:
                unique_spliced_barcoded += 1
                if barcode_info[2] >= self.good_barcode_score:
                    unique_spliced_good_barcoded += 1

        logger.info("Total reads %d , of them spliced %s" % (total_reads, spliced))
        logger.info("Assigned to any gene %d , of them barcoded %d %d" % (assigned, assigned_barcoded, assigned_good_barcoded))
        logger.info("Uniquely assigned %d , of them barcoded %d %d" % (unique, unique_barcoded, unique_good_barcoded))
        logger.info("Uniquely assigned and spliced %d , of them barcoded %d %d" %
                    (unique_spliced, unique_spliced_barcoded, unique_spliced_good_barcoded))
        print("\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d" % (total_reads, spliced, assigned, assigned_barcoded,
                                                                assigned_good_barcoded, unique, unique_barcoded,
                                                                unique_good_barcoded, unique_spliced,
                                                                unique_spliced_barcoded, unique_spliced_good_barcoded))


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", type=str, help="read - barcode - UMI table", required=True)
    parser.add_argument("--read_assignments", "-r", type=str, help="IsoQuant read assignments", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    stat_counter = GeneBarcodeStats(args)
    stat_counter.count_stats()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)