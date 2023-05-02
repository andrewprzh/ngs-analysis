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
import logging

from kmer_indexer import KmerIndexer
from common import *
from reports import *
from detect_double_barcodes import BarcodeCaller

logger = logging.getLogger('BarcodeCaller')


class PRCDetector:
    UPS_PRIMER = "AAGCAGTGGTATCAACGCAGAGT"
    PCR_PRIMER = "CTACACGACGCTCTTCCGATCT"
    UPS_PRIMER_REV = reverese_complement(UPS_PRIMER)
    TERMINAL_MATCH_DELTA = 3

    def __init__(self):
        self.pcr_primer_indexer = KmerIndexer([PRCDetector.PCR_PRIMER], kmer_size=6)
        self.ups_primer_indexer = KmerIndexer([PRCDetector.UPS_PRIMER], kmer_size=6)
        self.ups_primer_rev_indexer = KmerIndexer([PRCDetector.UPS_PRIMER_REV], kmer_size=6)


    def find_barcode_umi(self, read_id, sequence):
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id, sequence):
        ups_occurrences = self.ups_primer_rev_indexer.get_occurrences(sequence)
        ups_start, ups_end = detect_exact_positions(sequence, 0, len(sequence),
                                                    self.ups_primer_rev_indexer.k, self.UPS_PRIMER_REV,
                                                    ups_occurrences, min_score=15,
                                                    start_delta=self.TERMINAL_MATCH_DELTA)

        if not ups_start:
            ups_end = -1
            ups_start = len(sequence) - 1

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence)
        primer_start, primer_end = detect_exact_positions(sequence, 0, ups_start,
                                                          self.pcr_primer_indexer.k, self.PCR_PRIMER,
                                                          primer_occurrences, min_score=15,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if not primer_start:
            primer_start = -1

        ups1_start = -1
        if primer_start == -1:
            ups1_occurrences = self.ups_primer_indexer.get_occurrences(sequence)
            ups1_start, ups1_end = detect_exact_positions(sequence, 0, ups_start,
                                                          self.ups_primer_indexer.k, self.UPS_PRIMER,
                                                          ups1_occurrences, min_score=15,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
            if not ups1_start:
                ups1_start = -1

        return BarcodeDetectionResult(read_id, ups_end, primer_start, ups1_start)


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM)", required=True)
    parser.add_argument("--output", "-o", type=str, help="output file", required=True)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    set_logger(logger)
    barcode_detector = PRCDetector()
    barcode_caller = BarcodeCaller(args.output, barcode_detector)
    barcode_caller.process(args.input)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
