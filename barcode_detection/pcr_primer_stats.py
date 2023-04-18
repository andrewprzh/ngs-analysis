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

    def __init__(self):
        self.pcr_primer_indexer = KmerIndexer([PRCDetector.PCR_PRIMER], kmer_size=6)
        self.ups_primer_indexer = KmerIndexer([PRCDetector.UPS_PRIMER], kmer_size=6)

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

    def _detect_exact_positions(self, sequence, start, end, kmer_size, pattern, pattern_occurrences,
                                min_score=0, start_delta=-1, end_delta=-1):
        if not pattern_occurrences:
            return None, None
        potential_start = start + min(pattern_occurrences[pattern][2])
        potential_start = max(start, potential_start - kmer_size)
        potential_end = start + max(pattern_occurrences[pattern][2])
        potential_end = min(end, potential_end + 2 * kmer_size)
        start_pos, end_pos, pattern_start, pattern_end = \
            align_pattern_ssw(sequence, potential_start, potential_end, pattern, min_score)

        if start_pos is None:
            return None, None

        if start_delta > 0 and pattern_start > start_delta:
            return None, None
        if end_delta > 0 and len(pattern) - pattern_end - 1 > end_delta:
            return None, None
        return start_pos, end_pos

    def _find_barcode_umi_fwd(self, read_id, sequence):
        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:polyt_start + 1])
        primer_start, primer_end = self._detect_exact_positions(sequence, 0, polyt_start + 1,
                                                                self.pcr_primer_indexer.k, self.PCR_PRIMER,
                                                                primer_occurrences, min_score=5,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is None:
            return BarcodeDetectionResult(read_id, polyt_start)
        logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))

        linker_occurrences = self.linker_indexer.get_occurrences(sequence[primer_end + 1:polyt_start + 1])
        linker_start, linker_end = self._detect_exact_positions(sequence, primer_end + 1, polyt_start + 1,
                                                                self.linker_indexer.k, self.LINKER,
                                                                linker_occurrences, min_score=11,
                                                                start_delta=self.TERMINAL_MATCH_DELTA,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)
        if linker_start is None:
            return BarcodeDetectionResult(read_id, polyt_start, primer_end)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        potential_barcode = sequence[primer_end + 1:linker_start] + \
                            sequence[linker_end + 1:linker_end + self.RIGHT_BC_LENGTH + 2]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode, min_score=12)
        if barcode is None:
            return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = primer_end + 1 + (linker_end - linker_start + 1) + bc_end + 1
        potential_umi_end = polyt_start - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        logger.debug("Potential UMI: %s" % potential_umi)

        umi = None
        good_umi = False
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_score, umi_start, umi_end = \
                find_candidate_with_max_score_ssw(matching_umis, potential_umi, min_score=7)
            logger.debug("Found UMI %s %d-%d" % (umi, umi_start, umi_end))

        if not umi :
            umi = potential_umi
            if self.UMI_LENGTH - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LENGTH + self.UMI_LEN_DELTA and \
                    all(x != "T" for x in umi[-self.NON_T_UMI_BASES:]):
                good_umi = True

        if not umi:
            return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end, barcode, bc_score)
        return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end,
                                      barcode, bc_score, umi, good_umi)


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
