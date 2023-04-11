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
import gzip
from collections import defaultdict
from traceback import print_exc

import pysam
from Bio import SeqIO
import logging

from kmer_indexer import KmerIndexer
from common import *

logger = logging.getLogger('BarcodeCaller')


class DoubleBarcodeDetector:
    LINKER = ""
    PCR_PRIMER = "CTACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 6
    RIGHT_BC_LENGTH = 8
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    UMI_LENGTH = 12

    def __int__(self, joint_barcode_list, umi_list=None):
        self.pcr_primer_indexer = KmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=7)
        self.linker_indexer = KmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=7)
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)
        self.umi_set = set(umi_list)
        if umi_list:
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)

    def find_barcode_umi(self, sequence):
        barcode, umi = self._find_barcode_umi_fwd(sequence)
        if barcode is not None:
            return barcode, umi

        rev_seq = reverese_complement(sequence)
        return self._find_barcode_umi_fwd(rev_seq)

    def _detect_exact_positions(self, sequence, start, end, kmer_size, pattern, pattern_occurrences):
        if not pattern_occurrences:
            return None, None
        potential_start = min(pattern_occurrences[pattern][2])
        potential_start = max(0, potential_start - kmer_size)
        potential_end = max(pattern_occurrences[pattern][2])
        potential_end = min(len(sequence), potential_end + 2 * kmer_size)
        pattern_start, pattern_end = align_pattern(sequence, potential_start, potential_end, pattern)
        return start + pattern_start, start + pattern_end

    def _find_barcode_umi_fwd(self, sequence):
        polyt_start = find_polyt_start(sequence)
        if polyt_start == -1:
            return None, None

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:polyt_start + 1])
        primer_start, primer_end = self._detect_exact_positions(sequence, 0, polyt_start + 1,
                                                                self.pcr_primer_indexer.k, self.PCR_PRIMER,
                                                                primer_occurrences)
        if primer_start is None:
            return None, None

        linker_occurrences = self.linker_indexer.get_occurrences(sequence[primer_end + 1:polyt_start + 1])
        linker_start, linker_end = self._detect_exact_positions(sequence, primer_end + 1, polyt_start + 1,
                                                                self.linker_indexer.k, self.LINKER,
                                                                linker_occurrences)
        if linker_start is None:
            return None, None

        potential_barcode = sequence[primer_end + 1:linker_start] + \
                            sequence[linker_end + 1:linker_end + self.RIGHT_BC_LENGTH + 2]
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_start, bc_end = find_candidate_with_max_score(matching_barcodes, potential_barcode)
        if barcode is None:
            return None, None

        potential_umi_start = primer_end + (linker_start - linker_end + 1) + bc_end + 1
        potential_umi_end = polyt_start - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        umi = None, None
        if self.umi_list:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_start, umi_end = find_candidate_with_max_score(matching_umis, potential_umi)
        if not umi and len(potential_umi) == self.UMI_LENGTH:
            umi = potential_umi

        return barcode, umi


class BarcodeCaller:
    def __int__(self, , output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_table = output_table

    def process(self, input_file):
        logger.info("Processing " + input_file)
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

        handle = input_file
        if low_ext in ['.gz', '.gzip']:
            handle = gzip.open(input_file, "rt")
            input_file = fname
            fname, outer_ext = os.path.splitext(os.path.basename(input_file))
            low_ext = outer_ext.lower()

        if low_ext in ['.fq', '.fastq']:
            self._process_fastx(SeqIO.parse(handle, "fastq"))
        elif low_ext in ['.fa', '.fasta']:
            self._process_fastx(SeqIO.parse(handle, "fasta"))
        elif low_ext in ['.bam', '.sam']:
            self._process_bam(pysam.AlignmentFile(input_file, "rb"))
        else:
            logger.error("Unknown file format " + input_file)

    def _process_fastx(self, read_handler):
        for r in read_handler:
            read_id = r.id
            seq = str(r.seq)
            self._process_read(read_id, seq)

    def _process_bam(self, read_handler):
        for r in read_handler:
            read_id = r.query_name
            seq = r.query_sequence
            self._process_read(read_id, seq)

    def _process_read(self, read_id, read_sequence):
        self.barcode_detector.find_barcode_umi(read_sequence)




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

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
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
