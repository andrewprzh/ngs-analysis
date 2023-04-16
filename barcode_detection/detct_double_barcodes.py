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
    LINKER = "CTCTTCAGCGTTCCC"
    PCR_PRIMER = "CTACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 6
    RIGHT_BC_LENGTH = 8
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    UMI_LENGTH = 12
    TERMINAL_MATCH_DELTA = 2

    def __init__(self, joint_barcode_list, umi_list=None):
        self.pcr_primer_indexer = KmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=5)
        logger.info("Barcodes %d" % len(joint_barcode_list))
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)
        self.umi_set = None
        if umi_list:
            self.umi_set =  set(umi_list)
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)

    def find_barcode_umi(self, sequence):
        barcode, umi = self._find_barcode_umi_fwd(sequence)
        if barcode is not None:
            return barcode, umi

        rev_seq = reverese_complement(sequence)
        return self._find_barcode_umi_fwd(rev_seq)

    def _detect_exact_positions(self, sequence, start, end, kmer_size, pattern, pattern_occurrences,
                                min_score=0, start_delta=-1, end_delta=-1):
        if not pattern_occurrences:
            return None, None
        logger.info(pattern)
        logger.info(pattern_occurrences)
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

    def _find_barcode_umi_fwd(self, sequence):
        polyt_start = find_polyt_start(sequence)
        if polyt_start == -1:
            return None, None

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:polyt_start + 1])
        primer_start, primer_end = self._detect_exact_positions(sequence, 0, polyt_start + 1,
                                                                self.pcr_primer_indexer.k, self.PCR_PRIMER,
                                                                primer_occurrences, min_score=5,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is None:
            return None, None
        logger.info("PRIMER: %d-%d" % (primer_start, primer_end))

        linker_occurrences = self.linker_indexer.get_occurrences(sequence[primer_end + 1:polyt_start + 1])
        linker_start, linker_end = self._detect_exact_positions(sequence, primer_end + 1, polyt_start + 1,
                                                                self.linker_indexer.k, self.LINKER,
                                                                linker_occurrences, min_score=11,
                                                                start_delta=self.TERMINAL_MATCH_DELTA,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)
        if linker_start is None:
            return None, None
        logger.info("LINKER: %d-%d" % (linker_start, linker_end))

        potential_barcode = sequence[primer_end + 1:linker_start] + \
                            sequence[linker_end + 1:linker_end + self.RIGHT_BC_LENGTH + 2]
        logger.info("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10)
        barcode, bc_start, bc_end = find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode)
        if barcode is None:
            return None, None
        logger.info("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = primer_end + 1 + (linker_end - linker_start + 1) + bc_end + 1
        potential_umi_end = polyt_start - 1
        logger.info("UMI coordinates: %d - %d" % (potential_umi_start, potential_umi_end))
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        logger.info("Potential UMI: %s" % potential_umi)
        umi = None
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_start, umi_end = find_candidate_with_max_score_ssw(matching_umis, potential_umi)
        if not umi and len(potential_umi) == self.UMI_LENGTH:
            umi = potential_umi

        return barcode, umi


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_file = open(output_table, "w")

    def __del__(self):
        self.output_file.close()

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
        logger.info("==== %s ====" % read_id)
        barcode, umi = self.barcode_detector.find_barcode_umi(read_sequence)
        if barcode:
            self.output_file.write("%s\t%s\t%s\n" % (read_id, barcode, umi))


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


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
    parser.add_argument("--barcodes", "-b", type=str, help="barcode whitelist", required=True)
    parser.add_argument("--umi", "-u", type=str, help="potential UMIs")
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM)", required=True)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    set_logger(logger)
    barcodes = load_barcodes(args.barcodes)
    barcode_detector = DoubleBarcodeDetector(barcodes)
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
