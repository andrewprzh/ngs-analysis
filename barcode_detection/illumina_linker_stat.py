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
import gzip
from traceback import print_exc
import itertools
import shutil
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

import pysam
from Bio import SeqIO
import logging

from kmer_indexer import KmerIndexer
from common import *
from reports import *

logger = logging.getLogger('BarcodeCaller')


READ_CHUNK_SIZE = 10000


class DoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    KMER_SIZE = 5

    TERMINAL_MATCH_DELTA = 1
    STRICT_TERMINAL_MATCH_DELTA = 0

    def __init__(self, joint_barcode_list):
        self.pcr_primer_indexer = KmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=5)
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)

        self.barcodes = set(joint_barcode_list)
        self.left_barcodes = set()
        self.right_barcodes = set()
        for b in joint_barcode_list:
            self.left_barcodes.add(b[:self.LEFT_BC_LENGTH])
            self.right_barcodes.add(b[self.LEFT_BC_LENGTH:])
        self.left_kmer_dist = defaultdict(int)
        self.right_kmer_dist = defaultdict(int)

    def find_barcode_umi(self, read_id, sequence):
        f_linker, f_left, f_right, f_exact, f_bc = self._find_barcode_umi_fwd(read_id, sequence)
        rev_seq = reverese_complement(sequence)
        r_linker, r_left, r_right, r_exact, r_bc = self._find_barcode_umi_fwd(read_id, rev_seq)
        if f_linker:
            return f_linker, f_left, f_right, f_exact, f_bc
        return r_linker, r_left, r_right, r_exact, r_bc

    def _find_barcode_umi_fwd(self, read_id, sequence):
        # look for linker in the entire read
        linker_occurrences = self.linker_indexer.get_occurrences(sequence)
        linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                          self.linker_indexer.k, self.LINKER,
                                                          linker_occurrences, min_score=14,
                                                          start_delta=self.TERMINAL_MATCH_DELTA,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return False, None, None, None, None
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        barcode_start = linker_start - self.LEFT_BC_LENGTH
        if barcode_start < 0:
            barcode_start = 0
        barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1
        if barcode_end > len(sequence):
            barcode_end = len(sequence)

        left_bc = sequence[barcode_start:linker_start]
        left_found = left_bc in self.left_barcodes
        right_bc = sequence[linker_end + 1:barcode_end]
        right_found = right_bc in self.right_barcodes

        for i in range(0, len(left_bc) - self.KMER_SIZE + 1, 1):
            self.left_kmer_dist[left_bc[i:i + self.KMER_SIZE]] += 1
        for i in range(0, len(right_bc) - self.KMER_SIZE + 1, 1):
            self.right_kmer_dist[right_bc[i:i + self.KMER_SIZE]] += 1

        potential_barcode = left_bc + right_bc
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=len(potential_barcode) - 1, score_diff=2)
        if not left_found: left_bc = None
        if not right_found: right_bc = None
        barcode_seq = potential_barcode if potential_barcode in self.barcodes else None
        return True, left_bc, right_bc, barcode_seq, barcode


class IlluminaReadStats:
    def __init__(self):
        self.read_count = 0
        self.linker_count = 0
        self.left_bc_count = 0
        self.right_bc_count = 0
        self.exact_match = 0
        self.whole_bc_count = 0

    def add_read(self, linker_found, left_bc, right_bc, bc_exact, bc_found):
        self.read_count += 1
        if linker_found:
            self.linker_count += 1
        if left_bc:
            self.left_bc_count += 1
        if right_bc:
            self.right_bc_count += 1
        if bc_exact:
            self.exact_match += 1
        if bc_found:
            self.whole_bc_count += 1

    def __str__(self):
        return "Total reads:\t%d\n" \
               "Linker found:\t%d\nLeft barcode detected:\t%d\nRight barcode detected:\t%d\n" \
               "Exact barcode matches:\t%d\nWhole barcode detected:\t%d\n" % \
            (self.read_count, self.linker_count, self.left_bc_count, self.right_bc_count, self.exact_match, self.whole_bc_count)


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_table = output_table
        self.read_stat = IlluminaReadStats()
        self.output_file = open(output_table, "w")

    def __del__(self):
        # logger.info("\n%s" % str(self.read_stat))
        stat_out = open(self.output_table + ".stats", "w")
        stat_out.write(str(self.read_stat))
        stat_out.close()
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
        logger.info("Finished " + input_file)

    def _process_fastx(self, read_handler):
        counter = 0
        for r in read_handler:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.id
            seq = str(r.seq)
            self._process_read(read_id, seq)

    def _process_bam(self, read_handler):
        counter = 0
        for r in read_handler:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.query_name
            seq = r.query_sequence
            self._process_read(read_id, seq)

    def _process_read(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        linker_found, left_bc_found, right_bc_found, exact_match, bc_found = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        self.read_stat.add_read(linker_found, left_bc_found, right_bc_found, exact_match, bc_found)
        self.output_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (read_id, str(linker_found), str(left_bc_found), str(right_bc_found), str(exact_match), str(bc_found)))

    def process_chunk(self, read_chunk):
        for read_id, seq in read_chunk:
            self._process_read(read_id, seq)


def print_dict(dict):
    kv_list = sorted(list(dict.items()), key=lambda x:x[1], reverse=True)
    print(kv_list[:20])


def process_single_thread(args):
    barcodes = load_barcodes(args.barcodes)
    barcode_detector = DoubleBarcodeDetector(barcodes)
    barcode_caller = BarcodeCaller(args.output, barcode_detector)
    barcode_caller.process(args.input)
    print_dict(barcode_detector.left_kmer_dist)
    print_dict(barcode_detector.right_kmer_dist)


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
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM)", required=True)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    set_logger(logger)
    process_single_thread(args)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
