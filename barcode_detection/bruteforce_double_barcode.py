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

from common import *


logger = logging.getLogger('BarcodeCaller')


READ_CHUNK_SIZE = 10000


class ReadStatsShort:
    def __init__(self):
        self.read_count = 0
        self.linker_count = 0
        self.bc_count = 0

    def add_read(self, linker_found, bc):
        self.read_count += 1
        if linker_found:
            self.linker_count += 1
        if bc is not None:
            self.bc_count += 1

    def __str__(self):
        return "Total reads:\t%d\n" \
               "Linker found:\t%d\nBarcode detected:\t%d\n" % \
            (self.read_count, self.linker_count, self.bc_count)


class DoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    LEFT_BC_LENGTH = 2
    RIGHT_BC_LENGTH = 12
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH

    def __init__(self, joint_barcode_list):
        self.barcode_set = set(joint_barcode_list)

    def find_barcode_umi(self, read_id, sequence):
        linker_found, barcode = self._find_barcode_umi_fwd(sequence)
        if linker_found:
            return linker_found, barcode

        rev_seq = reverese_complement(sequence)
        return self._find_barcode_umi_fwd(rev_seq)

    def _find_barcode_umi_fwd(self, sequence):
        pos = sequence.find(self.LINKER)
        if pos == -1:
            return False, None

        bc_start = max(0, pos - self.LEFT_BC_LENGTH)
        barcode = sequence[bc_start:pos]
        linker_end = pos + len(self.LINKER)
        bc_end = min(len(sequence), linker_end + 1 + self.RIGHT_BC_LENGTH)
        barcode += sequence[linker_end + 1:bc_end]
        if len(barcode) != self.BC_LENGTH or barcode not in self.barcode_set:
            return True, None
        return True, barcode


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_table = output_table
        self.output_file = open(output_table, "w")
        self.read_stat = ReadStatsShort()

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
#            if counter % 100 == 0:
#                sys.stdout.write("Processed %d reads\r" % counter)
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
        linker_found, barcode = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        self.output_file.write("%s\t%s\t%s\n" % (read_id, linker_found, barcode))
        self.read_stat.add_read(linker_found, barcode)

    def process_chunk(self, read_chunk):
        for read_id, seq in read_chunk:
            self._process_read(read_id, seq)


def fastx_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        current_chunk.append((r.id, str(r.seq)))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def bam_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        current_chunk.append((r.query_name, r.query_sequence))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def process_chunk(barcodes, read_chunk, output_file, num):
    output_file += "_" + str(num)
    barcode_detector = DoubleBarcodeDetector(barcodes)
    barcode_caller = BarcodeCaller(output_file, barcode_detector)
    barcode_caller.process_chunk(read_chunk)
    return output_file


def process_single_thread(args):
    barcodes = load_barcodes(args.barcodes)
    barcode_detector = DoubleBarcodeDetector(barcodes)
    barcode_caller = BarcodeCaller(args.output, barcode_detector)
    barcode_caller.process(args.input)


def process_in_parallel(args):
    barcodes = load_barcodes(args.barcodes)
    logger.info("Loaded %d barcodes" % len(barcodes))

    input_file = args.input
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
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb"))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    while os.path.exists(tmp_dir):
        tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    os.makedirs(tmp_dir)

    barcode_calling_gen = (
        process_chunk,
        itertools.repeat(barcodes),
        read_chunk_gen,
        itertools.repeat(os.path.join(tmp_dir, "bc")),
        itertools.count(start=0, step=1),
    )

    with ProcessPoolExecutor(max_workers=args.threads) as proc:
        output_files = proc.map(*barcode_calling_gen, chunksize=1)
    outf = open(args.output, "w")
    stat_dict = defaultdict(int)
    for tmp_file in output_files:
        shutil.copyfileobj(open(tmp_file, "r"), outf)
        for l in open(tmp_file + ".stats", "r"):
            v = l.strip().split("\t")
            if len(v) != 2:
                continue
            stat_dict[v[0]] += int(v[1])

    logger.info(stat_dict)
    shutil.rmtree(tmp_dir)
    logger.info("Done")


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
    parser.add_argument("--threads", "-t", type=int, help="threads to use", default=8)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    set_logger(logger)
    if args.threads == 1:
        process_single_thread(args)
    else:
        process_in_parallel(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
