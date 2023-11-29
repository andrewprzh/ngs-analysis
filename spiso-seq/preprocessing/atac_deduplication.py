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
from collections import defaultdict
import logging

import pysam

logger = logging.getLogger('ATACFilter')


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def overlaps_at_least(range1, range2, delta=0):
    cutoff = min([delta, range1[1] - range1[0] + 1, range2[1] - range2[0] + 1])
    overlap = min(range1[1], range2[1]) + 1 - max(range1[0], range2[0])
    return overlap >= cutoff


def load_barcodes(read2barcode_file, read_column=0, barcode_column=1, delim='\t'):
    read2barcode_dict = {}
    for l in open(read2barcode_file):
        if l.startswith("#"): continue
        v = l.strip().split(delim)
        barcode = v[barcode_column]
        if barcode == "*": continue
        read_id = v[read_column]
        if read_id.startswith("@"):
            read_id = read_id[1:]
        read2barcode_dict[read_id] = barcode

    logger.info("Loaded %d barcoded reads" % len(read2barcode_dict))
    return read2barcode_dict


class IntervalWithAlignments:
    def __init__(self, alignment):
        self.chr_id = alignment.reference_name
        self.start = alignment.reference_start
        self.end = alignment.reference_end
        self.alignment_list = [alignment]

    def next_alignment_overlaps(self, alignment):
        return (self.chr_id == alignment.reference_name and
                overlaps((self.start, self.end), (alignment.reference_start, alignment.reference_end)))

    def extend(self, alignment):
        self.end = max(self.end, alignment.reference_end)
        self.alignment_list.append(alignment)


def filter_atac(bam_file, read2barcode_dict, out_bam):
    bam = pysam.AlignmentFile(bam_file, "rb")
    outf = pysam.AlignmentFile(out_bam, "wb", template=bam)
    logger.info("Filtering %s, writing to %s" % (bam_file, out_bam))

    current_intervals = {}
    in_reads = 0
    barcoded_reads = 0
    out_reads = 0
    for a in bam:
        if a.reference_id == -1 or a.is_supplementary or a.is_secondary:
            continue
        in_reads += 1
        read_id = a.query_name
        if read_id not in read2barcode_dict:
            continue
        barcoded_reads += 1
        barcode = read2barcode_dict[read_id]
        if barcode not in current_intervals:
            current_intervals[barcode] = IntervalWithAlignments(a)
            continue

        current_interval = current_intervals[barcode]
        if current_interval.next_alignment_overlaps(a):
            current_interval.extend(a)
        else:
            outf.write(current_interval.alignment_list[0])
            out_reads += 1
            current_intervals[barcode] = IntervalWithAlignments(a)

    # dump all the remaining alignments
    for current_interval in current_intervals.values():
        out_reads += 1
        outf.write(current_interval.alignment_list[0])

    logger.info("Done. Processed %d primary alignments (%d barcoded), save %d alignments" % (in_reads, barcoded_reads, out_reads))
    outf.close()


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output BAM file", required=True)
    parser.add_argument("--input", "-i", type=str, help="input BAM file", required=True)
    parser.add_argument("--read2barcode", type=str, help="read to barcode table (TSV)", required=True)
    parser.add_argument("--read_column", type=int, default=0, help="read column in read2barcode table [0]")
    parser.add_argument("--barcode_column", type=int, default=1, help="barcode column in read2barcode table [1]")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    barcodes = load_barcodes(args.read2barcode, args.read_column, args.barcode_column)
    filter_atac(args.input, barcodes, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

