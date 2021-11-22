#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from Bio import SeqIO
from Bio import Seq
import pysam
from traceback import print_exc
import numpy


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output FASTQ file", default="ont_truncated.fastq")
    parser.add_argument("--info", type=str, help="TSV barcode detection", required=True)
    parser.add_argument("--ont", "-i", type=str, help="input file with ONT sequences (BAM)", required=True)
    parser.add_argument("--read_length", type=int, help="illumina read length [76]", default=76)
    parser.add_argument("--mean_frag", type=int, help="mean fragment length [340]", default=340)
    parser.add_argument("--stdev_frag", type=int, help="standard deviation of the fragment length [90]", default=90)
    parser.add_argument("--seed", type=int, help="seed [11]", default=11)

    args = parser.parse_args()
    return args


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def read_all_info(info_file):
    print("Reading read info from " + info_file)
    read_dict = {}
    for l in open(info_file):
        vals = l.strip().split()
        if len(vals) < 3:
            continue
        read_id = vals[0]
        if read_id[0] == '@':
            read_id = read_id[1:]

        read_dict[read_id] = 1 if vals[2] == 'rev' else 0
    return read_dict


def truncate_reads(args, read_dict):
    print("Extracting sequences from " + args.ont)
    total_sequences = 0
    with open(args.output, 'w') as outf:
        seq_records = []
        for r in pysam.AlignmentFile(args.ont, 'rb'):
            if r.id not in read_dict:
                continue
            assigned_strand = read_dict[r.id]
            mapped_strand = (r.flag & 16) >> 4
            seq = r.query_sequence
            quality = ''.join(map(lambda x: chr(x+33), r.query_qualities)) if r.query_qualities else ''
            tlen = numpy.random.normal(args.mean_frag, args.stdev_frag)
            if assigned_strand ^ mapped_strand:
                end = min(len(seq), tlen)
                start = end - args.read_length
                seq = seq[start:end]
                if r.query_qualities:
                    qual = quality[start:end]
            else:
                start = max(0, len(seq) - tlen)
                end = start + args.read_length
                seq = seq[start:end]
                if r.query_qualities:
                    qual = quality[start:end]
            total_sequences += 1
            seq_records.append(SeqIO.SeqRecord(seq=Seq.Seq(seq), id=r.id + "_%d" % tlen, description="", name="",
                                               letter_annotations=qual))

            if len(seq_records) > 10000:
                SeqIO.write(seq_records, outf, 'fastq')
                seq_records = []


def main():
    args = parse_args()
    random.seed(args.seed)
    read_dict = read_all_info(args.info)
    truncate_reads(args, read_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
