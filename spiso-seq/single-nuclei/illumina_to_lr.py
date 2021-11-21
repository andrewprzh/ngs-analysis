#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gffutils
from Bio import SeqIO
import gzip
from traceback import print_exc
from collections import defaultdict


MAX_LONG_READS = 10


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", required=True)
    parser.add_argument("--all_info", "-a", type=str, help="read IDs", required=True)
    parser.add_argument("--long_reads", "-l", type=str, help="input file with ONT/PB sequences", required=True)
    parser.add_argument("--illumina", "-i", type=str, help="input file with Illumina sequences", required=True)
    parser.add_argument("--illumina_column", type=int, help="illumina read id column in all info", default=2)
    parser.add_argument("--lr_column", type=int, help="long read id column in all info", default=0)

    args = parser.parse_args()
    return args


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def read_all_info(args):
    read_dict = {}
    max_col = max(args.illumina_column, args.lr_column) + 1
    for l in open(args.all_info):
        vals = l.strip().split()
        if len(vals) < max_col:
            continue
        sr_read_id = vals[args.illumina_column]
        lr_read_id = vals[args.lr_column]
        read_dict[lr_read_id] = sr_read_id
    return read_dict


def select_reads(args, read_dict):
    selected_ids = {}

    print("Extracting long reads from " + args.long_reads)
    for r in SeqIO.parse(args.long_reads, 'fastq'):
        if r.id in read_dict:
            fprefix = r.id.replace('/', '_')
            fname = os.path.join(args.output, fprefix + '.fasta')
            SeqIO.write([r], fname, 'fasta')
            for sr_id in read_dict[r.id]:
                selected_ids[sr_id] = fprefix
        if len(selected_ids) > MAX_LONG_READS:
            break

    print("Extracting long reads from " + args.long_reads)
    selected_records = defaultdict(list)
    for r in SeqIO.parse(args.illumina, 'fastq'):
        if r.id in selected_ids:
            selected_records[selected_ids[r.id]].append(r)

    print("Writing short reads")
    for fprefix in selected_records.keys():
        fname = os.path.join(args.output, fprefix + 'sr.fastq')
        SeqIO.write(selected_records[fname], fname, 'fastq')

    print("Done")


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
