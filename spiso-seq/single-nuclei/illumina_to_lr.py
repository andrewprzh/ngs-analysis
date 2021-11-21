#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import subprocess
import sys
import argparse
import pysam
from Bio import SeqIO
from traceback import print_exc
from collections import defaultdict
import numpy


MAX_LONG_READS = 10000


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
    print("Reading all info from " + args.all_info)
    read_dict = defaultdict(list)
    max_col = max(args.illumina_column, args.lr_column) + 1
    for l in open(args.all_info):
        vals = l.strip().split()
        if len(vals) < max_col:
            continue
        sr_read_id = vals[args.illumina_column]
        lr_read_id = vals[args.lr_column]
        read_dict[lr_read_id].append(sr_read_id)
    return read_dict


def select_reads(args, read_dict):
    selected_ids = {}
    lr_count = 0
    print("Extracting long reads from " + args.long_reads)
    for r in SeqIO.parse(args.long_reads, 'fastq'):
        if r.id in read_dict:
            lr_count += 1
            fprefix = r.id.replace('/', '_')
            fname = os.path.join(args.output, fprefix + '.fasta')
            if os.path.exists(fname):
                continue
            SeqIO.write([r], fname, 'fasta')
            lr_count += 1
            for sr_id in read_dict[r.id]:
                selected_ids[sr_id] = fprefix
        if lr_count > MAX_LONG_READS:
            break

    print("Extracting short reads from " + args.illumina)
    #print(selected_ids)
    selected_records = defaultdict(list)
    for r in SeqIO.parse(args.illumina, 'fastq'):
        if r.id in selected_ids:
            selected_records[selected_ids[r.id]].append(r)

    print("Writing short reads")
    total_reads = 0
    total_files = 0
    for fprefix in selected_records.keys():
        fname = os.path.join(args.output, fprefix + '.sr.fastq')
        total_files += 1
        total_reads += len(selected_records[fprefix])
        if os.path.exists(fname):
            continue
        SeqIO.write(selected_records[fprefix], fname, 'fastq')

    print("Saved %d reads into %d files " % (total_reads, total_files))
    return selected_records.keys()


def map_reads(args, file_prefixes):
    print("Running minimap2 for %d sets" % len(file_prefixes))
    for fprefix in file_prefixes:
        sr = os.path.join(args.output, fprefix + '.sr.fastq')
        lr = os.path.join(args.output, fprefix + '.fasta')
        bam = os.path.join(args.output, fprefix + '.bam')
        if os.path.exists(bam):
            continue
        os.system("minimap2 %s %s -t 4 -ax sr 2> /dev/null | samtools sort -o %s" % (lr, sr, bam))
        os.system("samtools index " + bam)


def count_tlen(args, file_prefixes):
    tlens = []
    print("Running position analysis for %d sets" % len(file_prefixes))
    for fprefix in file_prefixes:
        bam = os.path.join(args.output, fprefix + '.bam')
        bam_file = pysam.AlignmentFile(bam, 'r')
        long_read_len = bam_file.lengths[0]
        for a in bam_file:
            if a.flag & 16 == 0:
                tlens.append(long_read_len - a.reference_start)
            else:
                tlens.append(a.reference_start + a.query_length)

    with open(os.path.join(args.output, 'stats.tsv'), 'w') as outf:
        outf.write("%s\t%s\n\n" % (str(numpy.mean(tlens)), str(numpy.std(tlens))))
        tlen_hist, tlen_bins = numpy.histogram(tlens, bins=[i for i in range(max(tlens) + 1)])
        outf.write('\t'.join([str(x) for x in tlen_hist]) + '\n')
        outf.write('\t'.join([str(x) for x in tlen_bins]) + '\n')


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    read_dict = read_all_info(args)
    file_prefixes = select_reads(args, read_dict)
    map_reads(args, file_prefixes)
    count_tlen(args, file_prefixes)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
