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
from Bio import Seq
from Bio import SeqIO
from traceback import print_exc
from collections import defaultdict
import numpy


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", required=True)
    parser.add_argument("--all_info", "-a", type=str, help="read IDs tab file", required=True)
    parser.add_argument("--barcode_info", "-b", type=str, help="read barcode info", required=True)

    parser.add_argument("--long_reads", "-l", type=str, help="input file with ONT/PB sequences (FASTQ/BAM)", required=True)
    parser.add_argument("--illumina", "-i", type=str, help="input file with Illumina sequences (FASTQ/BAM)", required=True)
    parser.add_argument("--illumina_column", type=int, help="illumina read id column in all info", default=1)
    parser.add_argument("--lr_column", type=int, help="long read id column in all info", default=0)
    parser.add_argument("--max_long_reads", type=int, help="max long reads to extract", default=10000)
    parser.add_argument("--minimap2", type=str, help="path to minimap2 (expected to be in $PATH if not set", default='')

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
        if sr_read_id[0] == '@':
            sr_read_id = sr_read_id[1:]
        lr_read_id = vals[args.lr_column]
        if lr_read_id[0] == '@':
            lr_read_id = lr_read_id[1:]
        read_dict[lr_read_id].append(sr_read_id)
    return read_dict


def read_barcode_positions(info_file):
    print("Reading barcodes from " + info_file)
    read_dict = {}
    for l in open(info_file):
        vals = l.strip().split()
        if len(vals) < 4:
            continue
        read_id = vals[0]
        if read_id[0] == '@':
            read_id = read_id[1:]

        read_dict[read_id] = int(vals[3])
    print('Loaded %d barcode offsets ' % len(read_dict))
    return read_dict


def select_reads(args, read_dict):
    selected_ids = {}
    lr_count = 0
    print("Extracting long reads from " + args.long_reads)
    if args.long_reads.endswith('bam'):
        for r in pysam.AlignmentFile(args.long_reads, 'rb'):
            r_id = r.query_name
            if r_id in read_dict:
                lr_count += 1
                fprefix = r_id.replace('/', '_')
                fname = os.path.join(args.output, fprefix + '.fasta')
                if not os.path.exists(fname):
                    seq_record = SeqIO.SeqRecord(seq=Seq.Seq(r.query_sequence), id=r_id,  description="", name="")
                    SeqIO.write([seq_record], fname, 'fasta')
                for sr_id in read_dict[r_id]:
                    selected_ids[sr_id] = fprefix
            if lr_count > args.max_long_reads:
                break
    else:
        for r in SeqIO.parse(args.long_reads, 'fastq'):
            if r.id in read_dict:
                lr_count += 1
                fprefix = r.id.replace('/', '_')
                fname = os.path.join(args.output, fprefix + '.fasta')
                if not os.path.exists(fname):
                    SeqIO.write([r], fname, 'fasta')
                for sr_id in read_dict[r.id]:
                    selected_ids[sr_id] = fprefix
            if lr_count > args.max_long_reads:
                break

    print("Extracting short reads from " + args.illumina)
    #print(selected_ids)
    selected_records = defaultdict(list)

    if args.illumina.endswith('bam'):
        for r in pysam.AlignmentFile(args.illumina, 'rb'):
            r_id = r.query_name
            if r_id in selected_ids:
                seq_record = SeqIO.SeqRecord(seq=Seq.Seq(r.query_sequence), id=r_id, description="", name="",
                                             letter_annotations={'phred_quality': r.query_qualities})
                selected_records[selected_ids[r_id]].append(seq_record)
    else:
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
        if not os.path.exists(fname):
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
        minimap2 = os.path.join(args.minimap2, 'minimap2') if args.minimap2 else 'minimap2'
        os.system("%s %s %s -t 4 -ax sr 2> /dev/null | samtools sort -o %s" % (minimap2, lr, sr, bam))
        os.system("samtools index " + bam)


def count_tlen(args, file_prefixes, barcode_offsets=None):
    tlens = []
    print("Running position analysis for %d sets" % len(file_prefixes))
    for fprefix in file_prefixes:
        bam = os.path.join(args.output, fprefix + '.bam')
        bam_file = pysam.AlignmentFile(bam, 'r')
        ont_read_id = bam_file.get_reference_name(0)
        long_read_len = bam_file.lengths[0]
        for a in bam_file:
            if barcode_offsets is None or ont_read_id not in barcode_offsets:
                offset = 0
            else:
                offset = barcode_offsets[ont_read_id]
            if a.flag & 16 == 0:
                tlens.append(long_read_len - a.reference_start - offset)
            else:
                tlens.append(a.reference_start + a.query_length - offset)
    tlens = list(filter(lambda x: x < 1000, tlens))

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
    barcode_offsets = read_barcode_positions(args.barcode_info)
    file_prefixes = select_reads(args, read_dict)
    map_reads(args, file_prefixes)
    count_tlen(args, file_prefixes, barcode_offsets)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
