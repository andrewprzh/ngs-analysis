#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output dir name", required=True)
    parser.add_argument("--prefix", "-p", type=str, help="output prefix name", required=True)
    parser.add_argument("--fasta", "-f", type=str, help="input transcriptome", required=True)
    parser.add_argument("--bam_list", "-b", type=str, help="list of BAM files mapped to transcriptome", required=True)
    parser.add_argument("--no_secondary", help="ignore secondary alignments", action='store_true',  default=False)
    parser.add_argument("--cutoff", help="total counts cutoff", default=10)
    args = parser.parse_args()
    return args


def count_bams(bam_files, no_secondary=False):
    count_dict = {}
    for i, bam in enumerate(bam_files):
        for a in pysam.AlignmentFile(bam, "r"):
            if a.is_unmapped:
                continue
            if no_secondary and (a.is_secondary or a.is_supplementary):
                continue

            if a.reference_name not in count_dict:
                count_dict[a.reference_name] = [0 for j in range(len(bam_files))]
            count_dict[a.reference_name][i] += 1

    return count_dict


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    seq_dict = SeqIO.index(args.fasta, "fasta")
    bam_files = []
    for l in open(args.bam_list):
        bam_files.append(l.strip())

    count_dict = count_bams(bam_files, args.no_secondary)
    seq_records = []
    filtered_seq_records = []
    out_tsv = open(os.path.join(args.output, args.prefix + ".counts.tsv"), "w")
    out_tsv.write("#Gene\n" + "\t".join(bam_files) + "\n")
    for seq_id in sorted(count_dict.keys()):
        out_tsv.write(seq_id + "\t" + "\t".join([str(x) for x in count_dict[seq_id]]) + "\n")
        seq = seq_dict[seq_id]
        seq_records.append(seq)
        if sum(count_dict[seq_id]) >= args.cutoff:
            filtered_seq_records.append(seq)
    out_tsv.close()

    SeqIO.write(seq_records, os.path.join(args.output, args.prefix + ".expressed.fasta"), "fasta")
    SeqIO.write(filtered_seq_records, os.path.join(args.output, args.prefix + ".expressed10.fasta"), "fasta")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



