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
from traceback import print_exc
from Bio import SeqIO, SeqRecord, Seq

base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverse_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


def split_sequence(read_id, seq, qual, read_length, insert_size, min_len):
    start_pos = 0
    while start_pos < len(seq):
        end_pos = min(len(seq), start_pos + insert_size)
        if end_pos - start_pos < min_len:
            break
        r1 = seq[start_pos:start_pos+read_length]
        q1 = qual[start_pos:start_pos+read_length]
        r2 = reverse_complement(seq[end_pos-read_length:end_pos])
        q2 = qual[end_pos-read_length:end_pos]
        yield read_id + "_" + str(start_pos), r1, q1, r2, q2
        start_pos += insert_size


def process(input_file, output_prefix, read_length, insert_size, min_len):
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()
    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")

    left_reads = gzip.open(output_prefix + "_R1_001.fastq.gz", "wt")
    right_reads = gzip.open(output_prefix + "_R2_001.fastq.gz", "wt")
    for r in SeqIO.parse(handle, "fastq"):
        read_id = r.id
        seq = r.seq
        qual = r.letter_annotations['phred_quality']

        for read_pair in split_sequence(read_id, seq, qual, read_length, insert_size, min_len):
            lr = SeqRecord.SeqRecord(seq=Seq.Seq(read_pair[1]), id=read_pair[0] + "/1", description="",
                                     letter_annotations={'phred_quality': read_pair[2]})
            rr = SeqRecord.SeqRecord(seq=Seq.Seq(read_pair[3]), id=read_pair[0] + "/2", description="",
                                     letter_annotations={'phred_quality': read_pair[4]})
            SeqIO.write(lr, left_reads, "fastq")
            SeqIO.write(rr, right_reads, "fastq")

    left_reads.close()
    right_reads.close()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ)", required=True)
    parser.add_argument("--read_length", "-r", type=int, help="read length of output reads", default=150)
    parser.add_argument("--insert_size", type=int, help="insert size for output reads", default=300)
    parser.add_argument("--min_ont_len", type=int, help="read length of output reads", default=180)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    process(args.input, args.output, args.read_length, args.insert_size, args.min_ont_len)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
