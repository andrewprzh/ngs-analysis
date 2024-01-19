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
import pysam
from Bio import SeqIO


def process_fastx(read_handler):
    counter = 0
    lengths = []
    for r in read_handler:
        if counter % 100 == 0:
            sys.stdout.write("Processed %d reads\r" % counter)
        counter += 1
        lengths.append(len(str(r.seq)))

    return sorted(lengths)


def process_bam(read_handler):
    counter = 0
    lengths = []
    for r in read_handler:
        if counter % 100 == 0:
            sys.stdout.write("Processed %d reads\r" % counter)
        counter += 1
        lengths.append(len(r.query_sequence))

    return sorted(lengths)

def get_lengths(input_file):
    print("Processing " + input_file)
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    if low_ext in ['.fq', '.fastq']:
        lengths = process_fastx(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        lengths = process_fastx(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        lengths = process_bam(pysam.AlignmentFile(input_file, "rb"))
    else:
        print("Unknown file format " + input_file)
        return
    print("Finished " + input_file)
    return fname, lengths


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output name [output + .read_length.csv]")
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--delim", type=str, help="separator to use (T for tab, N for new line)", default=",")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if args.delim.lower() == "t":
        args.delim = "\t"
    elif args.delim.lower() == "n":
        args.delim = "\n"

    fname, lengths = get_lengths(args.input)
    if args.output:
        out_fname = args.output
    else:
        if args.delim == ",":
            out_fname = fname + ".read_lengths.csv"
        elif args.delim == "\t" or args.delim == "\n":
            out_fname = fname + ".read_lengths.tsv"
        else:
            out_fname = fname + ".read_lengths.txt"

    outf =  open(out_fname, "w")
    outf.write(fname + args.delim + args.delim.join(map(str, lengths)))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
