############################################################################
# Copyright (c) 2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
import pysam
import gffutils
from Bio import Seq
from traceback import print_exc
import random


nucls = ['A', 'C', 'G', 'T']


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for adding cellranger tags into BAM for simulated data.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output BAM", type=str, required=True)
    required_group.add_argument("--tsv", help="TSV with barcodes and UMIs", type=str, required=True)
    required_group.add_argument("--read_id_column", help="read id column in TSV (0-based)", type=int, default=0)
    required_group.add_argument("--barcode_column", help="barcode column in TSV (0-based)", type=int, default=4)
    required_group.add_argument("--umi_column", help="UMI column in TSV (0-based)", type=int, default=7)

    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def load_barcodes(args):
    read_bc_umi = {}
    print("Loading " + args.bam)
    bam_reads = set()
    for a in pysam.AlignmentFile(sys.argv[1], "r"):
        bam_reads.add(a.query_name)

    min_cols = max([args.read_id_column, args.barcode_column, args.umi_column]) + 1
    for l in open(args.tsv):
        v = l.strip().split('\t')
        if len(v) < min_cols:
            continue
        read_id = v[args.read_id_column]
        if read_id[0] == '@':
            read_id = read_id[1:]
        if read_id in bam_reads:
            read_bc_umi[read_id] = (v[args.barcode_column], v[args.umi_column])

    return read_bc_umi


def add_barcode_tags(args, read_bc_umi):
    inf = pysam.AlignmentFile(args.bam, "rb")
    outf = pysam.AlignmentFile(args.output, "wb", template=inf)
    print("Reading " + args.bam)

    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        read_id = read.query_name
        if read_id not in bc_umi_pair:
            print("Missing read" + read_id)
            continue
        bc_umi_pair = read_bc_umi[read_id]
        barcode = bc_umi_pair[0]
        umi = bc_umi_pair[1]

        # GX:Z:ENSG00000243485    GN:Z:MIR1302-2HG        fx:Z:ENSG00000243485
        # RE:A:E  xf:i:0  CR:Z:AAACGCTCACGCAGTC   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAACGCTCACGCAGTC-1 UR:Z:CACATTGTCTCT       UY:Z:FFFFFFFFFFFF
        read.set_tag("CR", barcode, value_type='Z')
        read.set_tag("CY", 'F' * len(barcode), value_type='Z')
        read.set_tag("CB", barcode + '-1', value_type='Z')
        read.set_tag("UB", umi, value_type='Z')
        read.set_tag("UR", umi, value_type='Z')
        read.set_tag("UY", 'F' * len(umi), value_type='Z')
        outf.write(read)

    outf.close()
    print("Saved to " + args.output)


def main():
    args = parse_args()
    read_bc_umi = load_barcodes(args)
    add_barcode_tags(args, read_bc_umi)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
