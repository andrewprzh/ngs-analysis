############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
import pysam
from traceback import print_exc

def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for splitting 10x BAM files by UMI UB tag.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file with UB tag", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str, default="./")
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def split_by_barcode(args):
    inf = pysam.AlignmentFile(args.bam, "rb")
    barcode_reads = {}
    for read in inf:
        if not read.has_tag("UB"):
            continue
        barcode = read.get_tag("UB")
        barcode = barcode.split(":")[-1]
        if barcode not in barcode_reads:
            barcode_reads[barcode] = []
        barcode_reads[barcode].append(read)
    inf.close()

    for barcode in barcode_reads.keys():
        bc_file = pysam.AlignmentFile(args.output_prefix + barcode + ".bam", "wb", template=inf)
        for read in barcode_reads[barcode]:
            bc_file.write(read)
        bc_file.close()

    for bc in barcode_reads.keys():
        barcode_file_name = args.output_prefix + bc + ".bam"
        pysam.sort("-o", 'tmp.bam', barcode_file_name)
        os.rename('tmp.bam', barcode_file_name)
        pysam.index(barcode_file_name)


def main():
    args = parse_args()
    split_by_barcode(args)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
