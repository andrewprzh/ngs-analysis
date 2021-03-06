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
    required_group.add_argument("--split_by", help="split by: UMI (will use UB tag) or BC (CB tag)", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str, default="./")
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)

    if args.split_by != "UMI" and args.split_by != "BC":
        print("Invalid split_by option")
        parser.print_help()
        exit(-1)

    return args


def split_by_barcode(args, get_tag_function):
    inf = pysam.AlignmentFile(args.bam, "rb")
    splitted_reads = {}
    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        tag = get_tag_function(read)
        if tag is None:
            continue
        if tag not in splitted_reads:
            splitted_reads[tag] = []
        splitted_reads[tag].append(read)
    print("Done                    ")
    inf.close()

    print("Dumping reads to files")
    for tag in splitted_reads.keys():
        bc_file = pysam.AlignmentFile(args.output_prefix + tag + ".bam", "wb", template=inf)
        for read in splitted_reads[tag]:
            bc_file.write(read)
        bc_file.close()

    print("Sorting files")
    for bc in splitted_reads.keys():
        tag_file_name = args.output_prefix + bc + ".bam"
        pysam.sort("-o", 'tmp.bam', tag_file_name)
        os.rename('tmp.bam', tag_file_name)
        pysam.index(tag_file_name)


def get_UMI(read):
    if not read.has_tag("UB"):
        return None
    umi = read.get_tag("UB")
    return umi.split(":")[-1]


def get_barcode(read):
    if not read.has_tag("CB"):
        return None
    barcode = read.get_tag("CB")
    return barcode.split(":")[-1][:-2]


def main():
    args = parse_args()
    if args.split_by == "UMI":
        split_by_barcode(args, get_UMI)
    elif args.split_by == "BC":
        split_by_barcode(args, get_barcode)
    else:
        print("Unknown split by option")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
