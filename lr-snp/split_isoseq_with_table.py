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
                                     description="A simple script for splitting Iso-seq BAM files using read ID table.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file with CB tag", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str, default="./")
    required_group.add_argument("--table", help="table with read ids, barcodes and UMIs", type=str)
    required_group.add_argument("--split_by", help="split by: UMI, BC (barcode) or GR (group)", type=str)

    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)

    if args.table is None:
        parser.print_help()
        exit(-1)

    return args


def get_read_table(table_file):
    read_table = {}
    for l in open(table_file):
        tokens = l.strip().split()
        if len(tokens) != 4:
            continue
        read_id = tokens[0]
        if read_id in read_table:
            print("Duplicated read id")
            continue

        read_table[read_id] = tokens[1:3]
    return read_table


def get_split_by_column_index(args):
    if args.split_by == "UMI":
        return 2
    elif args.split_by == "BC":
        return 0
    elif args.split_by == "GR":
        return 1
    else:
        print("Unsupported split_by option")
        return None


def split_by_barcode(args, read_table):
    inf = pysam.AlignmentFile(args.bam, "rb")
    splitted_reads = {}
    split_by_column_index = get_split_by_column_index(args)
    if split_by_column_index is None:
        return

    for read in inf:
        read_id = read.query_name
        if read_id not in read_table:
            continue

        property = read_table[read_id][split_by_column_index]
        if property not in splitted_reads:
            splitted_reads[property] = []
        splitted_reads[property].append(read)
    inf.close()

    for property in splitted_reads.keys():
        bc_file = pysam.AlignmentFile(args.output_prefix + property + ".bam", "wb", template=inf)
        for read in splitted_reads[property]:
            bc_file.write(read)
        bc_file.close()

    for bc in splitted_reads.keys():
        splitted_file_name = args.output_prefix + bc + ".bam"
        pysam.sort("-o", 'tmp.bam', splitted_file_name)
        os.rename('tmp.bam', splitted_file_name)
        pysam.index(splitted_file_name)


def main():
    args = parse_args()
    read_table = get_read_table(args.table)
    split_by_barcode(args, read_table)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
