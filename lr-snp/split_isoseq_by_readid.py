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
                                     description="A simple script for splitting Iso-seq BAM files using read IDs.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file with CB tag", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str, default="./")
    required_group.add_argument("--list", help="List with barcodes", type=str)
    required_group.add_argument("--split_by", help="split by: BC (barcode) or GR (group)", type=str)
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args




def split_by_barcode(args, barcodes):
    inf = pysam.AlignmentFile(args.bam, "rb")
    splitted_reads = {}

    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        read_id = read.query_name
        cell_type = read_id.split(":")[0].split('_')[1]
        barcode = read_id.split(":")[1]
        if len(barcodes) > 0 and barcode not in barcodes:
            continue

        property = (cell_type, barcode)
        if property not in splitted_reads:
            splitted_reads[property] = []
        splitted_reads[property].append(read)
    print("Done                    ")
    inf.close()

    print("Dumping reads to files")
    for property in splitted_reads.keys():
        cell_type = property[0]
        barcode = property[1]
        out_dir = os.path.join(args.output_prefix, cell_type)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        bc_file_name = os.path.join(out_dir, barcode +  ".bam")
        bc_file = pysam.AlignmentFile(bc_file_name, "wb", template=inf)
        for read in splitted_reads[property]:
            bc_file.write(read)
        bc_file.close()

        pysam.sort("-o", 'tmp.bam', bc_file_name)
        os.rename('tmp.bam', bc_file_name)
        pysam.index(bc_file_name)


def split_by_group(args):
    inf = pysam.AlignmentFile(args.bam, "rb")
    splitted_reads = {}

    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        read_id = read.query_name
        cell_type = read_id.split(":")[0].split('_')[1]
        property = cell_type
        if property not in splitted_reads:
            splitted_reads[property] = []
        splitted_reads[property].append(read)
    print("Done                    ")
    inf.close()

    print("Dumping reads to files")
    for property in splitted_reads.keys():
        cell_type = property
        bc_file_name = os.path.join(args.output_prefix, cell_type +  ".bam")
        bc_file = pysam.AlignmentFile(bc_file_name, "wb", template=inf)
        for read in splitted_reads[property]:
            bc_file.write(read)
        bc_file.close()

        pysam.sort("-o", 'tmp.bam', bc_file_name)
        os.rename('tmp.bam', bc_file_name)
        pysam.index(bc_file_name)


def main():
    args = parse_args()
    barcodes = set()
    if args.list is not None:
        for l in open(args.list):
            barcodes.add(l.split('/')[1].split('.')[0])

    if args.split_by == "BC":
        split_by_barcode(args, barcodes)
    elif args.split_by == "GR":
        split_by_group(args)
    else:
        print("Invalid split by option")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
