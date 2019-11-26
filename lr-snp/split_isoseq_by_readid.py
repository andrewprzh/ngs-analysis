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
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def split_by_barcode(args, barcodes):
    inf = pysam.AlignmentFile(args.bam, "rb")
    barcode_files = {}
    for read in inf:
        read_id = read.query_name
        cell_type = read_id.split(":")[0].split('_')[1]
#        if not cell_type.startswith("Gra"):
#            continue
        out_dir = args.output_prefix + cell_type
        barcode = read_id.split(":")[1]
        if len(barcodes) > 0 and barcode not in barcodes:
            continue
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        if (cell_type, barcode) not in barcode_files:
            barcode_files[(cell_type, barcode)] = pysam.AlignmentFile(os.path.join(out_dir, barcode +  ".bam") , "wb", template=inf)
        outf = barcode_files[(cell_type, barcode)]
        outf.write(read)

    inf.close()
    for bc_info in barcode_files.keys():
        barcode_files[bc_info].close()
        barcode_file_name =  os.path.join(args.output_prefix + bc_info[0], bc_info[1] + ".bam")
        pysam.sort("-o", 'tmp.bam', barcode_file_name)
        os.rename('tmp.bam', barcode_file_name)
        pysam.index(barcode_file_name)


def main():
    args = parse_args()
    barcodes = set()
    if args.list is not None:
        for l in open(args.list):
            barcodes.add(l.split('/')[1].split('.')[0])

    split_by_barcode(args, barcodes)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
