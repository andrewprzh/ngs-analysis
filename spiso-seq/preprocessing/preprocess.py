############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from traceback import print_exc
import argparse

from Bio import SeqIO

from common import *

DATA_TYPES = {'assembly' : 'gmap', '10x_assembly' : 'gmap', 'isoseq' : 'minimap', 'ccs' : 'star', 'ont' : 'minimap'}
ALIGNERS = set(['minimap', 'gmap', 'star'])

def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for spliced mapping of long reads/contigs using gmap/star/minimap2. "
                                                 "Uses slightly tuned set of parameters for each data type. "
                                                 "Correctly processes contigs obtained from 10x RNA data unsing rnaSPAdes.\n"
                                                 "Available data types: assembly, 10x_assembly (gmap will be used), "
                                                 "isoseq, ont (minimap will be used) and ccs (STAR long will be used). ")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument('--input', '-i', dest='fasta_files', nargs='+', help='list of FASTA / FASTQ files')
    required_group.add_argument("--index", "-x", help="reference genome index; should be consistent with selected aligner", type=str)
    required_group.add_argument("--data_type", "-d", help="assembly/10x_assembly/isoseq/ccs/ont", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    optional_group = parser.add_argument_group('optional parameters')
    optional_group.add_argument("--aligner", "-a", help="force to use aligner auto/minimap/star/gmap "
                                                        "(selected automatically based on data type)", type=str)
    optional_group.add_argument("--delim", help="delimeter between read id and barcode for 10x assembly "
                                                      "(default:_barcodeIDs_) ", type=str, default='_barcodeIDs_')
    optional_group.add_argument("--delim2", help="delimeter between barcodes "
                                                      "(default:comma) ", type=str, default=',')
    optional_group.add_argument("--threads", "-t", help="threads for aligner (16)", type=int, default=16)
    optional_group.add_argument("--max_len", help="skip contigs longer than (0), valid for 10x assembly only", type=int, default=16)


    args = parser.parse_args()

    if args.fasta_files is None or len(args.fasta_files) == 0 or args.index is None or args.output_prefix is None:
        print("ERROR: not enough parameters")
        parser.print_help()
        exit(-1)

    if args.data_type is None or args.data_type not in DATA_TYPES:
        print("ERROR: specify data type among assembly/10x_assembly/isoseq/ccs/ont")

    if args.aligner is None:
        args.aligner = DATA_TYPES[args.data_type]
        print(args.aligner + "will be used, make sure index is provided accordingly")
    elif args.aligner not in ALIGNERS:
        print("ERROR: specify aligner among assembly/10x_assembly/isoseq/ccs/ont")

    return args


def main():
    args = parse_args()

    if args.data_type == "10x_assembly":
        for inf in args.fasta_files:
            contigs_name, short_id_contigs_name = convert_fasta_with_barcodes(inf, args)
            align_fasta(short_id_contigs_name, args, contigs_name)
    else:
        for inf in args.fasta_files:
            align_fasta(inf, args)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



