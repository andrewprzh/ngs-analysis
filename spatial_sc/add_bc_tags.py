############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
import pysam
import gffutils
from traceback import print_exc
import random
from collections import defaultdict


nucls = ['A', 'C', 'G', 'T']


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for adding cellranger tags into BAM.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output BAM", type=str, required=True)
    required_group.add_argument("--barcodes", "-b", help="barcodes in TSV format", type=str, required=True)
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


# 82718f26-b81f-4561-ac75-59b987941431    ACGCGTTTAAGACG  ATCGGCTAG       14      True    -       82      40      49      66
# LH00376:28:22GVKLLT3:4:1101:2655:1016   41      -1      8       25      +       CAACTGGCCGGGTA  TACCGTCGT       13      False

def get_read_dict(inf, trusted_only=False):
    read_dict = {}
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split('\t')
        if v[-1] in ['True', 'False']:
            # old format
            if v[6] == '*' or v[7] == '*': continue
            if trusted_only and v[-1] != "True": continue
            read_id = v[0]
            read_dict[read_id]= (v[6], v[7])
            continue

        if v[1] == '*' or v[2] == '*': continue
        if trusted_only and v[4] != "True": continue
        read_id = v[0]
        read_dict[read_id] = (v[1], v[2])

    return read_dict


def add_barcodes_tags(args, read_dict):
    inf = pysam.AlignmentFile(args.bam, "rb")
    outf = pysam.AlignmentFile(args.output, "wb", template=inf)
    print("Reading " + args.bam)

    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        read_id = read.query_name
        if read_id not in read_dict:
            continue
        barcode = read_dict[read_id][0]
        umi = read_dict[read_id][1]
        spliced = any(x[0] == 3 for x in read.cigartuples)

        # GX:Z:ENSG00000243485    GN:Z:MIR1302-2HG        fx:Z:ENSG00000243485
        # RE:A:E  xf:i:0  CR:Z:AAACGCTCACGCAGTC   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAACGCTCACGCAGTC-1 UR:Z:CACATTGTCTCT       UY:Z:FFFFFFFFFFFF
        #read.set_tag("GX", gene_id, value_type='Z')
        #read.set_tag("GN", gene_name_dict[gene_id], value_type='Z')
        #read.set_tag("fx", gene_id, value_type='Z')
        read.set_tag("CR", barcode, value_type='Z')
        read.set_tag("CY", 'F' * len(barcode), value_type='Z')
        read.set_tag("CB", barcode + '-1', value_type='Z')
        read.set_tag("UB", umi, value_type='Z')
        read.set_tag("UR", umi, value_type='Z')
        read.set_tag("UY", 'F' * len(umi), value_type='Z')
        read.set_tag("SP", "spliced" if spliced else "unspliced", value_type='Z')
        outf.write(read)
    outf.close()
    print("Saved to " + args.output)

def main():
    args = parse_args()
    read_dict = get_read_dict(args.barcodes)
    add_barcodes_tags(args, read_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
