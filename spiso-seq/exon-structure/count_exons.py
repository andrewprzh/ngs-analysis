############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import argparse
import pysam
from assign_isoforms_to_barcodes import *
from traceback import print_exc


def read_property_table(property_file, column):
    for l in open(property_file):


class 
def get_property_from_read_id(read_id)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', metavar='BAM_FILE', nargs='+', type=str, help='sorted and indexed BAM file(s)')
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    parser.add_argument("--readmap", "-r", help="read property table (first column - read id, other - property)", type=str)
    parser.add_argument("--property_column", help="read property column", type=int, default=1)
    parser.add_argument("--keep_terminal", help="do not suppress terminal", action='store_true', default=False)

    args = parser.parse_args()
    return args

#Tune algorithm params
def set_codon_count_params(args):
    args.reads_cutoff = 0
    args.data_type = "long_reads"
    args.change_chr_prefix = False
    args.assign_codons_when_ambiguous = False
    args.consider_flanking_junctions = False
    args.junction_delta = 1
    args.count_isoform_stats = False
    args.exon_count_mode = True


def main():
    args = parse_args()
    set_codon_count_params(args)

    db_processor = GeneDBProcessor(args)
    db_processor.process_all_genes()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)