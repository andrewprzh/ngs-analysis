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


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', metavar='BAM_FILE', type=str,  help='sorted and indexed BAM file')
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
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