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
from Bio import SeqIO
from traceback import print_exc


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bam', nargs='+', type=str,  help='sorted and indexed BAM file(s), '
                                                            'each file will be treated as a separate sample')
    parser.add_argument('--fastq', nargs='+', type=str, help='input FASTQ file(s), '
                                                             'each file will be treated as a separate sample'
                                                             'reference genome should be provided when using raw reads')
    parser.add_argument('--sample_list', type=str, help='test file with input files, each line is treated as '
                                                        'a separate sample, but can contain multiple files')
    parser.add_argument('--labels', nargs='+', type=str, help='sample names to be used')
    parser.add_argument("--data_type", "-d", help="type of data to process, supported types are: "
                                                  "assembly, raw_long_reads, hq_long_reads, barcoded_short_reads", type=str)

    parser.add_argument("--genedb", "-g", help="gene database in gffutils .db format", type=str)
    parser.add_argument("--gtf", help="gene database in GTF/GFF format", type=str)
    parser.add_argument("--reference", help="reference genome in FASTA format,"
                                            "should be provided only when raw reads are used as an input", type=str)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically", type=str)
    parser.add_argument("--prefix", help="prefix for output files", type=str, default="")
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="1")
    parser.add_argument("--read_info", help="text file with tab-separated information about input reads, according to"
                                            "which counts are groupped, e.g. cell type, barcode, etc.", type=str)
    parser.add_argument("--aligner", help="force to use this alignment method, can be minimap2, star, gmap, hisat2", type=str)


    args = parser.parse_args()
    return args

#Tune algorithm params
def set_params(args):
    pass


def main():
    args = parse_args()
    set_params(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)