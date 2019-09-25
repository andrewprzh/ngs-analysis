############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
#import gffutils
import argparse
import pysam
from traceback import print_exc
from Bio import SeqIO
from lr_snp_caller import *


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple tool for calling non-germline SNPs from a set of BAM files. "
                                                 "BAM files need to be sorted and indexed using samtools. "
                                                 "Requires a reference genome in FASTA format. "
                                                 "To call SNPs only within gene regions provide annotation in GTF / GFF or "
                                                 "gff utils database format.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument('--bam', '-b', dest='bam_file', nargs='+', help='list of sorted and indexed BAM file')
    required_group.add_argument("--reference", "-r", help="reference genome used for alignment in FASTA format", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    optional_group = parser.add_argument_group('optional parameters')
    optional_group.add_argument("--genedb", help="gene database in gffutils db format", type=str)
    optional_group.add_argument("--gtf", help="gene database in GTF / GFF", type=str)
    optional_group.add_argument("--keep_germline", help="do not filter out germline SNPs [False]", type=bool, default=False)
    optional_group.add_argument("--min_freq", "-f", help="minimal SNP frequency within a sample, between 0.0 and 1.0 [0.2]", type=float, default=0.2)
    optional_group.add_argument("--min_freq_factor", "-m", help="minimal SNP frequency factor, > 1.0 [2.0]", type=float, default=2.0)
    optional_group.add_argument("--min_cov", "-c", help="minimal SNP read coverage depth within a sample, > 0 [50]", type=int, default=50)
    optional_group.add_argument("--out_format", help="output format: VCF or TSV [TSV]", type=str, default='TSV')

    args = parser.parse_args()

    if args.bam_file is None or len(args.bam_file) == 0 or args.reference is None or args.output_prefix is None:
        parser.print_help()
        exit(-1)
    return args

#Tune algorithm params
def set_params(args):
    if len(args.bam_file) == 1:
        print("WARN: Only a single BAM files is provided, will not filter out germline SNPs.")
        args.keep_germline = True
    if args.min_freq < 0 or args.min_freq > 1:
        raise Exception("ERROR: minimal SNP frequency should be between 0.0 and 1.0, but set to " + str(args.min_freq))
    if args.min_freq_factor < 1:
        raise Exception("ERROR: minimal SNP frequency factor should be larger than 1.0, but set to " + str(args.min_freq))
    if args.min_cov < 1:
        raise Exception("ERROR: minimal SNP coverage should be positive, but set to " + str(args.min_cov))
    if args.genedb or args.gtf:
        raise Exception("ERROR: gene SNPs are not supported yet")
    if args.out_format != 'TSV':
        raise Exception("ERROR: only TSV format is currently supported")
    args.window_lenth = 1000000

    for bam in args.bam_file:
        if not os.path.isfile(bam):
            raise Exception("BAM file " + bam + " does not exist")
        bamfile = pysam.AlignmentFile(bam, "rb")
        if not bamfile.has_index:
            raise Exception("BAM file " + bam + " is not indexed, run samtools index")


def main():
    args = parse_args()
    set_params(args)

    snp_caller = SNPCaller(args)
    snp_map = snp_caller.process()

    sample_names = map(lambda x: os.path.splitext(os.path.basename(x))[0], args.bam_file)
    #print_snp_map(snp_map, sample_names)

    snp_writer = SNPMapTSVWriter(args.output_prefix, sample_names, args)
    snp_writer.dump_to_file(snp_map)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)