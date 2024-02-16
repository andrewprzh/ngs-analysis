############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import argparse
import pysam
from traceback import print_exc
import random


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for subsampling BAM files.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output BAM", type=str, required=True)
    required_group.add_argument("--fraction", "-f", help="fraction of reads to keep (0,1)", type=float, default=0.1)
    required_group.add_argument("--seed", "-s", help="annotation in gffutils DB format", type=int, default=42)
    args = parser.parse_args(argv)

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def subsample(input_bam, output_bam, fraction):
    inf = pysam.AlignmentFile(input_bam, "rb")
    outf = pysam.AlignmentFile(output_bam, "wb", template=inf)
    print("Reading " + input_bam)

    count = 0
    subsample_count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        val = random.random()
        if val < fraction:
            outf.write(read)
            subsample_count += 1
    outf.close()
    print("Saved %d reads to %s" % (subsample_count, output_bam))


def main(argv):
    args = parse_args(argv)
    random.seed(args.seed)
    if args.fraction <=0 or args.fraction >= 1:
        print("Set fraction between 0 and 1")
        exit(-1)
    subsample(args.bam, args.output, args.fraction)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
