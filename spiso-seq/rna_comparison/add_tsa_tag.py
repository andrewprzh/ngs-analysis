#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc

import logging
import pysam

logger = logging.getLogger('IsoQuantQA')

def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def load_tags(inbam):
    tag_dict = {}
    for a in pysam.AlignmentFile(inbam, "rb"):
        tag_dict[a.query_name] = a.get_tag("ts")
    return tag_dict


def save_reads_to_bam(input_fname, output_fname, tags):
    bamfile_in = pysam.AlignmentFile(input_fname, "rb")
    bamfile_out = pysam.AlignmentFile(output_fname, "wb", template=bamfile_in)
    for alignment in bamfile_in:
        if alignment.query_name in tags:
            alignment.set_tag("ts", tags[alignment.query_name], value_type='A')
        bamfile_out.write(alignment)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output BAM name [default = input.tsa.bam", default="")
    parser.add_argument("--original_bam", "-g", required=True, type=str, help="BAM file with TS:A tag")
    parser.add_argument("--bam", "-b", required=True, type=str, help="input BAM file to fix")

    args = parser.parse_args()
    return args


def main():
    set_logger(logger)
    args = parse_args()
    tags = load_tags(args.original_bam)
    print(tags)
    if not args.output:
        base, ext = os.path.splitext(args.bam)
        args.output = base + ".tagged.bam"
    save_reads_to_bam(args.bam, args.output, tags)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



