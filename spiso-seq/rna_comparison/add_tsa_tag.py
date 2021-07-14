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
    tag_dict = defaultdict(dict)
    for a in pysam.AlignmentFile(inbam, "rb"):
        try:
            read_id = a.query_name
            if read_id in tag_dict[a.reference_name]:
                continue
            tag_dict[a.reference_name][read_id] = a.get_tag("ts")
        except KeyError:
            continue
    return tag_dict


def load_tags_from_tsv(intsv):
    tag_dict = defaultdict(dict)
    for l in open(intsv):
        if l.startswith('#'):
            continue
        t = l.split('\t')
        read_id = t[0]
        chr_id = t[1]
        strand = t[2]
        if read_id in tag_dict[chr_id]:
            continue
        if strand != '.':
            tag_dict[chr_id][read_id] = strand
    return tag_dict


def save_reads_to_bam(input_fname, output_fname, tags, invert_strand=False):
    bamfile_in = pysam.AlignmentFile(input_fname, "rb")
    bamfile_out = pysam.AlignmentFile(output_fname, "wb", template=bamfile_in)
    for alignment in bamfile_in:
        if alignment.reference_name in tags and alignment.query_name in tags[alignment.reference_name]:
            tag = tags[alignment.reference_name][alignment.query_name]
            if tag == '.':
                continue
            if invert_strand and alignment.is_reverse:
                tag = '-' if tag == '+' else '+'
            alignment.set_tag("ts", tag, value_type='A')
        bamfile_out.write(alignment)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output BAM name [default = input.tsa.bam", default="")
    parser.add_argument("--original_bam", "-g", type=str, help="BAM file with TS:A tag")
    parser.add_argument("--read_assignments", "-r", type=str, help="IsoQuant read assignments in TSV")
    parser.add_argument("--bam", "-b", required=True, type=str, help="input BAM file to fix")

    args = parser.parse_args()
    return args


def main():
    set_logger(logger)
    args = parse_args()
    logger.info("Loading tags")
    invert_strand = False
    if args.original_bam:
        tags = load_tags(args.original_bam)
    elif args.read_assignments:
        tags = load_tags_from_tsv(args.read_assignments)
        invert_strand = True
    else:
        logger.error("Tags were not provided")
        return -1
    if not args.output:
        base, ext = os.path.splitext(args.bam)
        args.output = base + ".tagged.bam"
    logger.info("Convertin bam to %s" % args.output)
    save_reads_to_bam(args.bam, args.output, tags, invert_strand)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



