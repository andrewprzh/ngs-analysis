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
from collections import namedtuple

import pysam
from Bio import SeqIO
import gffutils
from enum import Enum
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


def select_assignments(assignment_file, output_file, chr_id=None, assignment_type=None,
                       assignment_subtype=None, isoform_set=None, ignored_reads=None, forse_equal=True):
    outf = open(output_file, 'w')
    for l in open(assignment_file, 'r'):
        if l.startswith('#'):
            continue
        t = l.strip().split()
        read_id = t[0]
        chr = t[1]
        isoform_id = t[3]
        atype = t[5]
        asubtype = t[6]
        exons = t[7]

        if chr_id and chr != chr_id:
            continue
        if isoform_set and isoform_id not in isoform_set:
            continue
        if ignored_reads and read_id in ignored_reads:
            continue
        if assignment_type and atype != assignment_type:
            continue
        if assignment_subtype and asubtype.find(assignment_subtype) == -1:
            continue
        if forse_equal and read_id.find(isoform_id) == -1:
            continue

        outf.write("%s\t%s\t%s\n" % (exons, isoform_id, read_id))

    outf.close()


def get_set_id(inf, column=1):
    id_set = set()
    for l in open(inf, 'r'):
        t = l.strip().split()
        if len(t) <= column:
            continue
        id_set.add(t[column])
    return id_set


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", default="exmaples", required=True)
    parser.add_argument("--isoquant_assignments", "-i", type=str, help="path to IsoQuant output", required=True)
    parser.add_argument("--chr", type=str, help="chromosome id")
    parser.add_argument("--assignment_type", type=str, help="only use assignments of this type")
    parser.add_argument("--assignment_subtype", type=str, help="only use assignments of this subtype")
    parser.add_argument("--isoform_set", type=str, help="limit isoforms to this set (second column)")
    parser.add_argument("--ignored_reads", type=str, help="ignore these reads (third column)")
    parser.add_argument("--no_equal", help="do not force read to be originated from the assigned isoform",
                        action='store_true', default=False)

    args = parser.parse_args()
    return args


def main():
    set_logger(logger)
    args = parse_args()
    isoform_set = None
    if args.isoform_set:
        isoform_set = get_set_id(args.isoform_set, column=1)
    ignored_reads = None
    if args.ignored_reads:
        ignored_reads = get_set_id(args.ignored_reads, column=2)
    select_assignments(args.isoquant_assignments, args.output, args.chr, args.assignment_type, args.assignment_subtype,
                       isoform_set, ignored_reads, not args.no_equal)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
