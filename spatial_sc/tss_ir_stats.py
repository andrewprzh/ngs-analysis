#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
import os
import sys
import pysam
from collections import defaultdict
from enum import Enum
import gffutils


def read_ids(inf):
    total_reads = 0
    ids = set()
    for l in open(inf):
        total_reads += 1
        v = l.strip().split('\t')
        ids.add(v[0])
        if total_reads % 100000 == 0:
            sys.stdout.write("Processed " + str(total_reads) + " reads\r")

    print("Loaded %d read ids" % len(ids))
    return ids


def process_read_assignments(read_assignments, ids=None):
    stat_dict = defaultdict(int)
    processed_reads = set()
    for l in open(read_assignments):
        if l.startswith("#"): continue

        v = l.split("\t")
        read_id = v[0]
        if read_id in processed_reads or (ids is not None and read_id not in ids):
            continue

        processed_reads.add(read_id)
        attributes = v[6].split(",")
        if any(attr.startswith("tss_match") for attr in attributes):
            stat_dict["tss_match"] += 1
        if any(attr.startswith("tss_match_precise") for attr in attributes):
            stat_dict["tss_match_precise"] += 1

        if any(attr.startswith("incomplete_intron_retention") for attr in attributes):
            stat_dict["incomplete_intron_retention"] += 1
        if any(attr.startswith("intron_retention") for attr in attributes):
            stat_dict["intron_retention"] += 1
        if any(attr.startswith("unspliced_intron_retention") for attr in attributes):
            stat_dict["unspliced_intron_retention"] += 1

        if len(processed_reads) % 100000 == 0:
            sys.stdout.write("Processed " + str(len(processed_reads)) + " reads\r")

    print("Total reads\t%d" % len(processed_reads))
    for k in sorted(stat_dict.keys()):
        print("%s\t%d" % (k, stat_dict[k]))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--output", "-o", type=str, help="output folder (same as input file if not set), "
    #                                                      "or file name")
    parser.add_argument("--allinfo", type=str, help="allinfo file or read id list")
    parser.add_argument("--ra", type=str, help="read assignments", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    ids = None
    if args.allinfo:
        ids = read_ids(args.allinfo)

    process_read_assignments(args.ra, ids)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
