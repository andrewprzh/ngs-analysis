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


def read_reads(inf, col):
    read_set = set()
    warn_printed = False
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) > col:
            read_set.add(v[col])
        elif not warn_printed:
            print("Malformed line with %d columns %s" % (len(v), l))
            warn_printed = True
    return read_set


def read_called_barcodes(inf):
    read_set = set()
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) >= 2 and v[1] == "*":
            continue
        read_set.add(v[0])
    print("Loaded %d read ids" % len(read_set))
    return read_set


def has_intron(cigartules):
    return any(x[0] == 3 for x in cigartules)


def filter_reads(in_file_name, out_file_name, read_set=None, spliced=False):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    outf = pysam.AlignmentFile(out_file_name, "wb", template=inf)

    count = 0
    passed = 0

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if (not read_set or read.query_name in read_set) and (not spliced or has_intron(read.cigartuples)):
            outf.write(read)
            passed += 1

    print("Processed " + str(count) + " reads, written " + str(passed))
    inf.close()
    outf.close()
    pysam.index(out_file_name)




def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder (same as input file if not set), "
                                                         "or file name")
    parser.add_argument("--suffix", type=str, help="file suffix (works if output is a folder),"
                                                   " default: filtered", default="filtered")
    parser.add_argument("--list", "-l", type=str, help="read list")
    parser.add_argument("--read_column", type=int, help="read id column (0 based index)", default=0)
    parser.add_argument("--barcodes", default=False, action="store_true",
                        help="read list is barcode files, filter only barcoded ones")
    parser.add_argument("--spliced", default=False, action="store_true",
                        help="keep only spliced reads")

    parser.add_argument("--bam", "-b", nargs="+", type=str, help="BAM file ", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if not args.list:
        read_set = None
    elif args.barcodes:
        read_set = read_called_barcodes(args.list)
    else:
        read_set = read_reads(args.list, args.read_column)

    for in_bam in args.bam:
        print("Processing %s" % in_bam)
        if not args.output:
            args.output = os.path.dirname(in_bam)

        if os.path.isdir(args.output):
            out_bam = os.path.join(args.output, os.path.splitext(os.path.basename(in_bam))[0] + "." + args.suffix + ".bam")
        else:
            out_bam = args.output

        filter_reads(in_bam, out_bam, read_set, spliced=args.spliced)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
