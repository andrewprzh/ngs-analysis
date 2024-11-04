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


def amend_bam(in_file_name, out_file_name, main_bam):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    outf = pysam.AlignmentFile(out_file_name, "wb", template=main_bam)

    count = 0
    passed = 0

    main_header = main_bam.header

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.reference_name == "MT":
            new_reference_name = 'chrM'
        elif '.' not in read.reference_name:
            new_reference_name = 'chr' + read.reference_name
        else:
            new_reference_name = read.reference_name

        new_tid = main_header.get_tid(new_reference_name)
        if new_tid != -1:
            read.reference_id = new_tid
            read.reference_name = new_reference_name
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
    parser.add_argument("--master_bam", "-m", type=str, help="BAM to take header from", required=True)
    parser.add_argument("--bam", "-b", nargs="+", type=str, help="BAM file ", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    main_bam = pysam.AlignmentFile(args.master_bam)

    for in_bam in args.bam:
        print("Processing %s" % in_bam)
        if not args.output:
            args.output = os.path.dirname(in_bam)

        if os.path.isdir(args.output):
            out_bam = os.path.join(args.output, os.path.splitext(os.path.basename(in_bam))[0] + "." + args.suffix + ".bam")
        else:
            out_bam = args.output

        amend_bam(in_bam, out_bam, main_bam)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
