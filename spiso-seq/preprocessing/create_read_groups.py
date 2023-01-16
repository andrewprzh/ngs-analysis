############################################################################
# Copyright (c) 2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import copy
import os
import sys
import argparse
import pysam
import gffutils
from Bio import Seq
from traceback import print_exc
import random


nucls = ['A', 'C', 'G', 'T']


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    required_group.add_argument("--group_num", help="numbre of groups to split", type=int, default=2)

    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def create_group_files(args):
    inf = pysam.AlignmentFile(args.bam, "rb")
    grouped_out = pysam.AlignmentFile(args.output + ".RG.bam", "wb", template=inf)
    group_id_out = pysam.AlignmentFile(args.output + ".read_id.bam", "wb", template=inf)
    split_outs = []
    for i in range(args.group_num):
        split_outs.append(pysam.AlignmentFile(args.output + ".group%d.bam" % i, "wb", template=inf))
    tsv_out = open(args.output + ".groups.tsv", "w")

    print("Reading " + args.bam)
    processed_reads = set()

    for read in inf:
        read_id = read.query_name
        group_index = int(read_id[-2:], 16) % args.group_num
        read_group = "GR" + str(group_index)

        split_outs[group_index].write(read)
        if read_id not in processed_reads:
            tsv_out.write("%s\t%s\n" % (read_id, read_group))
            processed_reads.add(read_id)

        new_read = copy.deepcopy(read)
        new_read.query_name += "_" + read_group
        group_id_out.write(new_read)

        new_read = copy.deepcopy(read)
        new_read .set_tag("RG", read_group, value_type='Z')
        grouped_out.write(new_read)

    tsv_out.close()
    grouped_out.close()
    group_id_out.close()
    for f in split_outs:
        f.close()
    print("Processed %d reads" % len(processed_reads))
    print("Saved to " + args.output)


def main():
    args = parse_args()
    create_group_files(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
