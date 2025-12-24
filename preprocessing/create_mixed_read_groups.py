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

    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def create_group_files(args):
    inf = pysam.AlignmentFile(args.bam, "rb")
    split_outs = {}
    file_group_ids = inf.references
    for i in file_group_ids:
        split_outs[i] = pysam.AlignmentFile(args.output + ".group_%s.bam" % i, "wb", template=inf)
    tsv_out = open(args.output + ".groups.tsv", "w")
    csv_out = open(args.output + ".groups.csv", "w")


    print("Reading " + args.bam)
    processed_reads = set()

    for read in inf:
        read_id = read.query_name
        group_index = read.reference_name
        if not group_index: continue
        read_group = read_id.split('_')[0]

        new_read = copy.deepcopy(read)
        new_read.query_name += "_" + read_group
        new_read .set_tag("RG", read_group, value_type='Z')
        new_read .set_tag("CB", read_group, value_type='Z')
        read_id = new_read.query_name

        split_outs[group_index].write(new_read)
        if read_id not in processed_reads:
            tsv_out.write("%s\t%s\t%s\n" % (read_id, read_group, group_index))
            csv_out.write("%s,%s,%s\n" % (read_id, read_group, group_index))
            processed_reads.add(read_id)

    tsv_out.close()
    csv_out.close()

    for f in split_outs.values():
        f.close()
        pysam.index(f.filename.decode('utf-8'))
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
