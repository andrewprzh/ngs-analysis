############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import argparse
import pysam
from traceback import print_exc
from collections import defaultdict
import gffutils


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for subsampling BAM files.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--genedb", help="gffutils DB to convert ids (will not convert if not given)"
                                , type=str, required=False)
    required_group.add_argument("--isoform_col", help="isoform column in read id", type=int, default=0)
    required_group.add_argument("--output", "-o", help="output TSV", type=str, required=True)
    args = parser.parse_args(argv)

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def get_counts(input_bam, isoform_col):
    inf = pysam.AlignmentFile(input_bam, "rb")
    print("Reading " + input_bam)
    isoform_counts = defaultdict(int)
    count = 0

    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.is_secondary or read.is_supplementary:
            continue
        isoform_id = read.query_name.split("_")[isoform_col]
        isoform_counts[isoform_id] += 1

    return isoform_counts


def get_isoform_dict(genedb):
    db = gffutils.FeatureDB(genedb)
    isoform_dict = {}
    for t in db.features_of_type(featuretype=('transcript', 'mRNA')):
        isoform_dict[t.id.split('.')[0]] = t.id
    return isoform_dict


def dump_counts(count_dict, out_fname, isoform_dict):
    scale_factor = sum(count_dict.values()) / 1000000.0
    with open(out_fname, "w") as outf:
        for tid in sorted(count_dict.keys()):
            isofrom_id = isoform_dict[tid] if tid in isoform_dict else tid
            outf.write("%s\t%.2f\t%.6f\n" % (isofrom_id, count_dict[tid], count_dict[tid] / scale_factor))


def main(argv):
    args = parse_args(argv)
    isoform_dict = {}
    if args.genedb:
        isoform_dict = get_isoform_dict(args.genedb)
    count_dict = get_counts(args.bam, args.isoform_col)
    dump_counts(count_dict, args.output, isoform_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
