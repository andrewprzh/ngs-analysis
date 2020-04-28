############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import argparse
import pysam
from common import *
from traceback import print_exc


def generate_counts(db, count):
    rpm_map = {}
    for t in db.features_of_type('transcript'):
        transcript_len = 0
        for e in db.children(t, featuretype='exon', order_by='start'):
            transcript_len += e.end + 1 - e.start
        rpm_map[t.id] = float(count * 1000.0) / float(transcript_len)
    return rpm_map


def print_counts(rpm_map, counts, fname):
    rpm_sum = 0.0
    for t in rpm_map.keys():
        rpm_sum += rpm_map[t]
    scaling_factor = rpm_sum / 1000000.0

    outf = open(fname, "w")
    outf.write("target_id\test_counts\ttpm\n")
    for t in rpm_map:
        outf.write(t + "\t" + str(counts) + "\t" + "{0:.8f}".format(rpm_map[t] / scaling_factor) + "\n")
    outf.close()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output", "-o", help="output file", type=str)
    parser.add_argument("--counts", "-c", help="counts per each transcript", type=int, default=10)
    args = parser.parse_args()

    if args.genedb is None or args.output is None:
        print("Provide gene db and output")
        exit(-1)
    return args


def main():
    args = parse_args()

    db = gffutils.FeatureDB(args.genedb, keep_order=True)
    rpm_map = generate_counts(db, args.counts)
    print_counts(rpm_map, args.counts, args.output)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
