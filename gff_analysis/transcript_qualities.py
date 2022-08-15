############################################################################
# Copyright (c) 2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
import numpy


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for adding cellranger tags into BAM for simulated data.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--log", "-l", help="IsoQuant log file", type=str, required=True)
    required_group.add_argument("--gffcompare", "-g", help="gffcompare tracking file", type=str, required=True)
    args = parser.parse_args()

    return args


def parse_log(logfile):
    # 2022-07-11 21:45:35,844 - INFO - Novel model transcript_7209.chr6.nnic has quality 35.57
    qmap = {}
    for l in open(logfile):
        if l.find("has quality") == -1:
            continue

        v = l.strip().split()
        if len(v) < 5:
            continue
        t_id = v[-4]
        q = float(v[-1])
        qmap[t_id] = q
    return qmap


def parse_gffcompare(tracking_file):
    correct_transcripts = set()
    partial_transcripts = set()
    unmapped_transcripts = set()
    for l in open(tracking_file):
        if l.find('|2|') == -1:
            continue
        #TCONS_00000005  XLOC_000004     ENSMUSG00000054493.3|ENSMUST00000067599.2       =       q1:novel_gene_chr1_132|transcript_131.chr1.nnic|1|0.000000|0.000000|0.000000|2920
        v = l.strip().split('\t')
        status = v[3]
        t_id = v[4].split('|')[1]
        if status == '=':
            correct_transcripts.add(t_id)
        elif status == 'u':
            unmapped_transcripts.add(t_id)
        else:
            partial_transcripts.add(t_id)
    return correct_transcripts, partial_transcripts, unmapped_transcripts


def make_hist(qmap, id_set):
    values = []
    for t_id in id_set:
        if t_id in qmap:
            values.append(qmap[t_id])

    bins = [i * 2 for i in range(0, 31)]
    hist, bin_edges = numpy.histogram(values, bins)
    return hist


def main():
    args = parse_args()
    qmap = parse_log(args.log)
    for id_set in parse_gffcompare(args.gffcompare):
        print(make_hist(qmap, id_set))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
