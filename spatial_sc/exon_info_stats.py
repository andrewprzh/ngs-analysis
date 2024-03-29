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
import glob


class ExonInfo:
    def __init__(self, is_cds, contains_start, contains_stop, cds_overlap, whole_codon_count):
        self.is_cds = is_cds
        self.contains_start = contains_start
        self.contains_stop = contains_stop
        self.cds_overlap = cds_overlap
        self.whole_codon_count = whole_codon_count


def load_exon_info(exon_info_files):
    exon_info_dict = {}
    for l in open(exon_info_files):
        if l.startswith("#"): continue

        v = l.strip().split('\t')
        exon = v[0]
        exon_info_dict[exon] = ExonInfo(bool(v[1]), bool(v[2]), bool(v[3]), float(v[4]), bool(v[5]))
    return exon_info_dict


def count_stats(exon_info_dict, in_file):
    exon_count = 0
    cds_count = 0
    overlaps_cds = 0
    whole_codon_count = 0
    whole_codon_cds_count = 0
    start_count = 0
    stop_count = 0
    avg_cds_overlap = 0.0

    for l in open(in_file):
        if l.startswith("#"): continue
        exon_id = l.strip().split()[0]
        exon_info = exon_info_dict[exon_id]
        exon_count += 1
        if exon_info.is_cds:
            cds_count += 1
        if exon_info.whole_codon_count:
            whole_codon_count += 1
            if exon_info.is_cds:
                whole_codon_cds_count += 1
        if exon_info.contains_start:
            start_count += 1
        if exon_info.contains_stop:
            stop_count += 1
        if exon_info.cds_overlap > 0:
            overlaps_cds += 1
        avg_cds_overlap += exon_info.cds_overlap

    avg_cds_overlap /= exon_count
    return exon_count, cds_count, overlaps_cds, whole_codon_cds_count, whole_codon_cds_count, start_count, stop_count, avg_cds_overlap


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file name", required=True)
    parser.add_argument("--exon_info", type=str, help="TSV with exon info", required=True)
    parser.add_argument("--input", "-i", type=str, nargs='+',
                        help="one or more files/dirs with input exon lists", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    exon_info_dict = load_exon_info(args.exon_info)
    with open(args.output, "w") as outf:
        for inf in args.input:
            if os.path.isdir(inf):
                for f in glob.glob(inf):
                    if os.path.isfile(f):
                        name = os.path.basename(f)
                        outf.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n" % (name, *count_stats(exon_info_dict, f)))
            elif os.path.isfile(inf):
                name = os.path.basename(inf)
                outf.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n" % (name, *count_stats(exon_info_dict, inf)))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
