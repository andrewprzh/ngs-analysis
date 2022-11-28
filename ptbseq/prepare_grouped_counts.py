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


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", default="guides_output")
    parser.add_argument("--grouped_exon_counts", "-e", type=str, help="IsoQuant barcode grouped exon counts")
    parser.add_argument("--grouped_transcript_counts", "-t", type=str, help="IsoQuant barcode grouped isoform counts")
    parser.add_argument("--guide_info", "-g", type=str, help="Table with guides / factors", required=True)
    parser.add_argument("--guide_columns", type=str, help="comma-separated guide columns, first one is the guide column", default="1,4,2")
    parser.add_argument("--barcode2guide", "-b", type=str, help="Table with barcodes and guides", required=True)

    args = parser.parse_args()
    return args


def parse_barcode_table(guides, barcodes, columns):
    cols = list(map(int, columns.split(',')))
    guide_column = cols[0]
    min_cols = max(cols)
    cols = cols[1:]
    guide_dict = {}
    for l in open(guides):
        if l.startswith("#"):
            continue
        v = l.strip().split('\t')
        if len(v) <= min_cols:
            continue
        guide_dict[v[guide_column]] = [v[col].replace(" ", "_") for col in cols]

    barcode_dict = {}
    for l in open(barcodes):
        v = l.strip().split('\t')
        barcode = v[0]
        guide = "_".join(v[1].split("_")[:3])
        group = "CROPseq::" + "::".join(guide_dict[guide]) + "::" + guide + "::" + barcode
        barcode_dict[barcode] = group

    return barcode_dict


def convert_exon_counts(exon_counts, barcode_dict, output_file):
    outf = open(output_file, "w")
    for l in open(exon_counts):
        if l.startswith("#"):
            continue
        v = l.strip().split('\t')
        exon_id = "_".join([v[0], v[1], v[2], v[5], v[3]])
        barcode = v[6]
        inc_count = int(v[7])
        exc_count = int(v[8])
        if barcode in barcode_dict:
            group = barcode_dict[barcode]
            outf.write("%s\t%d\t%d\t%d\t%s\n" % (exon_id, inc_count, exc_count, inc_count + exc_count, group))
    outf.close()


def convert_isoform_counts(isoform_counts, barcode_dict, output_file):
    outf = open(output_file, "w")
    bc_list = []
    for l in open(isoform_counts):
        if l.startswith("#"):
            if not bc_list:
                bc_list = l.strip().split('\t')[1:]
            continue

        v = l.strip().split('\t')
        isoform_id = v[0]
        values = list(map(float, v[1:]))
        assert len(values) == len(bc_list)
        for i in range(len(bc_list)):
            barcode = bc_list[i]
            if barcode in barcode_dict and values[i] > 0:
                group = barcode_dict[barcode]
                outf.write("%s\t%.2f\t%s\n" % (isoform_id, values[i], group))
    outf.close()

def main():
    #set_logger(logger)
    args = parse_args()



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



