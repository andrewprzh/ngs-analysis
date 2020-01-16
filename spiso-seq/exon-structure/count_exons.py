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
from assign_isoforms_to_barcodes import *
from traceback import print_exc


def read_property_table(property_file, column):
    table = {}
    for l in open(property_file):
        tokens = l.strip().split()
        table[tokens[0].replace('@','')] = tokens[column]
    return table


class ReadIdPropertyGetter:
    def get_property(self, read_id):
        return


class TablePropertyGetter:
    def __init__(self, table):
        self.missed = 0
        self.table = table

    def get_property(self, read_id):
        if read_id not in self.table:
            self.missed += 1
            return "unknown"
        return self.table[read_id]


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', metavar='BAM_FILE', nargs='+', type=str, help='sorted and indexed BAM file(s)')
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    parser.add_argument("--readmap", "-r", help="read property table (first column - read id, other - property)", type=str)
    parser.add_argument("--property_column", help="0-based read property column index (0th column is read id)", type=int, default=2)
    parser.add_argument("--barcode_column", help="0-based barcode column index (0th column is read id)", type=int, default=1)
    parser.add_argument("--keep_terminal", help="do not suppress terminal", action='store_true', default=False)

    args = parser.parse_args()
    return args

#Tune algorithm params
def set_codon_count_params(args):
    args.reads_cutoff = 0
    args.data_type = "long_reads"
    args.change_chr_prefix = False
    args.assign_codons_when_ambiguous = False
    args.consider_flanking_junctions = False
    args.junction_delta = 1
    args.count_isoform_stats = False
    args.exon_count_mode = True


def main():
    args = parse_args()
    set_codon_count_params(args)

    property_getter = ReadIdPropertyGetter()
    args.read_info_map = None
    if args.readmap:
        table = read_property_table(args.readmap, args.property_column)
        property_getter = TablePropertyGetter(table)
        args.read_info_map = read_property_table(args.readmap, args.barcode_column)
        args.distinct_barcodes = float(len(set(args.read_info_map.values())))
        print("Barcode map loaded, total distinct barcodes " + str(args.distinct_barcodes))

    db_processor = GeneDBProcessor(args, property_getter)
    db_processor.process_all_genes()
    print(property_getter.missed)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
