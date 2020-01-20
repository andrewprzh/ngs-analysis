############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
from traceback import print_exc

def filter_exon_counts(args):
    outf = open(args.output, "w")
    for l in open(args.tsv):
        if l.startswith('#'):
            outf.write(l)
            continue

        tokens = l.strip().split()

        exon_type = tokens[5]
        if args.terminal <= 1 and exon_type.find('X') != -1:
            continue
        if args.terminal == 0 and exon_type.find('T') != -1:
            continue
        if args.multiple == 0 and exon_type.find('M') != -1:
            continue

        inclusion_rate = float(tokens[4])
        if inclusion_rate > args.max_inclusion_rate or inclusion_rate < 1.0 - max_inclusion_rate:
            continue
        gene_coverage = float(tokens[7])
        if gene_coverage < args.min_gene_coverage:
            continue
        outf.write(l)
        
    outf.close()

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsv', metavar='TSV_FILE', type=str, help='file with exon counts')
    parser.add_argument("--multiple", "-m", help="1 - keep exons common for different genes, default = 0 - discard common", type=int, default=0)
    parser.add_argument("--terminal", "-t", help="2 - keep all terminal exons, 1 - keep only terminal which are internal as well, default = 0 - discard all terminal", type=int, default=0)
    parser.add_argument("--max_inclusion_rate", "-i", help="max inclusion/exclusion rate, default = 1.0", type=float, default=1.0)
    parser.add_argument("--min_gene_coverage", help="min faction of cells for which this gene has a read, default = 0.0", type=float, default=0.0)
    parser.add_argument("--output", "-o", help="output file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    filter_exon_counts(args)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
