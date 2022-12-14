############################################################################
# Copyright (c) 2022 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
import glob
from traceback import print_exc
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="Collect delta-PSI table for the heatmap.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--input", help="output of exon testing, can be multiple files with wildcards", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output file", type=str, required=True)
    required_group.add_argument("--pval", help="p-value threshold (max)", type=float, default=1.0)
    required_group.add_argument("--fdr", help="FDR threshold (max)", type=float, default=1.0)
    required_group.add_argument("--dpsi", help="delta PST threshold (min)", type=float, default=0.3)

    args = parser.parse_args()

    return args


def process_single_file(inf, args, table):
    for l in open(inf):
        if l.startswith("GE"):
            continue
        v = l.strip().split('\t')
        pval = float(v[1])
        dpsi = float(v[2])
        fdr = float(v[5])
        if pval > args.pval or fdr > args.fdr:
            continue
        gene_name = v[6].split('-')[0]
        exon_id = v[0]
        table[gene_name][exon_id] = dpsi


def collect_table(args):
    table = defaultdict(dict)
    print("Reading input file(s) " + args.input)
    for f in glob.glob(args.input):
        process_single_file(f, args, table)

    all_exon_ids = set()
    for gene_name in table.keys():
        all_exon_ids.update(table[gene_name].keys())

    outf = open(args.output, "w")
    gene_names = sorted(table.keys())
    print("Total %d genes detected" % len(gene_names))
    outf.write("exon_id\t%s\n" % "\t".join(gene_names))
    exon_count = 0
    for exon_id in all_exon_ids:
        max_dpsi = 0
        dpsi_values = []
        for gene_name in gene_names:
            if exon_id in table[gene_name]:
                dpsi_values.append("%.3f" % table[gene_name][exon_id])
                max_dpsi = max(max_dpsi, abs(table[gene_name][exon_id]))
            else:
                dpsi_values.append("NA")
        if max_dpsi >= args.dpsi:
            outf.write("%s\t%s\n" % (exon_id, "\t".join(dpsi_values)))
            exon_count += 1
    print("Table with %d exons is written to %s" % (exon_count, args.output))
    outf.close()


def main():
    args = parse_args()
    collect_table(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
