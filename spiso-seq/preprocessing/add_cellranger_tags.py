############################################################################
# Copyright (c) 2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
import pysam
import gffutils
from traceback import print_exc

def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="A simple script for adding cellranger tags into BAM for simulated data.")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output BAM", type=str, required=True)
    required_group.add_argument("--genedb", "-g", help="annotation in gffutils DB format", type=str, required=True)
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def read_gene_names(genedb):
    gene_name_dict = {}
    print("Loading " + genedb)
    db = gffutils.FeatureDB(genedb)
    for g in db.features_of_type('gene'):
        gene_name_dict[g.id.split('.')[0]] = g['gene_name'][0]
    return gene_name_dict


def add_cellranger_tags(args, gene_name_dict):
    inf = pysam.AlignmentFile(args.bam, "rb")
    outf = pysam.AlignmentFile(args.output, "wb", template=inf)
    print("Reading " + args.bam)

    count = 0
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        read_id = read.query_name
        vals = read_id.split('_')
        if len(vals) < 4:
            continue
        gene_id = vals[0]
        barcode = vals[2]
        umi = vals[3]

        # GX:Z:ENSG00000243485    GN:Z:MIR1302-2HG        fx:Z:ENSG00000243485
        # RE:A:E  xf:i:0  CR:Z:AAACGCTCACGCAGTC   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAACGCTCACGCAGTC-1 UR:Z:CACATTGTCTCT       UY:Z:FFFFFFFFFFFF
        read.set_tag("GX", gene_id, value_type='Z')
        read.set_tag("GN", gene_name_dict[gene_id], value_type='Z')
        read.set_tag("fx", gene_id, value_type='Z')
        read.set_tag("CR", barcode, value_type='Z')
        read.set_tag("CY", 'F' * len(barcode), value_type='Z')
        read.set_tag("CB", barcode + '-1', value_type='Z')
        read.set_tag("UB", umi, value_type='Z')
        read.set_tag("UR", umi, value_type='Z')
        read.set_tag("UY", 'F' * len(umi), value_type='Z')
        outf.write(read)
    outf.close()
    print("Saved to " + args.output)


def main():
    args = parse_args()
    gene_name_dict = read_gene_names(args.genedb)
    add_cellranger_tags(args, gene_name_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
