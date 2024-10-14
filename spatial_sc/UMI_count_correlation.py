#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import pysam


# 82718f26-b81f-4561-ac75-59b987941431    ACGCGTTTAAGACG  ATCGGCTAG       14      True    -       82      40      49      66
# LH00376:28:22GVKLLT3:4:1101:2655:1016   41      -1      8       25      +       CAACTGGCCGGGTA  TACCGTCGT       13      False

def get_umi_dict(inf, bam=None, trusted_only=False):
    forbidden_reads = set()
    if bam:
        try:
            for r in pysam.AlignmentFile(bam, "rb").fetch("chrM"):
                forbidden_reads.add(r.query_name)
        except ValueError:
            try:
                for r in pysam.AlignmentFile(bam, "rb").fetch("MT"):
                    forbidden_reads.add(r.query_name)
            except ValueError:
                print("No chrM or MT found")

    discarded = 0
    umi_dict = defaultdict(set)
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split('\t')
        read_id = v[0]
        if read_id in forbidden_reads:
            discarded += 1
            continue
        if v[-1] in ['True', 'False']:
            # old format
            if v[6] == '*' or v[7] == '*': continue
            if trusted_only and v[-1] != "True": continue
            umi_dict[v[6]].add(v[7])
            continue

        if v[1] == '*' or v[2] == '*': continue
        if trusted_only and v[4] != "True": continue
        umi_dict[v[1]].add(v[2])

    print("Discarded %d reads from chrM" % discarded)
    return umi_dict


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file with UMI count tables", required=True)
    parser.add_argument("--input", "-i", nargs="+", type=str, help="input TSV(s) with barcoded reads", required=True)
    parser.add_argument("--bam", "-b", nargs="+", type=str, help="BAM(s) for filtering out chrM, same as barcoded reads (optional)")
    parser.add_argument("--only_trusted_umi", default=False, action="store_true", help="keep only spliced reads")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.bam and len(args.bam) != len(args.input):
        print("Unequal number of TSV and BAM files")
        exit(-1)

    dicts = []
    for i, tsv in enumerate(args.input):
        print("Loading %s" % tsv)
        dicts.append(get_umi_dict(tsv, args.bam[i], args.only_trusted_umi))

    print("Merging barcodes")
    all_barcodes = set()
    for d in dicts:
        all_barcodes.update(d.keys())

    print("Outputting results to %s" % args.output)
    with open(args.output, "w") as outf:
        outf.write("Barcode\t" + "\t".join(args.input) + "\n")
        for b in all_barcodes:
            counts = []
            for d in dicts:
                counts.append(len(d[b]))
            outf.write(b + "\t" + "\t".join(map(str, counts)) + "\n")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
