#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import gffutils
import numpy
from pyfaidx import Fasta


def find_cds(gffutils_db, t_id):
    t = gffutils_db[t_id]
    cds = []
    for e in gffutils_db.children(t, featuretype="CDS"):
        cds.append((e.start, e.end))
    if not cds:
        return "NC"

    cds_end_pos = cds[-1][1] if t.strand == "+" else cds[0][0]
    return "%s_%d_%d_%s" % (t.seqid, cds_end_pos, cds_end_pos, t.strand)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output file", required=True)
    parser.add_argument("--input", "-i", type=str, help="input TSV", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    random.seed(args.seed)
    gffutils_db = gffutils.FeatureDB(args.genedb)
    cds_dict = {}
    with open(args.output, "w") as outf:
        for l in open(args.input):
            if l.startswith("#"):
                outf.write(l)
                continue
            l = l.strip()
            t_id = l.split('\t')[1]
            if t_id in cds_dict:
                cds = cds_dict[t_id]
            else:
                cds = find_cds(gffutils_db, t_id)
                cds_dict[t_id] = cds
            outf.write(l + "\t" + str(cds) + "\n")



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
