#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
import logging

def load_allinfo(inf):
    print("Loading allinfo from %s" % inf)
    # layer -> gene -> (sample, read)
    layer_backets = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    count = 0
    for l in open(inf):
        v = l.strip().split("\t")
        layer = v[2]
        gene_id = v[1]
        sample = v[11]
        layer_backets[layer][gene_id][sample].append("\t".join(v[:-1]))
        count += 1
    print("Loaded %d reads and %d layer-age-conditions" % (count, len(layer_backets)))
    return layer_backets


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--allinfo", "-a", nargs='+', type=str, help="allinfo files", required=True)

    args = parser.parse_args()
    return args


def main():
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
