############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
from lr_snp_caller import *
from filter_snps import *


def common_snps(snp_storages):
    all_chromosomes = set()
    for storage in snp_storages:
        all_chromosomes.update(storage.snp_map.keys())

    total_counts = [0 for i in range(len(snp_storages))]
    intersect_map = {}
    stat_map = {}
    for chr_id in all_chromosomes:
        stat_map[chr_id] = {}
        all_positions = set()

        for i in range(len(snp_storages)):
            if chr_id in snp_storages[i].snp_map:
                all_positions.update(snp_storages[i].snp_map[chr_id].keys())
                total_counts[i] += len(snp_storages[i].snp_map[chr_id].keys())

        for pos in all_positions:
            stat_map[chr_id][pos] = [0 for i in range(len(snp_storages))]
            for i in range(len(snp_storages)):
                if chr_id in snp_storages[i].snp_map and pos in snp_storages[i].snp_map[chr_id]:
                    stat_map[chr_id][pos][i] = 1

            appearance_tuple = tuple(stat_map[chr_id][pos])
            if appearance_tuple not in intersect_map:
                intersect_map[appearance_tuple] = 0
            intersect_map[appearance_tuple] += 1

    return intersect_map, total_counts


def freq_stat(snp_storage):
    pass


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="Get stats for several SNP files ")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument('--tsv', dest='tsv_file', nargs='+', help='list of TSV files')
    #required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str)

    args = parser.parse_args()

    if args.tsv_file is None or len(args.tsv_file) == 0: # or args.output_prefix is None:
        parser.print_help()
        exit(-1)
    return args


# Tune algorithm params
def set_params(args):
    pass

def main():
    args = parse_args()
    # set_params(args)

    snp_storages = []
    for f in args.tsv_file:
        print("Reading from " + f)
        reader = TSVParser(f, [0, 1, 2], args, no_filter=True)
        snp_storages.append(SNPStorage())
        reader.fill_map(snp_storages[-1])

    print(common_snps(snp_storages))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
