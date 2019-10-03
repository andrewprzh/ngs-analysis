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
import numpy


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


def get_intersected_snps(snp_storages, storage_index = 0, cov_cutoff = 0):
    common_chromosomes = set(snp_storages[0].snp_map.keys())
    for i in range(1, len(snp_storages)):
        common_chromosomes =  common_chromosomes.intersection(set(snp_storages[i].snp_map.keys()))

    # chr:pos -> list of frequences
    snp_frequency_map = {}
    for chr_id in common_chromosomes:
        common_positions = set(snp_storages[0].snp_map[chr_id].keys())
        for i in range(1, len(snp_storages)):
            common_positions = common_positions.intersection(set(snp_storages[i].snp_map[chr_id].keys()))

        for pos in common_positions:
            if any(len(snp_storages[i].snp_map[chr_id][pos]) > 1 for i in range(len(snp_storages))):
                continue
            if any(snp_storages[i].snp_map[chr_id][pos][0].sample_coverage < cov_cutoff for i in range(len(snp_storages))):
                continue

            pos_id = chr_id + ":" + str(pos)
            snp_frequency_map[pos_id] = []
            sample_count = len(snp_storages[storage_index].snp_map[chr_id][pos][0].sample_coverage)
            for i in range(sample_count):
                snp_cov = snp_storages[storage_index].snp_map[chr_id][pos][0].sample_counts[i]
                total_cov = snp_storages[storage_index].snp_map[chr_id][pos][0].sample_coverage[i]
                snp_frequency_map[pos_id].append(float(snp_cov) / float(total_cov))

    return snp_frequency_map


def print_map(snp_frequency_map):
    for pos in sorted(snp_frequency_map.keys()):
        print(pos + '\t' + '\t'.join(map('{:.2f}'.format, snp_frequency_map[pos])))


def freq_stat(snp_frequency_map, sample_index):
    freqs = [snp_frequency_map[pos][sample_index] for pos in snp_frequency_map.keys()]
    print(numpy.histogram(freqs, bins = [0.1 * i for i in range(10)]))


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

    snp_freqs = get_intersected_snps(snp_storages, 2, 10)
    print_map(snp_freqs)
    for i in range(len(args.tsv_file)):
        freq_stat(snp_freqs, i)

    snp_freqs = get_intersected_snps(snp_storages, 2, 20)
    print_map(snp_freqs)
    for i in range(len(args.tsv_file)):
        freq_stat(snp_freqs, i)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
