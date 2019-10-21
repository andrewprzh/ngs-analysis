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


def get_intersected_snps(snp_storages, storage_index, cov_cutoff = 0, freq_cutoff = 0):
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
            if any(snp_storages[storage_index].snp_map[chr_id][pos][0].sample_coverage[i] < cov_cutoff for i in range(3)):
                continue

            pos_id = chr_id + ":" + str(pos)
            freq_list = []
            sample_count = len(snp_storages[storage_index].snp_map[chr_id][pos][0].sample_coverage)
            for i in range(sample_count):
                snp_cov = snp_storages[storage_index].snp_map[chr_id][pos][0].sample_counts[i]
                total_cov = snp_storages[storage_index].snp_map[chr_id][pos][0].sample_coverage[i]
                freq_list.append(float(snp_cov) / float(total_cov))

            if any(freq >= freq_cutoff for freq in freq_list) and  any(freq < 0.05 for freq in freq_list):
                snp_frequency_map[pos_id] = freq_list

    return snp_frequency_map


def print_map(snp_frequency_map):
    for pos in sorted(snp_frequency_map.keys()):
        print(pos + '\t' + '\t'.join(map('{:.2f}'.format, snp_frequency_map[pos])))


def freq_stat(snp_frequency_map):
    key = None
    values = []
    for i in range(3):
        freqs = [snp_frequency_map[pos][i] for pos in snp_frequency_map.keys()]
        v,k = numpy.histogram(freqs, bins = [0.1 * i for i in range(11)])
        key = k
        values.append(v)

    for i in range(len(v)):
        print('{:.2f}'.format(key[i]) + '\t' + '\t'.join(map(str, [values[j][i] for j in range(3)])))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="Get stats for several SNP files ")

    required_group = parser.add_argument_group('parameters')
    required_group.add_argument('--tsv', dest='tsv_file', nargs='+', help='list of TSV files')
    required_group.add_argument("--min_freq", help="absolute frequency cutoff, between 0.0 and 1.0 [0.0]", type=float, default=0.0)
    required_group.add_argument("--min_cov", help="absolute coverage cutoff, >= 0 [0]", type=int, default=0)
    required_group.add_argument("--tool_id", help="tools id taken for further analysis", type=int, default=2)

    args = parser.parse_args()

    if args.tsv_file is None or len(args.tsv_file) == 0: # or args.output_prefix is None:
        parser.print_help()
        exit(-1)
    return args


# Tune algorithm params
def set_params(args):
    if args.min_freq < 0 or args.min_freq > 1:
        raise Exception("ERROR: minimal SNP frequency should be between 0.0 and 1.0, but set to " + str(args.min_freq))
    if args.min_cov < 0:
        raise Exception("ERROR: minimal SNP coverage should be positive, but set to " + str(args.min_cov))

def main():
    args = parse_args()
    set_params(args)

    snp_storages = []
    for f in args.tsv_file:
        print("Reading from " + f)
        reader = TSVParser(f, [0, 1, 2], args, no_filter=True)
        snp_storages.append(SNPStorage())
        reader.fill_map(snp_storages[-1])

    print(common_snps(snp_storages))

    snp_freqs = get_intersected_snps(snp_storages, args.tool_id, args.min_cov, args.min_freq)
    print_map(snp_freqs)

    freq_stat(snp_freqs)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
