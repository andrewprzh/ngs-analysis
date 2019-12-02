############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
from Bio import SeqIO
from lr_snp_caller import *


class SNPStorage:
    def __init__(self):
        # chr_id -> map: position -> list of SelectedSNP
        self.snp_map = {}

    def add(self, chr_id, position, multi_snp):
        if chr_id not in self.snp_map:
            self.snp_map[chr_id] = {}
        if position not in self.snp_map[chr_id]:
            self.snp_map[chr_id][position] = [multi_snp]
        else:
            self.snp_map[chr_id][position].append(multi_snp)


class TSVParser:
    def __init__(self, infile, sample_ids, args):
        self.sample_start_column = 5
        self.infile = infile
        self.sample_ids = sample_ids
        if len(self.sample_ids) == 0:
            f = open(self.infile)
            header = f.readline()
            sample_count = (len(header.strip().split()) - 5) / 3
            self.sample_ids = range(sample_count)
        self.args = args
        self.no_filter = args is None or  args.no_cov_filter

    def fill_map(self, snp_storage):
        header = True
        for l in open(self.infile):
            if header:
                header = False
                continue

            tokens = l.strip().split('\t')
            total_cov = [int(tokens[self.sample_start_column + 3 * i]) for i in self.sample_ids]
            if not self.no_filter and any(cov < self.args.min_cov for cov in total_cov):
                continue
            freqs = [float(tokens[self.sample_start_column + 3 * i + 2]) for i in self.sample_ids]
            if not self.no_filter and max(freqs) < self.args.min_freq:
                continue
            snp_type = GERMLINE_SNP if not self.no_filter and is_germline(min(freqs), max(freqs), self.args) else SOMATIC_SNP
            snp_cov = [int(tokens[self.sample_start_column + 3 * i + 1]) for i in self.sample_ids]

            snp_storage.add(tokens[0], int(tokens[1]), SelectedSNP(tokens[2], tokens[3], snp_type, total_cov, snp_cov))


class VarScanParser:
    def __init__(self, infile, sample_ids, args):
        self.sample_start_column = 10
        self.infile = infile
        self.sample_ids = sample_ids
        self.args = args
        self.no_filter = args.no_cov_filter

    def fill_map(self, snp_storage):
        header = True
        for l in open(self.infile):
            if header:
                header = False
                continue

            tokens = l.strip().split()
            if len(tokens[3]) > 1:
                # ignore multiple variants so far
                continue
            sample_tokens = [tokens[self.sample_start_column + i] for i in self.sample_ids]
            snp_cov = []
            freqs = []
            total_cov = []
            for sample_record in sample_tokens:
                sample_data = sample_record.split(':')
                if len(sample_data) < 5:
                    total_cov.append(0)
                    snp_cov.append(0)
                    freqs.append(0.0)
                else:
                    total_cov.append(0 if sample_data[2] == '-' else int(sample_data[2]) + int(sample_data[3]))
                    snp_cov.append(0 if sample_data[3] == '-' else int(sample_data[3]))
                    freqs.append(0.0 if total_cov[-1] == 0 else float(snp_cov[-1]) / float(total_cov[-1]))

            if not self.no_filter and (any(cov < self.args.min_cov for cov in total_cov) or max(freqs) < self.args.min_freq):
                continue
            snp_type = GERMLINE_SNP if is_germline(min(freqs), max(freqs), self.args) else SOMATIC_SNP
            snp_storage.add(tokens[0], int(tokens[1]), SelectedSNP(tokens[2], tokens[3], snp_type, total_cov, snp_cov))


class VCFParser:
    def __init__(self, infile, sample_ids, args):
        self.sample_start_column = 9
        self.infile = infile
        self.sample_ids = sample_ids
        self.args = args
        self.no_filter = args.no_cov_filter
        self.format_set = False

    def set_sample_data_format(self, format_str):
        tokens = format_str.split(':')
        self.total_fields = len(tokens)
        self.ad_index = tokens.index('AD')
        self.dp_index = tokens.index('DP')
        self.format_set = True

    def fill_map(self, snp_storage):
        for l in open(self.infile):
            if l.startswith("#"):
                continue

            tokens = l.strip().split('\t')
            if len(tokens[4]) > 1:
                # ignore multiple variants so far
                continue

            if not self.format_set:
                self.set_sample_data_format(tokens[self.sample_start_column - 1])

            sample_tokens = [tokens[self.sample_start_column + i] for i in self.sample_ids]
            snp_cov = []
            freqs = []
            total_cov = []
            for sample_record in sample_tokens:
                sample_data = sample_record.split(':')
                if len(sample_data) < self.total_fields:
                    total_cov.append(0)
                    snp_cov.append(0)
                    freqs.append(0.0)
                else:
                    cov = sample_data[self.ad_index].split(',')
                    snp_cov.append(int(cov[1]))
                    dp = int(cov[0]) + int(cov[1])
                    if dp != int(sample_data[self.dp_index]):
                        sys.stderr.write("DP != sum(AD), line = " + l + '\n')
                    total_cov.append(dp)
                    freqs.append(0.0 if total_cov[-1] == 0 else float(snp_cov[-1]) / float(total_cov[-1]))

            if not self.no_filter and (any(cov < self.args.min_cov for cov in total_cov) or
                                       (max(freqs) < self.args.min_freq) and sum(freqs) < self.args.min_freq):
                continue
            snp_type = GERMLINE_SNP if is_germline(min(freqs), max(freqs), self.args) else SOMATIC_SNP
            chr_id = tokens[0]
            if chr_id.startswith('chr'):
                chr_id = chr_id[3:]
            snp_storage.add(chr_id, int(tokens[1]), SelectedSNP(tokens[3], tokens[4], snp_type, total_cov, snp_cov))


class SNPFilter:
    def __init__(self, args):
        self.args = args

    def filter_chromosome(self, chromosome_record, chr_map, filtered_chr_map):
        sorted_positions = sorted(chr_map.keys())
        for i in range(len(sorted_positions)):
            if i > 0 and sorted_positions[i] - sorted_positions[i - 1] <= self.args.min_distance_between_snps:
                continue
            if i < len(sorted_positions) - 1  and \
                    sorted_positions[i + 1] - sorted_positions[i] <= self.args.min_distance_between_snps:
                continue

            passed_snp_list = []
            position = sorted_positions[i]
            snp_region = chromosome_record.seq[max(0, position - 5):position + 5].upper()
            region_len = len(snp_region)

            poly_a_region = float(snp_region.count('A')) / float(region_len) >= self.args.poly_at_percentage
            poly_t_region = float(snp_region.count('T')) / float(region_len) >= self.args.poly_at_percentage
            for snp in chr_map[position]:
                if (poly_a_region and snp.alternative_nucl == 'A') or (poly_t_region and snp.alternative_nucl == 'T') \
                        or (snp.reference_nucl == 'A' and snp.alternative_nucl == 'G') or (snp.reference_nucl == 'T' and snp.alternative_nucl == 'C'):
                    # ignore X->A in polyA, X->T in polyT, A->G
                    continue
                if float(map(lambda x: x >= 1, snp.sample_coverage).count(True)) < self.args.min_frac * float(len(snp.sample_coverage)):
                    continue
                passed_snp_list.append(snp)
            if len(passed_snp_list) > 0:
                filtered_chr_map[position] = passed_snp_list

    def filter(self, snp_storage):
        filtered_storage = SNPStorage()
        ref_dict = SeqIO.to_dict(SeqIO.parse(self.args.reference, "fasta"))
        for chr_id in snp_storage.snp_map.keys():
            if chr_id not in ref_dict:
                continue
            filtered_storage.snp_map[chr_id] = {}
            self.filter_chromosome(ref_dict[chr_id], snp_storage.snp_map[chr_id], filtered_storage.snp_map[chr_id])

        return filtered_storage


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description="Filter out SNPs that are out of specified group, SNP clusters and SNPs in poly-A/T regions ")

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--reference", "-r", help="reference genome used for alignment in FASTA format", type=str)
    required_group.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    optional_group = parser.add_argument_group('optional parameters')
    optional_group.add_argument("--varscan", help="SNPs in VarScan format", type=str)
    optional_group.add_argument("--vcf", help="SNPs in VCF format (monovar/bcftools)", type=str)
    optional_group.add_argument("--tsv", help="SNPs in TSV format (this tool)", type=str)
    optional_group.add_argument("--no_filter", help="do not filter  SNP clusters and SNPs in poly-A/T regions", action='store_true', default=False)
    optional_group.add_argument("--no_cov_filter", help="do not use coverage filters", action='store_true', default=False)
    optional_group.add_argument("--min_freq", "-f", help="minimal SNP frequency within a sample, between 0.0 and 1.0 [0.2]", type=float, default=0.2)
    optional_group.add_argument("--min_freq_factor", "-m", help="minimal SNP frequency factor, > 1.0 [2.0]", type=float, default=2.0)
    optional_group.add_argument("--min_cov", "-c", help="minimal SNP read coverage depth within a sample, > 0 [50]", type=int, default=50)
    optional_group.add_argument("--sample_ids", help="comma-separated 0-based sample indices, e.g. 1,3,6", type=str, default='')
    optional_group.add_argument("--min_frac", help="minimal faction of samples for which SNP is covered, between 0.0 and 1.0 [0.0]", type=float, default=0.0)

    args = parser.parse_args()

    if args.reference is None or args.output_prefix is None:
        parser.print_help()
        exit(-1)
    return args


# Tune algorithm params
def set_params(args):
    if args.min_freq < 0 or args.min_freq > 1:
        raise Exception("ERROR: minimal SNP frequency should be between 0.0 and 1.0, but set to " + str(args.min_freq))
    if args.min_cov < 0:
        raise Exception("ERROR: minimal SNP coverage should be positive, but set to " + str(args.min_cov))
    if args.min_freq_factor < 1:
        raise Exception("ERROR: minimal SNP frequency factor should be larger than 1.0, but set to " + str(args.min_freq))

    if [args.varscan, args.tsv, args.vcf].count(None) != 2:
        raise Exception("ERROR: provide exactly one input file")

    args.min_distance_between_snps = 10
    args.poly_at_percentage = 0.5


def process_sample_ids(sample_ids_string):
    sample_ids = []
    tokens = sample_ids_string.split(',')
    for t in tokens:
        id_range = t.split('-')
        if len(id_range) == 1:
            sample_ids.append(int(id_range[0]))
        else:
            sample_ids.extend(range(int(id_range[0]), int(id_range[1]) + 1))

    return sample_ids


def main():
    args = parse_args()
    set_params(args)

    sample_ids = process_sample_ids(args.sample_ids)
    print(sample_ids)
    snp_writer = SNPMapTSVWriter(args.output_prefix, map(str, sample_ids))

    snp_storage = SNPStorage()
    reader = None
    if args.varscan:
        reader = VarScanParser(args.varscan, sample_ids, args)
    elif args.vcf:
        reader = VCFParser(args.vcf, sample_ids, args)
    elif args.tsv:
        reader = TSVParser(args.tsv, sample_ids, args)
    else:
        raise Exception("ERROR: no input")

    print("Reading from " + reader.infile)
    reader.fill_map(snp_storage)

    if args.no_filter:
        print("No filtering will be applied, writing to " + args.output_prefix)
        snp_writer.dump_to_file(snp_storage.snp_map)
    else:
        print("Filtering SNPs")
        snp_filter = SNPFilter(args)
        filtered_storage = snp_filter.filter(snp_storage)
        print("Writing to " + args.output_prefix)
        snp_writer.dump_to_file(filtered_storage.snp_map)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
