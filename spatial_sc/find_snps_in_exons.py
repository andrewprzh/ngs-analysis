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
import gzip
from traceback import print_exc
from collections import defaultdict


def load_snps(infile, skip_no_id=True):
    # chr1_14521_C_T_b38      chr1    14521   C       T       1       rs1378626194    1_14521_C_T_b37
    snp_dict = defaultdict(list)
    print("Loading SNPs from %s" % infile)
    for l in open(infile):
        v = l.strip().split('\t')
        id = v[0]
        chr_ids = v[1]
        pos = int(v[2])
        snp_id = v[6]
        if skip_no_id and snp_id == ".": continue

        snp_dict[chr_ids].append((id, pos, snp_id))

    print("Loaded %d SNPs" % sum([len(x) for x in snp_dict.values()]))
    return snp_dict


def load_snp_descriptions(infile):
    # rs11191560      NT5C2   NT5C2   www.ncbi.nlm.nih.gov/pubmed/28443625    Body_mass_index
    print("Loading SNP descriptions from %s" % infile)
    snp_description = {}
    for l in open(infile):
        v = l.strip().split('\t')
        snp_description[v[0]] = v[1:]
    print("Loaded %d SNPs" % len(snp_description))
    return snp_description


def filter_snps(snp_dict, snp_descriptions):
    new_snp_dict = defaultdict(list)
    print("Filtering SNPs")
    for chr_id in snp_dict:
        for snp in snp_dict[chr_id]:
            snp_id = snp[2]
            if snp_id in snp_descriptions:
                new_snp_dict[chr_id].append(snp)

    print("Kept %d SNPs" % sum([len(x) for x in new_snp_dict.values()]))
    return new_snp_dict


def load_exons(tab_file, min_fdr=0.0, max_fdr=1.0):
    # chr3_134547335_134547472_+      ENSG00000182923.18      -0.00874404     1
    print("Loading exons from %s, min FRD %.3f, max FDR %.3f " % (tab_file, min_fdr, max_fdr))
    exon_chr_dict = defaultdict(list)
    for l in open(tab_file):
        v = l.strip().split('\t')
        fdr = float(v[3])
        if fdr > max_fdr or fdr < min_fdr:
            continue

        gene_id = v[1]
        dpsi = float(v[2])
        exon_v = v[0].split('_')
        chr_id = exon_v[0]
        exon_start = int(exon_v[1])
        exon_end = int(exon_v[2])
        strand = exon_v[3]
        exon_chr_dict[chr_id].append((exon_start, exon_end, strand, gene_id, dpsi))

    print("Loaded %d exons" % sum([len(x) for x in exon_chr_dict.values()]))
    return exon_chr_dict


def find_snp(chr_id, snp, exon_chr_dict):
    for exon in exon_chr_dict[chr_id]:
        if exon[0] <= snp[1] <= exon[1]:
            return exon
    return None


def overlap_snps(snp_chr_dict, exon_chr_dict):
    overlapped_snps = []
    for chr_id in snp_chr_dict.keys():
        for snp in snp_chr_dict[chr_id]:
            exon = find_snp(chr_id, snp, exon_chr_dict)
            if exon:
                overlapped_snps.append((chr_id, snp, exon))

    return overlapped_snps


def dump_overlapped_snps(overlapped_snps, snp_description, out_fname):
    print("Saving %d SNPs to %s" % (len(overlapped_snps), out_fname))
    # chr_id exon_start exon_end strand gene_id dpsi snp_position snp_id snp_rs_id snp_description
    with open(out_fname, "w") as outf:
        outf.write("#chr_id\texon_start\texon_end\tstrand\tgene_id\tdpsi\tsnp_position\tsnp_id\tsnp_rs_id\tsnp_description\n")
        for record in overlapped_snps:
            chr_id = record[0]
            snp = record[1]
            exon = record[2]
            (exon_start, exon_end, strand, gene_id, dpsi) = exon
            (snp_name, pos, snp_id) = snp
            description = snp_description[snp_id] if snp_id in snp_description else []
            outf.write("%s\t%d\t%d\t%s\t%s\t%.3f\t%d\t%s\t%s\t%s\n" %
                       (chr_id, exon_start, exon_end, strand, gene_id, dpsi, pos, snp_name, snp_id, " ".join(description)))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--snp", "-s", type=str, required=True,
                        help="reference SNP files with coordinates (will use all SNPs from ref if not provided)")
    parser.add_argument("--snp_description", "-d", type=str, help="SNP description file")
    parser.add_argument("--exons", "-e", type=str,  nargs="+",
                        help="one or more TAB files with exons", required=True)
    parser.add_argument("--min_fdr", type=float, help="minimal FDR for exons (>=)", default=0.0)
    parser.add_argument("--max_fdr", type=float, help="maximal FDR for exons (<=)", default=1.0)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    snp_chr_dict = load_snps(args.snp)
    snp_description = {}
    if args.snp_description:
        snp_description = load_snp_descriptions(args.snp_description)
        snp_chr_dict = filter_snps(snp_chr_dict, snp_description)

    for exon_file in args.exons:
        exon_chr_dict = load_exons(exon_file, args.min_fdr, args.max_fdr)
        exon_name = os.path.splitext(os.path.basename(exon_file))[0]
        overlapped_snps = overlap_snps(snp_chr_dict, exon_chr_dict)
        print("%s, FDR in (%.3f, %.3f), overlapped SNPs: %d" % (exon_name, args.min_fdr, args.max_fdr, len(overlapped_snps)))
        dump_overlapped_snps(overlapped_snps, snp_description, args.output + ("%s_%.3f_%.3f" % (exon_name, args.min_fdr, args.max_fdr)) + ".tsv")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
