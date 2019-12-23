############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import argparse
import pysam
from traceback import print_exc

MIN_UMIS = 2
MIN_CELLS = 100
MIN_FRAC = 0.1

class ReadInfo:
    def __init__(self, UMI, isoform):
        self.UMI = UMI
        self.isoform = isoform


def add_to_map(map, key):
    if key not in map:
        map[key] = {}


def get_all_info_table(all_info_file):
    # cell_group -> {gene_id -> {barcode -> {read_id -> read_info}}}
    all_read_map = {}
    for l in open(all_info_file):
        tokens = l.strip().split()
        if len(tokens) < 9:
            continue
        isoform = tokens[6]
        if not isoform.startswith('ENSM'):
            continue
        read_id = tokens[0]
        gene_id = tokens[1]
        group_name = tokes[2]
        barcode = tokens[3]
        umi = tokens[4]

        add_to_map(all_read_map, group_name)
        add_to_map(all_read_map[group_name], gene_id)
        add_to_map(all_read_map[group_name][gene_id], barcode)
        if read_id in all_read_map[group_name][gene_id][barcode]:
            print("Duplicated read id " + read_id)
        all_read_map[group_name][gene_id][barcode][read_id] = ReadInfo(umi, isoform)

    print("Loaded " + str(len(all_read_map)) + " cell groups")
    return all_read_map


def process_gene(gene_map):
    isoform_map = {}
    for barcode in gene_map:
        umis = set([read_info.UMI for read_info in gene_map[barcode].values()])
        if len(umis) < MIN_UMIS:
            continue
        isoforms = set([read_info.isoform for read_info in gene_map[barcode].values()])
        if len(isoforms) != 1:
            continue
        isoform_id = list(isoforms)[0]
        if isoform_id not in isoform_map:
            isoform_map[isoform_id] = []
        isoform_map[isoform_id].append(barcode)
    return isoform_map


def has_diff_isoform(isoform_map):
    total_cells = sum(map(len, isoform_map.values()))
    total_isoforms = len(isoform_map.keys())
    if total_isoforms != 2 or total_cells < MIN_CELLS:
        return False

    sorted_isoforms = sorted(isoform_map.items(), reverse=True, key = lambda x: len(x[1]))
    if len(sorted_isoforms[1][1]) >= float(total_cells) * MIN_FRAC:
        return True
    return False


def process_table(all_read_map):
    for cell_type in all_read_map.keys():
        for gene_id in all_read_map[cell_type].keys():
            isoform_map = process_gene( all_read_map[cell_type][gene_id])
            if has_diff_isoform(isoform_map):
                sorted_isoforms = sorted(isoform_map.items(), reverse=True, key=lambda x: len(x[1]))
                print(cell_type + "\t" + gene_id + "\t" + sorted_isoforms[0][0] + "\t" + str(len(sorted_isoforms[0][1]))
                      + "\t" + sorted_isoforms[1][0] + "\t" + str(len(sorted_isoforms[1][1])))
                print("\t".join(sorted_isoforms[0][1]))
                print("\t".join(sorted_isoforms[1][1]))
                            

def main():
    if len(sys.argv) < 2:
        print("Provide read info table")
        sys.exit(-1)

    read_info_table = get_all_info_table(sys.argv)
    process_table(read_info_table)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)











