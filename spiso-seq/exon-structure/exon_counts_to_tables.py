############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
from traceback import print_exc

common_cell_types = ["Astro", "ExcitNeuron", "InhibNeuron", "Vasc"]

# returns exon_id -> {cell_group_id -> inclusion/exclusion counts}
def read_exon_counts(exon_count_file_name):
    exon_counts_map = {}
    for l in open(exon_count_file_name):
        tokens = l.strip().split()
        if tokens[0] not in exon_counts_map:
            exon_counts_map[tokens[0]] = {}
        if tokens[1] in exon_counts_map[tokens[0]]:
            print("Duplicating exon")
        exon_counts_map[tokens[0]][tokens[1]] = (int(tokens[2]), int(tokens[3]))

    return exon_counts_map


def process_exon(single_exon_counts):
    exon_count_tables = {}
    for cell_group in single_exon_counts.keys():
        cell_info = cell_group.split("_")
        cell_type = cell_info[1]

        cell_type_is_known = False
        for ct in common_cell_types:
            if cell_type.startswith(ct):
                cell_type = ct
                cell_type_is_known = True
                break
        if not cell_type_is_known:
            continue

        if cell_type not in exon_count_tables:
            exon_count_tables[cell_type] = {"P7Hipp" : [0, 0], "P7PFC" : [0, 0]}
        tissue = cell_info[0]
        print(cell_type, tissue)
        exon_count_tables[cell_type][tissue][0] += single_exon_counts[cell_group][0]
        exon_count_tables[cell_type][tissue][1] += single_exon_counts[cell_group][1]

    return exon_count_tables


def write_exon_count_table(exon_id, exon_count_tables, outf):
    for ct in exon_count_tables.keys():
        outf.write("====" + exon_id + "_" + ct + "\n")
        outf.write("\t".join(map(str, exon_count_tables[ct]["P7Hipp"])) + "\n")
        outf.write("\t".join(map(str, exon_count_tables[ct]["P7PFC"])) + "\n")


def process_all_counts(exon_counts_map, out_file_name):
    outf = open(out_file_name, "w")
    for exon_id in exon_counts_map.keys():
        exon_count_tables = process_exon(exon_counts_map[exon_id])
        write_exon_count_table(exon_id, exon_count_tables, outf)
    outf.close()


def main():
    if len(sys.argv) != 3:
        print("Usage: " + sys.argv[0] + " <exon counts file> <output>")
        return

    exon_count_map = read_exon_counts(sys.argv[1])
    process_all_counts(exon_count_map, sys.argv[2])


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


