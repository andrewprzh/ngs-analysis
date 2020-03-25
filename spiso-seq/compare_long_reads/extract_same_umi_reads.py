############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import copy
import argparse
from traceback import print_exc
from functools import partial


def get_read_info(inf, args):
    read_info_map = {}
    for l in open(inf):
        tokens = l.strip().split()
        if len(tokens) <= max(args.umi_column, args.barcode_column):
            continue
        read_id = tokens[0]
        if read_id.startswith('@'):
            read_id = read_id[1:]
        read_info_map[read_id] = (tokens[args.barcode_column], tokens[args.umi_column])

    return read_info_map


def invert_read_infos(read_info_map, inverted_read_info = {}, index = 0):
    for read_id in read_info_map.keys():
        info = read_info_map[read_id]
        barcode = info[0]
        umi = info[1]
        if barcode not in inverted_read_info:
            inverted_read_info[barcode] = {}
        if umi not in inverted_read_info[barcode]:
            inverted_read_info[barcode][umi] = ([], [])
        inverted_read_info[barcode][umi][index].append(read_id)
    return inverted_read_info


def intersect_read_infos(read_info_map1, read_info_map2):
    #barcode -> umi -> (read_ids, read_ids)
    inverted_read_info = invert_read_infos(read_info_map1)
    print("Inverted map 1")
    inverted_read_info = invert_read_infos(read_info_map2, inverted_read_info, 1)
    print("Inverted map 2")

    return inverted_read_info


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--read_info", nargs=2, help="file with UMIs and barcodes ", type=str)
    parser.add_argument("--barcode_column", help="barcode column", type=int, default=4)
    parser.add_argument("--umi_column", help="umi column", type=int, default=7)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    args = parser.parse_args()

    if args.read_info is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()

    read_info_map1 = get_read_info(args.read_info[0], args)
    print("Read " + str(len(read_info_map1)) + " reads from " + args.read_info[0])
    read_info_map2 = get_read_info(args.read_info[1], args)
    print("Read " + str(len(read_info_map2)) + " reads from " + args.read_info[1])

    inverted_read_info = intersect_read_infos(read_info_map1, read_info_map2)


    print("Printing equal Barcode-UMI pairs")
    found_in1 = 0
    found_in2 = 0
    more_that_one_read = 0
    found_in_both = 0
    outf = open(args.output_prefix + ".same_umi_reads.tsv", 'w')
    for barcode in inverted_read_info.keys():
        for umi in inverted_read_info[barcode].keys():
            read_ids = inverted_read_info[barcode][umi]
            if len(read_ids[0]) == 0:
                found_in2 += 1
            elif len(read_ids[1]) == 0:
                found_in1 += 1
            elif len(read_ids[0]) > 1 and len(read_ids[1]) > 1:
                more_that_one_read += 1
            else:
                found_in_both += 1
                outf.write(barcode + '\t' + umi + '\t' + read_ids[0][0] + '\t' + read_ids[1][0] + '\n')
    outf.close()

    print("Barcode-UMI pairs:")
    print("Found in both " + str(found_in_both))
    print("Found in 1 " + str(found_in1))
    print("Found in 2 " + str(found_in2))
    print("More than one read " + str(more_that_one_read))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
