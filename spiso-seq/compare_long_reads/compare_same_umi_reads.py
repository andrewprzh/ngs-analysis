############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import copy
import argparse
import numpy
import pysam
from traceback import print_exc
from functools import partial

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../exon-structure/'))
from common import *


POS_DIFF = 10000


def read_info(info_file):
    barcode_map1 = {}
    barcode_map2 = {}
    read_pairs = []
    for l in open(info_file):
        tokens = l.strip().split()
        if len(tokens) != 4:
            continue
        barcode_map1[tokens[2]] = (tokens[0], tokens[1])
        barcode_map2[tokens[3]] = (tokens[0], tokens[1])
        read_pairs.append((tokens[2], tokens[3]))
    return barcode_map1, barcode_map2, read_pairs


def read_bam_file(infile, barcode_map):
    alignment_map = {}
    for alignment in pysam.AlignmentFile(infile, "rb"):
        if alignment.reference_id == -1:
            continue
        read_id = alignment.query_name
        if read_id in barcode_map:
            if read_id not in alignment_map:
                alignment_map[read_id] = []
            alignment_map[read_id].append(alignment)

    return alignment_map



class AligmentComparator:

    def __init__(self, args):
        self.delta = args.delta

    def compare_junctions(self, junctions1, junctions2):
        pos1 = 0
        pos2 = 0

        features_present1 = [0 for i in range(0, len(junctions1))]
        features_present2 = [0 for i in range(0, len(junctions2))]

        while pos1 < len(junctions1) and pos2 < len(junctions2):
            if equal_ranges(junctions2[pos2], junctions1[pos1], self.delta):
                features_present1[pos1] = 1
                features_present2[pos2] = 1
                pos1 += 1
                pos2 += 1
            elif overlaps(junctions2[pos2], junctions1[pos1]):
                features_present1[pos1] = -1
                features_present2[pos2] = -1
                if (junctions1[pos1][1] < junctions2[pos2][1]):
                    pos1 += 1
                else:
                    pos2 += 1
            elif left_of(junctions2[pos2], junctions1[pos1]):
                if pos1 > 0:
                    features_present2[pos2] = -1
                pos2 += 1
            else:
                if pos2 > 0:
                    features_present1[pos1] = -1
                pos1 += 1

        if any(el == -1 for el in features_present1) or any(el == -1 for el in features_present2):
            return "different_junctions"
        elif all(el == 0 for el in features_present1):
            print("Empty")
            print(junctions1)
            print(junctions2)
            print(features_present1)
            print(features_present2)
            if all(el == 0 for el in features_present2):
                return "both_empty"
            else:
                return "first_empty"
        elif all(el == 0 for el in features_present2):
            print("Empty")
            print(junctions1)
            print(junctions2)
            print(features_present1)
            print(features_present2)
            return "second_empty"
        elif features_present1[0] == 0 and features_present1[-1] == 0:
            if features_present2[0] == 0 and features_present2[-1] == 0:
                print("Clipperd tails")
                print(junctions1)
                print(junctions2)
                print(features_present1)
                print(features_present2)
            return "first_longer"
        elif features_present2[0] == 0 and features_present2[-1] == 0:
            return "second_longer"
        elif (features_present2[0] == 0 and features_present1[-1] == 0) or (features_present1[0] == 0 and features_present2[-1] == 0):
            return "overlap"
        else:
            print("Unknown")
            print(junctions1)
            print(junctions2)
            print(features_present1)
            print(features_present2)
        return "unknown"


    def compare_aligments(self, blocks1, blocks2):
        if len(blocks1) == 1:
            if len(blocks2) == 1:
                return "both_non_spliced"
            else:
                return "first_non_spliced"
        elif len(blocks2) == 1:
            return "second_non_spliced"

        introns1 = junctions_from_blocks(blocks1)
        introns2 = junctions_from_blocks(blocks2)
        return self.compare_junctions(introns1, introns2)


    def compare_alignment_sets(self, alignments1, alignments2):
        if len(alignments1) == 0:
            if len(alignments2) == 0:
                return "both_empty"
            else:
                return "first_empty"
        elif len(alignments2) == 0:
            return "second_empty"

        for a1 in alignments1:
            for a2 in alignments2:
                if a1.reference_id == a2.reference_id:
                    blocks1 = concat_gapless_blocks(sorted(a1.get_blocks()), a1.cigartuples)
                    blocks2 = concat_gapless_blocks(sorted(a2.get_blocks()), a2.cigartuples)
                    region1 = (blocks1[0][0], blocks1[-1][1])
                    region2 = (blocks2[0][0], blocks2[-1][1])

                    if overlaps(region1, region2):
                        res = self.compare_aligments(blocks1, blocks2)
                        if a1.is_secondary or a1.is_supplementary or a2.is_secondary or a2.is_supplementary:
                            return "secondary_" + res
                        else:
                            return res
                    else:
                        return "same_chr_diff_genes"
        return "different_chr"


    def compare_reads(self, read_pairs, alignment_map1, alignment_map2, barcode_map1, barcode_map2):
        self.stats = {}
        for read_pair in read_pairs:
            read_id1 = read_pair[0]
            read_id2 = read_pair[1]
    #        print(read_id1 + " " + read_id2)
            if barcode_map1[read_id1] != barcode_map2[read_id2]:
                print("Unequal barcode/UMI!")
                continue

            if read_id2 not in alignment_map2:
                if read_id1 not in alignment_map1:
                    self.stats[barcode_map1[read_id1]] = "none_in_bam"
                else:
                    self.stats[barcode_map1[read_id1]] = "first_only"
            elif read_id1 not in alignment_map1:
                self.stats[barcode_map1[read_id1]] = "second_only"
            else:
                self.stats[barcode_map1[read_id1]] = self.compare_alignment_sets(alignment_map1[read_id1], alignment_map2[read_id2])


    def print_stats(self, out_file):
        self.aggr_stats = {}
        outf = open(out_file, "w")
        for k in self.stats.keys():
            v = self.stats[k]
            if v not in self.aggr_stats:
                self.aggr_stats[v] = 0
                self.aggr_stats[v] += 1
            outf.write(k[0] + "\t" + k[1] + "\t" + self.stats[k] + "\n")

        outf.close()
        for k,v in self.aggr_stats:
            print(k + "\t" + str(v))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bams", nargs=2, help="bam files", type=str)
    parser.add_argument("--read_info",  help="file read ids with the same UMI and barcode but mapped to different genes ", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    parser.add_argument("--delta",  help="delta ", type=int, default=0)

    args = parser.parse_args()

    if args.bams is None or args.read_info is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()
    barcode_map1, barcode_map2, read_pairs = read_info(args.read_info)
    print("Read " + len(read_pairs) + " barcode/UMI pairs")
    print("Collecting infor from first BAM file")
    alignment_map1 = read_bam_file(args.bams[0], barcode_map1)
    print("Read " + len(alignment_map1) + " entries")
    print("Collecting infor from second BAM file")
    alignment_map2 = read_bam_file(args.bams[1], barcode_map2)
    print("Read " + len(alignment_map2) + " entries")

    print("Counting stats")
    alignment_comparator = AligmentComparator(args)
    alignment_comparator.compare_reads(read_pairs, alignment_map1, alignment_map2, barcode_map1, barcode_map2)
    alignment_comparator.print_stats(args.output_prefix)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
