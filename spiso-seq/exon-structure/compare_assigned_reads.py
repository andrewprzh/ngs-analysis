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
from common import *
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

    found_in1 = 0
    found_in2 = 0
    found_in_both = 0
    for barcode in inverted_read_info.keys():
        for umi in inverted_read_info[barcode].keys():
            if len(inverted_read_info[barcode][umi][0]) == 0:
                found_in2 += 1
                # del inverted_read_info[barcode][umi]
            elif len(inverted_read_info[barcode][umi][1]) == 0:
                found_in1 += 1
                # del inverted_read_info[barcode][umi]
            else:
                found_in_both += 1

    return inverted_read_info, found_in_both, found_in1, found_in2


class IsoformAssignmentStat:
    def __init__(self):
        self.equal_isoforms = 0
        self.unequal_isoforms = 0
        self.only_first_assigned = 0
        self.only_second_assigned = 0
        self.both_unassigned = 0

    def increment(self, isoform1, isoform2):
        if isoform2 == isoform1 == 'None':
            self.both_unassigned += 1
        elif isoform2 == isoform1:
            self.equal_isoforms += 1
        elif isoform1 == 'None':
            self.only_second_assigned += 1
        elif isoform2 == 'None':
            self.only_first_assigned += 1
        else:
            self.unequal_isoforms += 1

    def to_str(self):
        return str(self.equal_isoforms) + '\t' + str(self.unequal_isoforms) + '\t' + str(self.equal_isoforms) + \
               '\t' + str(self.only_first_assigned) + '\t' + str(self.only_second_assigned) + '\t' + str(self.both_unassigned)


class AssignedReadsComparator:
    def __init__(self, read_info_map1, read_info_map2, intersected_inverted_read_info, args):
        self.read_info_map1 = read_info_map1
        self.read_info_map2 = read_info_map2
        self.intersected_inverted_read_info = intersected_inverted_read_info
        self.args = args

        self.more_than_one_read = 0

        self.equal_profiles = IsoformAssignmentStat()
        self.first_intron_longer = IsoformAssignmentStat()
        self.second_intron_longer = IsoformAssignmentStat()
        self.different_introns = IsoformAssignmentStat()
        self.contradictory_introns = IsoformAssignmentStat()
        # Equal introns
        self.first_exons_longer = IsoformAssignmentStat()
        self.second_exons_longer = IsoformAssignmentStat()
        self.different_exons = IsoformAssignmentStat()
        self.contradictory_exons = IsoformAssignmentStat()

    def read_assigned_reads_section(self, in_file):
        assigned_reads_map = {}
        l = in_file.readline()
        while l and not l.startswith("ENSMUS"):
            # m64055_200112_012317/103025403/ccs      None    Empty   [0, 0, 0, 0, 0, 0, 0, 0]        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            tokens = l.strip().split('\t')
            read_id = tokens[0]
            isoform_id = tokens[1]
            assignment_type = tokens[2]
            intron_proflie = list(map(int, tokens[3][1:-1].split(',')))
            exon_profile = list(map(int, tokens[4][1:-1].split(',')))
            assigned_reads_map[read_id] = (isoform_id, assignment_type, intron_proflie, exon_profile)
            l = in_file.readline()

        return assigned_reads_map, l.split()[0] if l is not None else None

    def compare_assignments(self, assignement1, assignemtn2):
        (isoform_id1, assignment_type1, intron_proflie1, exon_profile1) = assignement1
        (isoform_id2, assignment_type2, intron_proflie2, exon_profile2) = assignement2
        if assignment_type1 == 'Empty' and assignment_type2 == 'Empty':
            return

        if intron_proflie1 == intron_proflie2:
            if exon_profile1 == exon_profile2:
                self.equal_profiles.increment(isoform_id1, isoform_id2)
            elif is_subprofile(exon_profile1, exon_profile2):
                self.second_exons_longer.increment(isoform_id1, isoform_id2)
            elif  is_subprofile(exon_profile2, exon_profile1):
                self.first_exons_longer.increment(isoform_id1, isoform_id2)
            elif diff_only_present(exon_profile1, exon_profile2) == 0:
                self.different_exons.increment(isoform_id1, isoform_id2)
            else:
                self.contradictory_exons.increment(isoform_id1, isoform_id2)
        elif is_subprofile(intron_proflie1, intron_proflie2):
            self.second_intron_longer.increment(isoform_id1, isoform_id2)
        elif is_subprofile(intron_proflie2, intron_proflie1):
            self.first_intron_longer.increment(isoform_id1, isoform_id2)
        elif diff_only_present(intron_proflie1, intron_proflie2) == 0:
            self.different_introns.increment(isoform_id1, isoform_id2)
        else:
            self.contradictory_introns.increment(isoform_id1, isoform_id2)

    def compare_two_sets(self, assigned_reads_map1, assigned_reads_map2):
        for read_id1 in assigned_reads_map1:
            if read_id1 not in self.read_info_map1:
                continue
            (barcode, umi) = self.read_info_map1[read_id1]
            if len(self.intersected_inverted_read_info[barcode][umi][1]) == 0:
                continue

            if len(self.intersected_inverted_read_info[barcode][umi][0]) > 1 and len(self.intersected_inverted_read_info[barcode][umi][1]) > 1:
                self.more_than_one_read += 1
                continue

            read_id2 = self.intersected_inverted_read_info[barcode][umi][1][0]
            self.compare_assignments(assigned_reads_map1[read_id1], assigned_reads_map2[read_id2])

    def process(self):
        infile1 = open(self.args.assigned_reads_files[0], "r")
        infile2 = open(self.args.assigned_reads_files[1], "r")

        geneid1 = ""
        geneid2 = ""
        while geneid1 is not None and geneid2 is not None:
            print("Processing " + geneid1)
            assigned_reads_map1, geneid1 = self.read_assigned_reads_section(infile1)
            assigned_reads_map2, geneid2 = self.read_assigned_reads_section(infile2)

            self.compare_two_sets(assigned_reads_map1, assigned_reads_map2)

    def print_stat(self):
        print("Barcode-UMI pairs with > 1 read " +str(self.more_than_one_read))
        print("Profile comparison\tequal\tdiff\tassign1\tassign2\tunassigned")
        print("equal_profiles\t" + self.equal_profiles.to_str())
        print("first_intron_longer\t" + self.first_intron_longer.to_str())
        print("second_intron_longer\t" + self.second_intron_longer.to_str())
        print("different_introns\t" + self.different_introns.to_str())
        print("contradictory_introns\t" + self.contradictory_introns.to_str())
        print("first_exons_longer\t" + self.first_exons_longer.to_str())
        print("second_exons_longer\t" + self.second_exons_longer.to_str())
        print("different_exons\t" + self.different_exons.to_str())
        print("contradictory_exons\t" + self.contradictory_exons.to_str())


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assigned_reads_files", nargs=2, help="files with assigned reads", type=str)
    parser.add_argument("--read_info", nargs=2, help="file with UMIs and barcodes ", type=str)
    parser.add_argument("--barcode_column", help="barcode column", type=int, default=4)
    parser.add_argument("--umi_column", help="umi column", type=int, default=7)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    args = parser.parse_args()

    if args.assigned_reads_files is None or args.read_info is None:
        parser.print_help()
        exit(-1)

    return args


def main():
    args = parse_args()

    read_info_map1 = get_read_info(args.read_info[0], args)
    print("Read " + str(len(read_info_map1)) + " reads from " + args.read_info[0])
    read_info_map2 = get_read_info(args.read_info[1], args)
    print("Read " + str(len(read_info_map2)) + " reads from " + args.read_info[1])

    inverted_read_info, found_in_both, found_in1, found_in2 = intersect_read_infos(read_info_map1, read_info_map2)
    print("Barcode-UMI pairs:")
    print("Found in both " + str(found_in_both))
    print("Found in 1 " + str(found_in1))
    print("Found in 2 " + str(found_in2))

    print("Comparing equal Barcode-UMI pairs")
    comparator = AssignedReadsComparator(read_info_map1, read_info_map2, inverted_read_info, args)
    comparator.process()
    comparator.print_stat()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
