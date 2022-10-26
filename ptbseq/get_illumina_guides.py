#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
from Bio import SeqIO
import logging

logger = logging.getLogger('GuideCaller')


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", default="guides_output")
    parser.add_argument("--guide_ref", "-g", type=str, help="guides sequences in FASTA format", required=True)
    parser.add_argument("--guide_start", type=int, help="guide starting position in the reference, 1-based", default=1609)
    parser.add_argument("--guide_len", type=int, help="guide length", default=20)
    parser.add_argument("--bam", "-b", type=str, help="BAM file with Illumina reads", required=True)
    args = parser.parse_args()
    return args


# == range operations ==
def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def overlap_intervals(range1, range2):
    return max(range1[0], range2[0]), min(range1[1], range2[1])


def overlaps_at_least(range1, range2, delta=0):
    cutoff = min([delta, range1[1] - range1[0] + 1, range2[1] - range2[0] + 1])
    overlap = min(range1[1], range2[1]) + 1 - max(range1[0], range2[0])
    return overlap >= cutoff


def intersection_len(range1, range2):
    return max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)


class GuideCaller:
    def __init__(self, args):
        self.args = args
        self.minimal_overlap=10
        self.max_mismatches=0
        self.max_indels=0
        self.min_count=1
        self.guide_region = (args.guide_start - 1, args.guide_start - 2 + args.guide_len)

    def load_guides(self):
        guide_set = set()
        full_ref_sequences = {}
        for r in SeqIO.parse(self.args.guide_ref, "fasta"):
            guide_seq = str(r.seq)[self.guide_region[0]:self.guide_region[1] + 1]
            guide_set.add(guide_seq)
            full_ref_sequences[r.id] = str(r.seq)
        return full_ref_sequences, guide_set

    def get_guide_barcodes(self, in_bam, guide_id, reference_sequence):
        barcode_umis = defaultdict(set)
        barcode_counts = defaultdict(int)
        for a in in_bam.fetch(guide_id):
            if a.is_supplementary or a.is_secondary:
                continue
            # ref_region, 0-based closed interval
            if intersection_len((a.reference_start, a.reference_end), self.guide_region) < self.minimal_overlap:
                continue

            read_sequence = a.query_sequence
            aligned_pairs = a.get_aligned_pairs()
            start_index = -1
            end_index = -1
            for i, ap in enumerate(aligned_pairs):
                if ap[1] is None:
                    continue
                elif start_index == -1 and ap[1] >= self.guide_region[0]:
                    start_index = i
                elif ap[1] <= self.guide_region[1]:
                    end_index = i

            subs_count = 0
            indel_counts = 0
            read_guide_seq = ""
            for i in range(start_index, end_index + 1):
                ap = aligned_pairs[i]
                if ap[0] is None:
                    indel_counts += 1
                    continue

                read_nucl = read_sequence[ap[0]]
                read_guide_seq += read_nucl
                if ap[1] is None:
                    indel_counts += 1
                    continue

                ref_nucl = reference_sequence[ap[1]]
                if read_nucl != ref_nucl:
                    subs_count += 1

            if subs_count <= self.max_indels and indel_counts <= self.max_mismatches:
                try:
                    barcode = a.get_tag("CR")
                    umi = a.get_tag("CR")
                    barcode_umis[barcode].add(umi)
                except KeyError:
                    continue
            #elif guide_seq in guide_set:
            #    barcode_counts[guide_seq] += 1

        for bc in barcode_umis.keys():
            barcode_counts[bc] += len(barcode_umis[bc])
        return barcode_counts

    def collect_barcodes(self):
        barcode_to_guide = defaultdict(set)
        filtered_barcode_to_guide = defaultdict(set)
        guide_to_barcode = defaultdict(set)
        filtered_guide_to_barcode = defaultdict(set)
        full_ref_sequences, guide_set = self.load_guides()
        in_bam = pysam.AlignmentFile(self.args.bam, "rb")

        for ref_id in full_ref_sequences.keys():
            barcode_counts = self.get_guide_barcodes(in_bam, ref_id, full_ref_sequences[ref_id])
            print(barcode_counts)
            for bc in barcode_counts.keys():

                count = barcode_counts[bc]
                barcode_to_guide[bc].add((ref_id, count))
                guide_to_barcode[ref_id].add((bc, count))
                if count > self.min_count:
                    filtered_barcode_to_guide[bc].add((ref_id, count))
                    filtered_guide_to_barcode[ref_id].add((bc, count))

        print("Unfiltered barcodes: %d" % len(barcode_to_guide))
        for bc in barcode_to_guide:
            print("%s\t%s" % (bc, ", ".join(["%s,%d" % (x[0], x[1]) for x in barcode_to_guide[bc]])))

        print("Filtered barcodes: %d" % len(filtered_barcode_to_guide))
        for bc in filtered_barcode_to_guide:
            print("%s\t%s" % (bc, ", ".join(["%s,%d" % (x[0], x[1]) for x in filtered_barcode_to_guide[bc]])))

        print("Unfiltered guides: %d" % len(guide_to_barcode))
        for guide in guide_to_barcode:
            print("%s\t%s" % (guide, ", ".join(["%s,%d" % (x[0], x[1]) for x in guide_to_barcode[guide]])))
        print("Filtered guides: %d" % len(filtered_guide_to_barcode))
        for guide in filtered_guide_to_barcode:
            print("%s\t%s" % (guide, ", ".join(["%s,%d" % (x[0], x[1]) for x in filtered_guide_to_barcode[guide]])))


def main():
    #set_logger(logger)
    args = parse_args()
    guide_caller = GuideCaller(args)
    guide_caller.collect_barcodes()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



