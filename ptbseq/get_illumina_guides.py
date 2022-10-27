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
    parser.add_argument("--guide_start", type=int, help="guide starting position in the reference, 1-based", default=1610)
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
        self.min_count=2
        self.check_alignment = True
        self.guide_region = (args.guide_start - 1, args.guide_start - 2 + args.guide_len)

    def load_guides(self):
        guide_set = set()
        full_ref_sequences = {}
        for r in SeqIO.parse(self.args.guide_ref, "fasta"):
            guide_seq = str(r.seq)[self.guide_region[0]:self.guide_region[1] + 1]
            guide_set.add(guide_seq)
            full_ref_sequences[r.id] = str(r.seq)
        return full_ref_sequences, guide_set

    def verify_alignment(self, alignment, reference_sequence):
        # ref_region, 0-based closed interval
        if intersection_len((alignment.reference_start, alignment.reference_end), self.guide_region) < self.minimal_overlap:
            return False

        read_sequence = alignment.query_sequence
        aligned_pairs = alignment.get_aligned_pairs()
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

        return subs_count <= self.max_indels and indel_counts <= self.max_mismatches

    def get_guide_barcodes(self, in_bam, guide_id, reference_sequence):
        barcode_umis = defaultdict(set)
        barcode_counts = defaultdict(int)
        for a in in_bam.fetch(guide_id):
            if a.is_supplementary or a.is_secondary:
                continue

            if self.check_alignment and not self.verify_alignment(a, reference_sequence):
                continue
            try:
                barcode = a.get_tag("CR")
                umi = a.get_tag("UR")
                barcode_umis[barcode].add(umi)
            except KeyError:
                continue
            #elif guide_seq in guide_set:
            #    barcode_counts[guide_seq] += 1

        #print(barcode_umis)
        for bc in barcode_umis.keys():
            barcode_counts[bc] += len(barcode_umis[bc])
        return barcode_counts

    def filter_barcodes(self, barcode_to_guide):
        filtered_barcode_to_guide = defaultdict(set)
        for bc in barcode_to_guide.keys():
            if len(barcode_to_guide[bc]) <= 1:
                filtered_barcode_to_guide[bc].update(barcode_to_guide[bc])
                continue
            top_count = max(x[1] for x in barcode_to_guide[bc])
            cutoff = 1 if top_count == 1 else max(self.min_count, top_count / 2.0)
            for guide_pair in barcode_to_guide[bc]:
                if guide_pair[1] >= cutoff:
                    filtered_barcode_to_guide[bc].add(guide_pair)
        return filtered_barcode_to_guide

    def reverse_map(self, barcode_to_guide):
        guide_to_barcode = defaultdict(set)
        for bc in barcode_to_guide.keys():
            for guide_pair in barcode_to_guide[bc]:
                guide_to_barcode[guide_pair[0]].add((bc, guide_pair[1]))
        return guide_to_barcode

    def dict_stats(self, dict_with_counts, key_label="barcodes", value_label="guides"):
        print("Total %s: %d" % (key_label, len(dict_with_counts)))
        counts = defaultdict(int)
        for k in dict_with_counts.keys():
            counts[len(dict_with_counts[k])] += 1
        print("# distinct %s\t# %s" % (value_label, key_label))
        for k in sorted(counts.keys()):
            print("%d\t%d" % (k, counts[k]))

    def collect_barcodes(self):
        barcode_to_guide = defaultdict(set)
        full_ref_sequences, guide_set = self.load_guides()

        in_bam = pysam.AlignmentFile(self.args.bam, "rb")

        for ref_id in full_ref_sequences.keys():
            barcode_counts = self.get_guide_barcodes(in_bam, ref_id, full_ref_sequences[ref_id])
            for bc in barcode_counts.keys():
                count = barcode_counts[bc]
                barcode_to_guide[bc].add((ref_id, count))

        print("Unfiltered barcodes: %d" % len(barcode_to_guide))
        self.dict_stats(barcode_to_guide)
        filtered_barcode_to_guide = self.filter_barcodes(barcode_to_guide)
        print("Filtered barcodes: %d" % len(filtered_barcode_to_guide))
        self.dict_stats(filtered_barcode_to_guide)

        guide_to_barcode = self.reverse_map(barcode_to_guide)
        print("Unfiltered guides: %d" % len(guide_to_barcode))
        self.dict_stats(guide_to_barcode)
        filtered_guide_to_barcode = self.reverse_map(filtered_barcode_to_guide)
        print("Filtered guides: %d" % len(filtered_guide_to_barcode))
        self.dict_stats(filtered_guide_to_barcode)


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



