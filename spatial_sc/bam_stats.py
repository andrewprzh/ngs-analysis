#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc
import os
import sys
import pysam
from collections import defaultdict
from enum import Enum
import gffutils


def read_called_barcodes(inf):
    total_reads = 0
    barcode_dict = {}
    for l in open(inf):
        total_reads += 1
        v = l.strip().split('\t')
        if len(v) >= 2 and v[1] == "*":
            continue
        barcode_dict[v[0]] = v[1]
        if total_reads % 100000 == 0:
            sys.stdout.write("Processed " + str(total_reads) + " reads\r")

    print("Loaded %d read ids" % len(barcode_dict))
    return barcode_dict


def load_intron_chr_dict(genedb_path):
    gene_db = gffutils.FeatureDB(genedb_path)
    intron_chr_dict = defaultdict(set)
    total_genes = 0

    for g in gene_db.features_of_type('gene'):
        total_genes += 1
        if total_genes % 100000 == 0:
            sys.stdout.write("Loaded " + str(total_genes) + " genes\r")
        for t in gene_db.children(g, featuretype=('transcript', 'mRNA')):
            exons = []
            for e in gene_db.children(t, featuretype='exon', order_by='start'):
                exons.append((e.start, e.end))

            for intron in junctions_from_blocks(exons):
                intron_chr_dict[g.seqid].add(intron)

    return intron_chr_dict


def has_intron(cigartules):
    return any(x[0] == 3 for x in cigartules)


class CigarEvent(Enum):
    match = 0
    insertion = 1
    deletion = 2
    skipped = 3
    soft_clipping = 4
    hard_clipping = 5
    padding = 6
    seq_match = 7
    seq_mismatch = 8

    @classmethod
    def get_match_events(cls):
        return {cls.match, cls.seq_match, cls.seq_mismatch}

    @classmethod
    def get_ins_del_match_events(cls):
        return {cls.match, cls.insertion, cls.deletion, cls.seq_match, cls.seq_mismatch}


def get_read_blocks(ref_start, cigar_tuples):
    read_pos = 0
    ref_pos = ref_start + 1
    cigar_index = 0
    current_ref_block_start = None
    current_read_block_start = None
    current_cigar_block_start = None
    has_match = False
    ref_blocks = []
    cigar_blocks = []
    read_blocks = []

    while cigar_index < len(cigar_tuples):
        cigar_event = CigarEvent(cigar_tuples[cigar_index][0])
        event_len = cigar_tuples[cigar_index][1]

        if current_ref_block_start is None and cigar_event in CigarEvent.get_ins_del_match_events():
            # init new block from match
            current_ref_block_start = ref_pos
            current_read_block_start = read_pos
            current_cigar_block_start = cigar_index
            if cigar_event == CigarEvent.insertion:
                read_pos += event_len
            elif cigar_event == CigarEvent.deletion:
                ref_pos += event_len
            else:
                read_pos += event_len
                ref_pos += event_len
                has_match = True
        # found intron, add current block
        elif cigar_event in CigarEvent.get_match_events():
            read_pos += event_len
            ref_pos += event_len
            has_match = True
        elif cigar_event == CigarEvent.insertion:
            read_pos += event_len
        elif cigar_event == CigarEvent.deletion:
            ref_pos += event_len
        elif cigar_event == CigarEvent.skipped:
            if current_ref_block_start:
                if has_match:
                    ref_blocks.append((current_ref_block_start, ref_pos - 1))
                    read_blocks.append((current_read_block_start, read_pos - 1))
                    cigar_blocks.append((current_cigar_block_start, cigar_index - 1))
                has_match = False
                current_ref_block_start = None
                current_read_block_start = None
                current_cigar_block_start = None
            ref_pos += event_len
        elif cigar_event == CigarEvent.soft_clipping:
            if current_ref_block_start:
                if has_match:
                    ref_blocks.append((current_ref_block_start, ref_pos - 1))
                    read_blocks.append((current_read_block_start, read_pos - 1))
                    cigar_blocks.append((current_cigar_block_start, cigar_index - 1))
                has_match = False
                current_ref_block_start = None
                current_read_block_start = None
                current_cigar_block_start = None
            read_pos += event_len

        cigar_index += 1

    if current_ref_block_start and has_match:
        ref_blocks.append((current_ref_block_start, ref_pos - 1))
        read_blocks.append((current_read_block_start, read_pos - 1))
        cigar_blocks.append((current_cigar_block_start, cigar_index - 1))

    return ref_blocks, read_blocks, cigar_blocks


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def correct_bam_coords(blocks):
    return list(map(lambda x: (x[0] + 1, x[1]), blocks))


def get_introns(read):
    exons, _, _ = get_read_blocks(read.reference_start, read.cigartuples)
    return junctions_from_blocks(exons)


def count_known(introns, chr_id, intron_chr_dict):
    known_introns = 0
    if intron_chr_dict:
        for i in introns:
            known_introns += 1 if i in intron_chr_dict[chr_id] else 0

    return known_introns


def count_stats(in_file_name, barcode_dict=None, intron_chr_dict=None):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    total_reads = 0
    stat_dict = defaultdict(int)
    intron_barcode_dict = defaultdict(int)
    known_intron_barcode_dict  = defaultdict(int)

    for read in inf:
        total_reads += 1

        if read.reference_id == -1:
            mapped = "unmapped"
        elif not read.is_secondary and not read.is_supplementary:
            mapped = "primary"
        else:
            mapped = "non-primary"

        is_barcoded = barcode_dict is not None and read.query_name in barcode_dict
        if is_barcoded:
            barcoded = "barcoded"
        else:
            barcoded = "no barcode"

        introns = get_introns(read)
        is_spliced = len(introns) > 0
        chr_id = read.reference_name
        known_introns = count_known(introns, chr_id, intron_chr_dict)
        all_introns = len(introns)
        if is_spliced:
            if known_introns == all_introns:
                spliced = "spliced_known"
            elif known_introns == 0:
                spliced = "spliced_all_novel"
            else:
                spliced = "spliced_novel"
        else:
            spliced = "unspliced"

        if is_barcoded and is_spliced:
            barcode = barcode_dict[read.query_name]
            for i in introns:
                intron_barcode_dict[(barcode, chr_id, i)] += 1
                if i in intron_chr_dict[chr_id]:
                    known_intron_barcode_dict[(barcode, chr_id, i)] += 1

        stat_dict[(mapped, barcoded, spliced)] += 1

        if total_reads % 100000 == 0:
            sys.stdout.write("Processed " + str(total_reads) + " reads\r")

    print("Processed %d reads" % total_reads)
    inf.close()
    for k in sorted(stat_dict.keys()):
        print("%s\t%d" % (k, stat_dict[k]))

    print("Barcode-intron pairs: %d" % len(intron_barcode_dict))
    print("Barcode-known-intron pairs: %d" % len(known_intron_barcode_dict))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("--output", "-o", type=str, help="output folder (same as input file if not set), "
    #                                                      "or file name")
    parser.add_argument("--barcodes", type=str, help="read barcodes")
    parser.add_argument("--genedb", type=str, help="gffutils genedb (for intron stats)")
    parser.add_argument("--bam", "-b", nargs="+", type=str, help="BAM file ", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    barcode_dict = None
    if args.barcodes:
        print("Loading barcodes from %s" % args.barcodes)
        barcode_dict = read_called_barcodes(args.barcodes)
        print("Loaded %d barcodes" % len(barcode_dict))

    intron_chr_dict = None
    if args.genedb:
        print("Loading introns from %s" % args.genedb)
        intron_chr_dict = load_intron_chr_dict(args.genedb)
        print("Loaded %d introns" % sum([len(x) for x in intron_chr_dict.values()]))

    for in_bam in args.bam:
        print("Processing %s" % in_bam)
        count_stats(in_bam, barcode_dict, intron_chr_dict)
        print("End processing %s" % in_bam)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
