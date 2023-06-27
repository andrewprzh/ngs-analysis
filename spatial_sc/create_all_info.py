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
from traceback import print_exc
from collections import defaultdict
import logging


logger = logging.getLogger('AllInfo')


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


class AllInfoGenerator:
    def __init__(self, barcodes):
        self.barcode_dict = self.load_barcodes(barcodes)
        self.good_barcode_score = 13

    def load_barcodes(self, in_file):
        # afaf0413-4dd7-491a-928b-39da40d68fb3    99      56      66      83      +       GCCGATACGCCAAT  CGACTGAAG       13      True
        logger.info("Loading barcodes from " + in_file)
        barcode_dict = {}
        for l in open(in_file):
            v = l.strip().split("\t")
            if len(v) < 10: continue
            barcode = v[6]
            if barcode == "*":
                continue
            umi = v[7]
            score = int(v[8])
            trusted_umi = v[9]
            barcode_dict[v[0]] = (barcode, umi, score, trusted_umi)
        logger.info("Loaded %d barcodes " % len(barcode_dict))
        return barcode_dict

    def convert_to_allinfo(self, assignment_file, output_file):
        outf = open(output_file, "w")
        logger.info("Converting " + assignment_file)
        processed_reads = set()
        assignment_count = 0
        # 1251f521-d4c2-4dcc-96d7-85070cc44e12    chr1    -       ENST00000477740.5       ENSG00000238009.6       ambiguous       ism_internal    112699-112804,120721-120759     PolyA=False
        for l in open(assignment_file):
            if l.startswith("#"): continue
            v = l.strip().split("\t")
            read_id = v[0]
            chr_id = v[1]
            strand = v[2]
            gene_id = v[4]
            assignment_type = v[5]
            if not assignment_type.startswith("unique"):
                if read_id in processed_reads:
                    continue
                processed_reads.add(read_id)
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))
            intron_blocks = junctions_from_blocks(exon_blocks)
            exons_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (chr_id, e[0], e[1], strand) for e in exon_blocks])
            introns_str = ";%;" + ";%;".join(["%s_%d_%d_%s" % (chr_id, e[0], e[1], strand) for e in intron_blocks])

            barcode = "*"
            umi = "*"
            if read_id in self.barcode_dict:
                score = self.barcode_dict[read_id][2]
                if score >= self.good_barcode_score:
                    barcode = self.barcode_dict[read_id][0]
                    umi = self.barcode_dict[read_id][1]
            cell_type = "None"
            read_type = "known" if assignment_type.startswith("unique") else "novel"

            polyA = "NoPolyA"
            TSS = "NoTSS"
            matching_events = v[6]
            if "tss_match" in matching_events:
                tss_pos = exon_blocks[-1][1] if strand == "-" else exon_blocks[0][0]
                TSS = "%s_%d_%d_%s" % (chr_id, tss_pos, tss_pos, strand)
            if "correct_polya" in matching_events:
                polyA_pos = exon_blocks[-1][1] if strand == "+" else exon_blocks[0][0]
                polyA = "%s_%d_%d_%s" % (chr_id, polyA_pos, polyA_pos, strand)
            outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % (read_id, gene_id, cell_type, barcode, umi,
                                                                         introns_str, TSS, polyA, exons_str, read_type,
                                                                         len(intron_blocks)))
            assignment_count += 1

        logger.info("Converted %d assignments " % assignment_count)
        outf.close()
        return assignment_count


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", type=str, help="read - barcode - UMI table", required=True)
    parser.add_argument("--read_assignments", "-r", type=str, help="IsoQuant read assignments", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)
    convertor = AllInfoGenerator(args.barcodes)
    # convertor.good_barcode_score = 13
    # convertor.convert_to_allinfo(args.read_assignmendt, args.output + ".score13.tsv")
    convertor.good_barcode_score = 14
    convertor.convert_to_allinfo(args.read_assignments, args.output + ".score14.tsv")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
