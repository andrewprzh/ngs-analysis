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
import editdistance

from kmer_indexer import KmerIndexer


logger = logging.getLogger('UMIFilter')


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


class MoleculeInfo:
    def __init__(self, read_id, umi, exon_blocks):
        self.read_id = read_id
        self.umi = umi
        self.exon_blocks = exon_blocks


class UMIFilter:
    def __init__(self, args):
        self.max_edit_distance = args.min_distance
        self.barcode_dict = self.load_barcodes(args.barcodes, not args.untrusted_umis)
        self.discarded = 0


    def load_barcodes(self, in_file, trusted_umis=False):
        # afaf0413-4dd7-491a-928b-39da40d68fb3    99      56      66      83      +       GCCGATACGCCAAT  CGACTGAAG       13      True
        logger.info("Loading barcodes from " + in_file)
        barcode_dict = {}
        for l in open(in_file):
            v = l.strip().split("\t")
            if len(v) < 10: continue
            barcode = v[6]
            if barcode == "*":
                continue
            if trusted_umis and v[9] != "True":
                continue
            umi = v[7]
            barcode_dict[v[0]] = (barcode, umi)
        logger.info("Loaded %d barcodes " % len(barcode_dict))
        return barcode_dict

    def _process_duplicates(self, molecule_list):
        if len(molecule_list) <= 1:
            logger.debug("Unique " + molecule_list[0].read_id)
            return [x.read_id for x in molecule_list]

        resulting_reads = []
        umi_dict = defaultdict(list)
        umi_indexer = KmerIndexer([], kmer_size=3)
        for m in molecule_list:
            occurrences = umi_indexer.get_occurrences(m.umi, hits_delta=1)
            similar_umi = None
            best_dist = 100
            # choose the best UMI among checked
            for occ in occurrences.keys():
                ed = editdistance.eval(occ, m.umi)
                if ed <= self.max_edit_distance and ed < best_dist:
                    similar_umi = occ
                    best_dist = ed

            if similar_umi is None:
                # in none is added, add to indexer
                umi_dict[m.umi].append(m)
                umi_indexer.append(m.umi)
            else:
                umi_dict[similar_umi].append(m)

        for umi in umi_dict.keys():
            if len(umi_dict[umi]) == 1:
                resulting_reads.append(umi_dict[umi][0].read_id)
                continue

            best_read = umi_dict[umi][0]
            logger.debug("Selecting from:")
            for m in umi_dict[umi]:
                logger.debug("%s %s" % (m.read_id, m.umi))
                if len(m.exon_blocks) > len(best_read.exon_blocks):
                    best_read = m
                elif len(m.exon_blocks) == len(best_read.exon_blocks) and \
                        m.exon_blocks[-1][1] - m.exon_blocks[0][0] > \
                        best_read.exon_blocks[-1][1] - best_read.exon_blocks[0][0]:
                    best_read = m
            logger.debug("Selected %s %s" % (best_read.read_id, best_read.umi))
            self.discarded += len(umi_dict[umi]) - 1

            resulting_reads.append(best_read.read_id)

        return resulting_reads

    def _process_gene(self, gene_dict):
        resulting_reads = []
        for barcode in gene_dict:
            resulting_reads += self._process_duplicates(gene_dict[barcode])
        return resulting_reads

    def _process_chunk(self, gene_barcode_dict, outf):
        read_count = 0
        for gene_id in gene_barcode_dict:
            read_list = self._process_gene(gene_barcode_dict[gene_id])
            read_count += len(read_list)
            for read_id in read_list:
                outf.write(read_id + "\n")
        return read_count

    def process(self, assignment_file, output_file):
        outf = open(output_file, "w")
        # 1251f521-d4c2-4dcc-96d7-85070cc44e12    chr1    -       ENST00000477740.5       ENSG00000238009.6       ambiguous       ism_internal    112699-112804,120721-120759     PolyA=False
        gene_barcode_dict = defaultdict(lambda: defaultdict(list))
        read_to_gene = {}
        current_interval = (-1, -1)
        current_chr = None
        read_count = 0
        for l in open(assignment_file):
            if l.startswith("#"): continue
            v = l.strip().split("\t")
            read_id = v[0]
            chr_id = v[1]
            gene_id = v[4]
            if gene_id == "." or read_id not in self.barcode_dict:
                continue
            if read_id in read_to_gene and read_to_gene[read_id] != gene_id:
                # ignore ambiguous gene assignments
                continue
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))

            read_interval = (exon_blocks[0][0], exon_blocks[-1][1])
            if current_chr != chr_id or not overlaps(current_interval, read_interval):
                read_count += self._process_chunk(gene_barcode_dict, outf)
                if current_chr != chr_id:
                    logger.info("Processing chromosome " + chr_id)
                current_chr = chr_id
                current_interval = read_interval
                gene_barcode_dict.clear()
                read_to_gene.clear()

            barcode, umi = self.barcode_dict[read_id]
            gene_barcode_dict[gene_id][barcode].append(MoleculeInfo(read_id, umi, exon_blocks))
            read_to_gene[read_id] = gene_id
            current_interval = (current_interval[0], max(current_interval[1], read_interval[1]))

        read_count += self._process_chunk(gene_barcode_dict, outf)
        outf.close()
        logger.info("Saved %d reads to %s" % (read_count, output_file))
        logger.info("Discarded %d" % self.discarded)


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
    parser.add_argument("--min_distance", type=int, help="minimal edit distance for UMIs to be considered distinct", default=3)
    parser.add_argument("--untrusted_umis", action="store_true", help="keep only trusted UMIs", default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    umi_filter = UMIFilter(args)
    umi_filter.process(args.read_assignments, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
