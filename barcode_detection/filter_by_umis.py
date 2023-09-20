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
        self.args = args
        self.max_edit_distance = args.min_distance
        self.disregard_length_diff = args.disregard_length_diff
        self.barcode_dict = self.load_barcodes(args.barcodes)
        self.discarded = 0
        self.no_barcode = 0
        self.no_gene = 0
        self.total_reads = 0
        self.ambiguous = 0
        self.monoexonic = 0
        self.duplicated_molecule_counts = defaultdict(int)

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
            if not self.args.untrusted_umis and v[9] != "True":
                continue
            umi = v[7]
            barcode_dict[v[0]] = (barcode, umi)
        logger.info("Loaded %d barcodes " % len(barcode_dict))
        return barcode_dict

    def _find_similar_umi(self, umi, umi_indexer):
        if self.max_edit_distance == -1:
            return None if umi_indexer.empty() else umi_indexer.seq_list[0]
        occurrences = umi_indexer.seq_list
        similar_umi = None
        best_dist = 100
        # choose the best UMI among checked
        for occ in occurrences:
            if self.max_edit_distance == 0:
                if self.disregard_length_diff:
                    similar, ed = occ == umi, 0
                elif len(occ) < len(umi):
                    similar, ed = occ in umi, abs(len(occ) - len(umi))
                else:
                    similar, ed =  umi in occ, abs(len(occ) - len(umi))
            else:
                ed = editdistance.eval(occ, umi)
                if not self.disregard_length_diff:
                    ed -= abs(len(occ) - len(umi))
                similar, ed = ed <= self.max_edit_distance, ed

            if similar and ed < best_dist:
                similar_umi = occ
                best_dist = ed

        return similar_umi

    def _process_duplicates(self, molecule_list):
        if len(molecule_list) <= 1:
            self.duplicated_molecule_counts[1] += 1
            logger.debug("Unique " + molecule_list[0].read_id)
            return [x.read_id for x in molecule_list]

        resulting_reads = []
        umi_dict = defaultdict(list)
        umi_indexer = KmerIndexer([], kmer_size=3)
        for m in molecule_list:
            similar_umi = self._find_similar_umi(m.umi, umi_indexer)
            if similar_umi is None:
                # if none is added, add to indexer
                umi_dict[m.umi].append(m)
                umi_indexer.append(m.umi)
            else:
                umi_dict[similar_umi].append(m)

        if len(umi_dict) > 1 and random.random() < 0.001:
            for umi in umi_dict.keys():
                logger.info(umi)
                for m in umi_dict[umi]:
                    umi2 = m.umi
                    ed = editdistance.eval(umi, umi2)
                    if not self.disregard_length_diff:
                        ed -= abs(len(umi) - len(umi2))
                    logger.info("  %s: %d" % (umi2, ed))

        for umi in umi_dict.keys():
            duplicate_count = len(umi_dict[umi])
            self.duplicated_molecule_counts[duplicate_count] += 1
            if duplicate_count == 1:
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
            logger.info("+ Checking barcode %s +" % barcode)
            resulting_reads += self._process_duplicates(gene_dict[barcode])
        return resulting_reads

    def _process_chunk(self, gene_barcode_dict, outf):
        read_count = 0
        for gene_id in gene_barcode_dict:
            logger.info("== Processing gene %s ==" % gene_id)
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
            assignment_type = v[5]
            self.total_reads += 1
            if gene_id == ".":
                self.no_gene += 1
                continue
            if read_id not in self.barcode_dict:
                self.no_barcode += 1
                continue
            if assignment_type == "ambiguous":
                self.ambiguous += 1
                if self.args.only_unique: continue
            if read_id in read_to_gene and read_to_gene[read_id] != gene_id:
                if assignment_type != "ambiguous":
                    self.ambiguous += 1
                # ignore ambiguous gene assignments
                continue
            exon_blocks_str = v[7]
            exon_blocks = list(map(lambda x: tuple(map(int, x.split('-'))), exon_blocks_str.split(',')))
            if len(exon_blocks) == 1:
                self.monoexonic += 1
                if self.args.only_spliced: continue

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
        logger.info("Total assignments processed %d" % self.total_reads)
        logger.info("Unspliced %d %s" % (self.monoexonic, "(discarded)" if self.args.only_spliced else ""))
        logger.info("Ambiguous %d %s" % (self.ambiguous, "(discarded)" if self.args.only_unique else ""))
        logger.info("No barcode detected %d" % self.no_barcode)
        logger.info("No gene assigned %d" % self.no_gene)
        logger.info("Discarded as duplicates %d" % self.discarded)
        logger.info("Duplicate count histogram:")

        outf = self.args.output + ".counts"
        with open(outf, "w") as count_hist_file:
            for k in sorted(self.duplicated_molecule_counts.keys()):
                count_hist_file.write("%d\t%d\n" % (k, self.duplicated_molecule_counts[k]))


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
    parser.add_argument("--min_distance", type=int, help="minimal edit distance for UMIs to be considered distinct;"
                                                         "length difference is added to this values by default;"
                                                         "0 for equal UMIs, -1 for keeping only a single gene-barcode "
                                                         "read", default=2)
    parser.add_argument("--untrusted_umis", action="store_true", help="keep only trusted UMIs", default=False)
    parser.add_argument("--only_spliced", action="store_true", help="keep only spliced reads", default=False)
    parser.add_argument("--only_unique", action="store_true", help="keep only non-ambiguous reads", default=False)
    parser.add_argument("--disregard_length_diff", action="store_true", help="do not account for length difference "
                                                                             "when comapring UMIs", default=False)


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
