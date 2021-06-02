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
import gffutils
from enum import Enum
import logging

logger = logging.getLogger('IsoQuantQA')

def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def get_original_isoform(read_id):
    return read_id.split('_')[1]


class TrackingData:
    def __init__(self, tracking_file):
        logger.info("Reading tracking data from %s" % tracking_file)
        self.all_ref_isoforms = set()
        self.unmapped_isoforms = set()
        self.duplicate_isoforms = {}
        self.incorrect_isoforms = []
        all_ref_tracking = defaultdict(list)
        for l in open(tracking_file):
            t = l.strip().split()
            assignment_type = t[3]
            isoform_id = t[4].split(':', 2)[1].split('|')[1]
            if assignment_type == 'u':
                self.unmapped_isoforms.add(isoform_id)
                continue
            ref_isoform = t[2].split('|')[1]
            all_ref_tracking[ref_isoform].append(isoform_id)
            if assignment_type != '=':
                self.incorrect_isoforms.append((assignment_type, ref_isoform, isoform_id))

        self.all_ref_isoforms = all_ref_tracking.keys()
        for ref_iso in all_ref_tracking.keys():
            isoform_list = all_ref_tracking[ref_iso]
            if len(isoform_list) > 1:
                self.duplicate_isoforms[ref_iso] = isoform_list

        self.incorrect_isoforms = sorted(self.incorrect_isoforms )


class AssignmentData:
    def __init__(self, isoquant_output_prefix):
        self.isoform_to_read = defaultdict(list)
        self.assigned_reads = defaultdict(dict)
        self.parse_assignments(isoquant_output_prefix + "read_assignments.tsv")
        self.parse_isoform_map(isoquant_output_prefix + "transcript_models_reads.tsv")

    def parse_assignments(self, tsv_file):
        logger.info("Reading assignments from %s" % tsv_file)
        for l in open(tsv_file):
            if l.startswith('#'):
                continue
            tokens = l.strip().split()
            read_id = tokens[0]
            original_isoform_id = get_original_isoform(read_id)
            info = '\t'.join(tokens[4:8])
            self.assigned_reads[original_isoform_id][read_id] = info

    def parse_isoform_map(self, tsv_file):
        logger.info("Reading isoform to read id map from %s" % tsv_file)
        for l in open(tsv_file):
            if l.startswith('#'):
                continue
            tokens = l.strip().split()
            read_id = tokens[0]
            isoform_id = tokens[1]
            if isoform_id != '*':
                self.isoform_to_read[isoform_id].append(read_id)


class GeneDBHandler:
    def __init__(self, gene_db, output_prefix, dbname, complete_genedb=False):
        self.gene_db = gene_db
        if not gene_db.lower().endswith('db'):
            self.gene_db = output_prefix + '.' + dbname + '.db'
            if not os.path.exists(self.gene_db):
                self.gtf2db(gene_db, self.gene_db, complete_genedb)
        self.db = gffutils.FeatureDB(self.gene_db, keep_order=True)
        self.isoform_to_exon = defaultdict(str)
        self.parse_db()

    def gtf2db(self, gtf, db, complete_db=False):
        logger.info("Converting gene annotation file to .db format (takes a while)...")
        gffutils.create_db(gtf, db, force=True, keep_order=True, merge_strategy='merge',
                           sort_attribute_values=True, disable_infer_transcripts=complete_db,
                           disable_infer_genes=complete_db)
        logger.info("Gene database written to " + db)

    def parse_db(self):
        logger.info("Loading gene database")
        for t in self.db.features_of_type(featuretype=('transcript', 'mRNA')):
            exon_strs = []
            for e in self.db.children(t, featuretype=('exon')):
                exon_strs.append("%d-%d" % (e.start, e.end))
            self.isoform_to_exon[t.id] = ','.join(exon_strs)
        logger.info("Gene database loaded with %d transcripts" % len(self.isoform_to_exon))


class StatCounter:
    def __init__(self, ref_db, isoquant_db, isoquant_data, gff_compare_data, output_prefix):
        self.missed_isoforms_file = output_prefix + ".missed_isoforms.tsv"
        self.duplicated_isoforms_file = output_prefix + ".duplicated_isoforms.tsv"
        self.incorrect_isoforms_file = output_prefix + ".incorrect_isoforms.tsv"
        self.unmapped_isoforms_file = output_prefix + ".unmapped_isoforms.tsv"
        self.ref_db = ref_db
        self.isoquant_db = isoquant_db
        self.isoquant_data = isoquant_data
        self.gff_compare_data = gff_compare_data

    def process_all(self):
        self.process_incorrect()
        self.process_missing()
        self.process_unmapped()
        self.process_duplicated()

    def process_incorrect(self):
        logger.info("Saving incorrect isoforms to %s" % self.incorrect_isoforms_file)
        with open(self.incorrect_isoforms_file, 'w') as outf:
            for info_tuple in self.gff_compare_data.incorrect_isoforms:
                (assignment_type, ref_isoform, isoform_id) = info_tuple
                outf.write('\n= Incorrect isoform, assignment type %s\n' % assignment_type)
                outf.write('%s\t%s\n' % (ref_isoform, isoform_id))
                outf.write(self.ref_db.isoform_to_exon[ref_isoform] + '\n')
                outf.write(self.isoquant_db.isoform_to_exon[isoform_id] + '\n')
                outf.write("Reads contributed: %d\n" % len(self.isoquant_data.isoform_to_read[isoform_id]))
                for read_id in self.isoquant_data.isoform_to_read[isoform_id]:
                    original_isoform = get_original_isoform(read_id)
                    outf.write(read_id + '\t' + self.isoquant_data.assigned_reads[original_isoform][read_id] + '\n')

    def process_missing(self):
        logger.info("Saving missing isoforms to %s" % self.missed_isoforms_file)
        with open(self.missed_isoforms_file, 'w') as outf:
            for ref_isoform_id in self.ref_db.isoform_to_exon.keys():
                if ref_isoform_id in self.gff_compare_data.all_ref_isoforms:
                    continue
                outf.write('\n= Missed isoform %s\n' % ref_isoform_id)
                outf.write(self.ref_db.isoform_to_exon[ref_isoform_id] + '\n')
                all_reads = self.isoquant_data.assigned_reads[ref_isoform_id]
                outf.write("Original reads : %d\n" % len(all_reads))
                for read_id in all_reads:
                    outf.write(read_id + '\t' + all_reads[read_id] + '\n')

    def process_unmapped(self):
        logger.info("Saving unmapped isoforms to %s" % self.unmapped_isoforms_file)
        with open(self.unmapped_isoforms_file, 'w') as outf:
            for isoform_id in self.gff_compare_data.unmapped_isoforms:
                outf.write('\n= Unmapped isoform %s\n' % isoform_id)
                outf.write(self.isoquant_db.isoform_to_exon[isoform_id] + '\n')
                outf.write("Reads contributed: %d\n" % len(self.isoquant_data.isoform_to_read[isoform_id]))
                for read_id in self.isoquant_data.isoform_to_read[isoform_id]:
                    original_isoform = get_original_isoform(read_id)
                    outf.write(read_id + '\t' + self.isoquant_data.assigned_reads[original_isoform][read_id] + '\n')

    def process_duplicated(self):
        logger.info("Saving duplicated isoforms to %s" % self.duplicated_isoforms_file)
        with open(self.duplicated_isoforms_file, 'w') as outf:
            for ref_isoform_id in self.gff_compare_data.duplicate_isoforms:
                outf.write('\n= Duplicated reference isoform %s\n' % ref_isoform_id)
                outf.write(self.isoquant_db.isoform_to_exon[ref_isoform_id] + '\n')
                for i, isoform_id in enumerate(self.gff_compare_data.duplicate_isoforms[ref_isoform_id]):
                    outf.write("> Isoform #%d, %s, reads contributed: %d\n" %
                               (i, isoform_id, len(self.isoquant_data.isoform_to_read[isoform_id])))
                    for read_id in self.isoquant_data.isoform_to_read[isoform_id]:
                        original_isoform = get_original_isoform(read_id)
                        outf.write(read_id + '\t' + self.isoquant_data.assigned_reads[original_isoform][read_id] + '\n')


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", default="gff_stats")
    parser.add_argument("--isoquant_prefix", "-i", type=str, help="path to IsoQuant output")
    parser.add_argument("--gffcompare_tracking", "-t", type=str, help="gffcompare tracking output")
    parser.add_argument("--genedb", "-g", type=str, help="reference gene database")
    parser.add_argument("--isoquantdb", "-q", type=str, help="IsoQuant gene database used for gffcompare")


    args = parser.parse_args()
    if not args.isoquant_prefix or not args.gffcompare_tracking or not args.genedb or not args.isoquantdb:
        parser.print_usage()
        exit(-1)
    return args


def main():
    set_logger(logger)
    args = parse_args()
    ref_db = GeneDBHandler(args.genedb, args.output, 'reference', True)
    isoquant_db = GeneDBHandler(args.isoquantdb, args.output, 'isoquant', False)
    isoquant_data = AssignmentData(args.isoquant_prefix)
    gff_compare_data = TrackingData(args.gffcompare_tracking)
    stat_counter = StatCounter(ref_db, isoquant_db, isoquant_data, gff_compare_data, args.output)
    stat_counter.process_all()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



