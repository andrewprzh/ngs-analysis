############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gffutils
from functools import partial

from src.common import *

logger = logging.getLogger('IsoQuant')


# storage for feature profiles of all known isoforms of a gene or a set of overlapping genes
class FeatureProfiles:
    profiles = {}
    features = []

    def __init__(self):
        self.profiles = {}
        self.features = []

    def set_features(self, features):
        self.features = features

    def set_profiles(self, transcript_id, transcript_features, comaprator):
        self.profiles[transcript_id] = [-1] * len(self.features)
        pos = 0

        for feature in transcript_features:
            while pos < len(self.features) and not comaprator(feature, self.features[pos]):
                pos += 1
            while pos < len(self.features) and comaprator(feature, self.features[pos]):
                self.profiles[transcript_id][pos] = 1
                pos += 1


# All gene(s) information
class GeneInfo:
    # chr_bam_prefix: additional string used when bam files were aligned to different reference that has difference in chromosome names (e.g. 1 and chr1)
    def __init__(self, gene_db_list, db):
        # gffutils main structure
        self.db = db
        # list of genes in cluster
        self.gene_db_list = gene_db_list
        # gene region
        self.chr_id, self.start, self.end = self.get_gene_region()

        # profiles for all known isoforoms
        self.intron_profiles = FeatureProfiles()
        self.exon_profiles = FeatureProfiles()
        self.split_exon_profiles = FeatureProfiles()
        self.ambiguous_isoforms = set()

        all_isoforms_introns, all_isoforms_exons = self.set_introns_and_exons()
        self.split_exon_profiles.set_features(self.split_exons(self.exon_profiles.features))

        self.set_junction_profiles(all_isoforms_introns, all_isoforms_exons)
        self.detect_ambiguous()

    # return region of overlapping gene set
    def get_gene_region(self):
        start = self.gene_db_list[0].start
        end = self.gene_db_list[-1].end
        chr_id = self.gene_db_list[0].seqid

        for gene_db in self.gene_db_list:
            if start > gene_db.start:
                start = gene_db.start
            if end < gene_db.end:
                end = gene_db.end

        return chr_id, start, end

    # assigns an ordered list of all known exons and introns to self.exons and self.introns
    # returns 2 maps, isoform id -> intron / exon list
    def set_introns_and_exons(self):
        # dictionary: isoform id -> ordered list of intron coordinates
        all_isoforms_introns = {}
        # dictionary: isoform id -> ordered list of exon coordinates
        all_isoforms_exons = {}

        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype='transcript', order_by='start'):
                all_isoforms_exons[t.id] = []
                for e in self.db.children(t, order_by='start'):
                    if e.featuretype == 'exon':
                        all_isoforms_exons[t.id].append((e.start, e.end))

                all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])

        introns = set()
        exons = set()
        for i in all_isoforms_introns.keys():
            introns.update(all_isoforms_introns[i])
        for i in all_isoforms_exons.keys():
            exons.update(all_isoforms_exons[i])

        self.intron_profiles.set_features(sorted(list(introns)))
        self.exon_profiles.set_features(sorted(list(exons)))

        return all_isoforms_introns, all_isoforms_exons

    # split exons into non-overlapping covering blocks
    def split_exons(self, exons):
        exon_starts = sorted(map(lambda x:x[0], exons))
        exon_ends = sorted(map(lambda x: x[1], exons))

        current_state = 0
        starts_pos = 0
        ends_pos = 0
        last_border = -1
        exon_blocks = []

        while starts_pos < len(exon_starts):
            if exon_starts[starts_pos] <= exon_ends[ends_pos]:
                # do not consider the same exon star
                if starts_pos == 0 or exon_starts[starts_pos] > exon_starts[starts_pos - 1]:
                    cur_border = exon_starts[starts_pos]
                    if last_border != -1 and current_state > 0:
                        exon_blocks.append((last_border, cur_border - 1))
                    last_border = cur_border
                current_state += 1
                starts_pos += 1
            else:
                if last_border == -1 or current_state == 0:
                    print("Error, exon ends before the start")

                if ends_pos == 0 or  exon_ends[ends_pos] >  exon_ends[ends_pos - 1]:
                    cur_border = exon_ends[ends_pos]
                    exon_blocks.append((last_border, cur_border))
                    last_border = cur_border + 1
                current_state -= 1
                ends_pos += 1

        while ends_pos < len(exon_ends):
            if ends_pos == 0 or exon_ends[ends_pos] > exon_ends[ends_pos - 1]:
                cur_border = exon_ends[ends_pos]
                exon_blocks.append((last_border, cur_border))
                last_border = cur_border + 1
            current_state -= 1
            ends_pos += 1

        if current_state != 0:
            logger.critical("Unequal number of starts and ends")

        if exon_blocks != sorted(exon_blocks):
            logger.critical("Somehow block are unsorted")

        return exon_blocks

    # calculate junction profiles for known isoforms
    def set_junction_profiles(self, all_isoforms_introns, all_isoforms_exons):
        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype='transcript', order_by='start'):
                # setting up intron profiles for current isoform
                self.intron_profiles.set_profiles(t.id, all_isoforms_introns[t.id], partial(equal_ranges, delta=0))
                # setting up exon profiles for current isoform
                self.exon_profiles.set_profiles(t.id, all_isoforms_exons[t.id], partial(equal_ranges, delta=0))
                # setting up split exon profiles for current isoform
                self.split_exon_profiles.set_profiles(t.id, all_isoforms_exons[t.id], contains)

    # detect isoforms which are exact sub-isoforms of others
    def detect_ambiguous(self):
        for t in self.intron_profiles.profiles.keys():
            intron_profile = self.intron_profiles.profiles[t]
            exon_profile = self.split_exon_profiles.profiles[t]
            for t2 in self.intron_profiles.profiles.keys():
                if t == t2:
                    continue
                if is_subprofile(intron_profile, self.intron_profiles.profiles[t2]) and \
                        is_subprofile(exon_profile, self.exon_profiles.profiles[t2]):
                    self.ambiguous_isoforms.add(t)
                    break