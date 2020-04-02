############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import copy
import gffutils
import argparse
import pysam
import glob
from common import *
from traceback import print_exc
from functools import partial

# global params, to be fixed
RESOLVE_AMBIGUOUS = False
DEDUCE_CODONS_FROM_CDS = False
READS_CUTOFF = 10
MIN_CODON_COUNT = 2
ASSIGN_CODONS_WHEN_AMBIGUOUS = False
CONSIDER_FLANKING_JUNCTIONS = False
JUNCTION_DELTA = 1
LR_JUNCTION_DELTA = 2
COUNT_ISOFORM_STATS = True
# merge --- merge overlaping genes and treat as one
# separate --- count start/stop codons independently for each gene
# ignore_overlaps --- do not output overlapping genes at all
# shared_exons --- merge ones with at least 2 shared exons
CODON_OUTPUT = "separate"
WRITE_CODON_COORDINATES = True

# global variables for carrying out the stats
global_assignment_map = {}
global_unassignable_set = set()


def add_to_global_stats(read_id, mathched_isoforms):
    if read_id not in global_assignment_map:
        global_assignment_map[read_id] = set()
    for i in mathched_isoforms:
        global_assignment_map[read_id].add(i)

DEBUG = False
def print_debug(s):
    if DEBUG:
        print(s)

# class for saving all the stats
class ReadAssignmentStats:
    low_covered = 0
    uniquely_assigned = 0
    unique_extra_exon = 0
    unique_extra_intros = 0
    assigned_to_ncrna = 0
    contradictory = 0
    empty = 0
    ambiguous = 0
    ambiguous_codon_assigned = 0
    ambiguous_subisoform_assigned = 0
    ambiguous_unassignable = 0

    correctly_assigned = 0
    unassigned = 0
    mismapped = 0
    unmapped = 0
    empty_bc = 0
    incorrectly_assigned_nc = 0
    unassigned_nc = 0
    incorrectly_assigned_same_gene = 0
    incorrectly_assigned_other_gene = 0
    unassignable = 0


    def __init__(self):
        self.low_covered = 0
        self.uniquely_assigned = 0
        self.unique_extra_exon = 0
        self.unique_extra_intros = 0
        self.assigned_to_ncrna = 0
        self.contradictory = 0
        self.empty = 0
        self.ambiguous = 0
        self.ambiguous_codon_assigned = 0
        self.ambiguous_subisoform_assigned = 0
        self.ambiguous_unassignable = 0

        self.correctly_assigned = 0
        self.unassigned = 0
        self.mismapped = 0
        self.unmapped = 0
        self.empty_bc = 0
        self.incorrectly_assigned_nc = 0
        self.unassigned_nc = 0
        self.incorrectly_assigned_same_gene = 0
        self.incorrectly_assigned_other_gene = 0
        self.unassignable = 0

    def merge(self, stat):
        self.low_covered += stat.low_covered
        self.uniquely_assigned += stat.uniquely_assigned
        self.unique_extra_exon += stat.unique_extra_exon
        self.unique_extra_intros += stat.unique_extra_intros
        self.assigned_to_ncrna += stat.assigned_to_ncrna
        self.contradictory += stat.contradictory
        self.empty += stat.empty
        self.ambiguous += stat.ambiguous
        self.ambiguous_codon_assigned += stat.ambiguous_codon_assigned
        self.ambiguous_subisoform_assigned += stat.ambiguous_subisoform_assigned
        self.ambiguous_unassignable += stat.ambiguous_unassignable

        self.correctly_assigned += stat.correctly_assigned
        self.unassigned += stat.unassigned
        self.mismapped += stat.mismapped
        self.unmapped += stat.unmapped
        self.empty_bc += stat.empty_bc
        self.incorrectly_assigned_nc += stat.incorrectly_assigned_nc
        self.unassigned_nc += stat.unassigned_nc
        self.incorrectly_assigned_same_gene += stat.incorrectly_assigned_same_gene
        self.incorrectly_assigned_other_gene += stat.incorrectly_assigned_other_gene
        self.unassignable += stat.unassignable

    def isoform_stats(self):

        total = self.correctly_assigned + self.unassigned + self.mismapped + self.unmapped +\
                self.incorrectly_assigned_same_gene + self.incorrectly_assigned_other_gene +  self.empty_bc+ \
                self.incorrectly_assigned_nc + self.unassigned_nc + self.unassignable
        s = "\nTotal\tcorrect\twrong_same\twrong_other\tunassigned\tmismapped\tunmapped\tempty\t\twrong_nc\tunassigned_nc\tunassignable\n"
        return s + "%d\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d" % \
            (total, self.correctly_assigned, self.incorrectly_assigned_same_gene, self.incorrectly_assigned_other_gene,
             self.unassigned, self.mismapped, self.unmapped, self.empty_bc, self.incorrectly_assigned_nc,
             self.unassigned_nc, self.unassignable)

    def to_str(self):
        total_bc = self.low_covered + self.uniquely_assigned + self.unique_extra_exon + self.unique_extra_intros +\
                   self.assigned_to_ncrna + self.contradictory + self.empty + self.ambiguous + self.ambiguous_codon_assigned + \
                   self.ambiguous_subisoform_assigned + self.ambiguous_unassignable
        s = "\nTotal\t\tlow_covered\tunique\tunique_ee\tunique_ei\t\tncrna\t\tcontradictory\tempty\t\tambiguous\tambiguous_codon\tambiguous_assigned\tunassignable\n"
        return s + "%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t\t%d" % \
            (total_bc, self.low_covered, self.uniquely_assigned, self.unique_extra_exon, self.unique_extra_intros,
             self.assigned_to_ncrna, self.contradictory, self.empty, self.ambiguous, self.ambiguous_codon_assigned,
             self.ambiguous_subisoform_assigned, self.ambiguous_unassignable)


# class for storing support vectors for known features (junctions)
class FeatureVector:
    profile = []
    reads = 0
    check_flanking = False
    fill_gaps = True
    block_comparator = None
    
    def __init__(self, num, check_flanking, fill_gaps, block_comparator, ignore_flanking_blocks = True):
        self.profile = [0 for i in range(0, num)]  
        self.reads = 0
        self.check_flanking = check_flanking
        self.fill_gaps = fill_gaps
        self.block_comparator = block_comparator
        self.ignore_flanking_blocks = ignore_flanking_blocks

    # update vector using features from alignment
    def add_from_blocks(self, read_features, known_features): 
        read_pos = 0
        ref_pos = 0
        
        #print read_features
        #print known_features

        features_present = [0 for i in range(0, len(known_features) + 2)]
        starting_matched = False
        terminating_matched = False

        if self.check_flanking:
            if len(read_features) > 0 and len(known_features) > 0  and left_of(read_features[0], known_features[0]):
                features_present[0] = 1

            if len(read_features) > 0 and len(known_features) > 0 and left_of(known_features[-1], read_features[-1]):
                features_present[-1] = 1

        while read_pos < len(read_features) and ref_pos < len(known_features):
            while read_pos < len(read_features) and left_of(read_features[read_pos], known_features[ref_pos]):
                read_pos += 1
            if read_pos == len(read_features):
                break

            while ref_pos < len(known_features) and left_of(known_features[ref_pos], read_features[read_pos]):
                ref_pos += 1
            if ref_pos == len(known_features):
                break

            if self.block_comparator(known_features[ref_pos], read_features[read_pos]):
                if read_pos == 0:
                    starting_matched = True
                elif read_pos ==  len(read_features) - 1:
                    terminating_matched = True
                features_present[ref_pos + 1] = 1
                ref_pos += 1
            elif overlaps(known_features[ref_pos], read_features[read_pos]):
                if 0 < read_pos < len(read_features) - 1:
                    features_present[ref_pos + 1] = -1
                elif read_pos == 0 and equal_right_border(known_features[ref_pos], read_features[read_pos]):
                    features_present[ref_pos + 1] = -1
                elif read_pos == len(read_features) - 1 and equal_left_border(known_features[ref_pos], read_features[read_pos]):
                    features_present[ref_pos + 1] = -1
                ref_pos += 1
            elif left_of(known_features[ref_pos], read_features[read_pos]):
                ref_pos += 1
            else:
                read_pos +=1

        # filling gaps between included / excluded features
        if self.fill_gaps:
            self.fill_gaps_in_profile(features_present)

        # if ignoring inexact 1st and last blocks, mark overlapping blocks with 0 (as unknown) in case not exact match is found
        if self.ignore_flanking_blocks and not starting_matched:
            for ref_pos in range(len(known_features)):
                if overlaps(known_features[ref_pos], read_features[0]):
                    features_present[ref_pos + 1] = 0
        if self.ignore_flanking_blocks and not terminating_matched:
            for ref_pos in range(len(known_features)):
                if overlaps(known_features[ref_pos], read_features[-1]):
                    features_present[ref_pos + 1] = 0

        self.reads += 1
        for i in range(0, len(self.profile)):
            self.profile[i] += features_present[i]

    def fill_gaps_in_profile(self, features_present):
        start = 0
        while start < len(features_present) and features_present[start] == 0:
            start += 1

        end = len(features_present) - 1
        while end > 0 and features_present[end] == 0:
            end -= 1

        for i in range(start, end):
            if features_present[i] == 0:
                features_present[i] = -1

# Two feature vectors (intons + exons) + support information for a read/barcode mapping
class ReadMappingInfo:
    read_id = ""
    total_reads = 0
    junctions_counts = None
    exons_counts = None

    def __init__(self, read_id, introns_count, exons_count,
                 check_flanking = CONSIDER_FLANKING_JUNCTIONS,
                 exon_counting_mode = False, delta = 1):
        self.read_id = read_id
        self.exon_counting_mode = exon_counting_mode
        self.total_reads = 0
        self.junctions_counts = \
            FeatureVector(introns_count, check_flanking = check_flanking, fill_gaps = True, block_comparator = equal_ranges)
        self.exons_counts = FeatureVector(exons_count, check_flanking = False, fill_gaps = True,
                                          block_comparator = overlaps, ignore_flanking_blocks = True) \
            if self.exon_counting_mode else \
            FeatureVector(exons_count, check_flanking = check_flanking, fill_gaps = True, block_comparator = overlaps)

    def concat_gapless_blocks(self, blocks, cigar_list):
        cigar_index = 0
        block_index = 0
        resulting_blocks = []

        current_block = None
        deletions_before_block = 0

        while cigar_index < len(cigar_list) and block_index < len(blocks):
            #init new block
            if current_block is None:
                #init new block from match
                if cigar_list[cigar_index][0] == 0:
                    current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                    deletions_before_block = 0
                    block_index += 1
                # keep track of deletions before matched block
                elif cigar_list[cigar_index][0] == 2:
                    deletions_before_block = cigar_list[cigar_index][1]
            # found intron, add current block
            elif cigar_list[cigar_index][0] == 3:
                resulting_blocks.append(current_block)
                current_block = None
            # add deletion to block
            elif cigar_list[cigar_index][0] == 2:
                current_block = (current_block[0], current_block[1] + cigar_list[cigar_index][1])
            # found match - merge blocks
            elif cigar_list[cigar_index][0] == 0:
                if abs(current_block[1] - blocks[block_index][0]) > 1:
                    print("Distant blocks")
                    print(current_block, blocks[block_index])
                current_block = (current_block[0], blocks[block_index][1])

                block_index += 1
            cigar_index += 1

        if current_block is not None:
            resulting_blocks.append(current_block)

        return resulting_blocks

    def add_read(self, alignment, known_introns, known_exons):
        self.total_reads += 1
        #converting to 1-based coordinates
        #second coordinate is not converted since alignment block is end-exclusive, i.e. [x, y)
        blocks = map(lambda x: (x[0] + 1, x[1]), self.concat_gapless_blocks(alignment.get_blocks(), alignment.cigartuples))

#        if alignment.query_name == 'e44d07f3-1abc-41f8-918a-8547c40d511a' or alignment.query_name == 'b6922288-7e4c-4aec-a6dd-a440fc19837f':
#            print(blocks)

        if len(blocks) >= 2 and not self.exon_counting_mode:
            read_junctions = junctions_from_blocks(blocks)
            self.junctions_counts.add_from_blocks(read_junctions, known_introns)
            print_debug("ID " + self.read_id + "\n" + str(read_junctions))
        self.exons_counts.add_from_blocks(blocks, known_exons)



# storage for feature profiles of all known isoforms of a gene or a set of overlapping genes
class IsoformProfileStorage:
    intron_profiles = {}
    exon_profiles = {}
    empty = set()
    ambiguous = set()

    def __init__(self):
        self.intron_profiles = {}
        self.exon_profiles = {}
        self.ambiguous = set()

    # detect isoforms which are exact sub-isoforms of others
    def detect_ambiguous(self):
        for t in self.intron_profiles:
            intron_profile = self.intron_profiles[t]
            exon_profile = self.exon_profiles[t]
            for t2 in self.intron_profiles:
                if is_subprofile(intron_profile, self.intron_profiles[t2]) and \
                        is_subprofile(exon_profile, self.exon_profiles[t2]):
                    self.ambiguous.add(t)
                    #print("Unassignable " + t)
                    #print(isoform_profile)
                    #print(self.isoform_profiles[t2])
                    break
        for i in self.ambiguous:
            global_unassignable_set.add(i)
        #print(self.ambiguous)


# All gene(s) information
class GeneInfo:
    # list of genes in cluster
    gene_db_list = []
    # gene region
    chr_id = None
    start = 0
    end = 0
    # gffutils main structure
    db = None

    # profiles for isoforms having CDS
    coding_rna_profiles = IsoformProfileStorage()
    # profiles for all known isoforoms
    all_rna_profiles = IsoformProfileStorage()
    # sorted list of all known introns (pairs of coordinates)
    introns = []
    # sorted list of exons or non-overlapping exon regions (depends on split_exons option)
    exons = []

    # args: parameters of the tool
    # chr_bam_prefix: additional string used when bam files were aligned to different reference that has difference in chromosome names (e.g. 1 and chr1)
    def __init__(self, gene_db_list, db, chr_bam_prefix = "", split_exons = True):
        self.db = db
        self.gene_db_list = gene_db_list
        self.chr_id, self.start, self.end = self.get_gene_region()
        self.chr_id = chr_bam_prefix + self.chr_id

        self.introns = []
        self.exons = []
        self.coding_rna_profiles = IsoformProfileStorage()
        self.all_rna_profiles = IsoformProfileStorage()

        self.exon_to_geneid = {}
        self.terminal_exons = set()
        self.terminal_only_exons = set()
        i_introns, i_exons = self.set_introns_and_exons(True)
        self.need_to_split_exons = split_exons
        if self.need_to_split_exons:
            self.exons = self.split_exons(self.exons)
        print_debug(self.exons)

        self.set_junction_profiles(self.coding_rna_profiles, i_introns, i_exons, False)
        self.set_junction_profiles(self.all_rna_profiles, i_introns, i_exons, True)

        self.coding_rna_profiles.detect_ambiguous()
        self.all_rna_profiles.detect_ambiguous()
        # print("Gene has " + str(len(self.all_rna_profiles.ambiguous)) + " ambiguous isoforms")

    # return start-stop codon pair for a known isoform
    # return None when codon could not be found
    def get_codon_pair(self, transcript):
        start_codon = None
        stop_codon = None
        for s in self.db.children(transcript, featuretype='start_codon', order_by='start'):
            start_codon = s.start
            if s.strand == "+":
                break
        for s in self.db.children(transcript, featuretype='stop_codon', order_by='start'):
            stop_codon = s.start
            if s.strand == "-":
                break

        if not DEDUCE_CODONS_FROM_CDS:
            return start_codon, stop_codon

        if start_codon is None:
            for s in self.db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    start_codon = s.start
                    break
                else:
                    start_codon = s.end
        if stop_codon is None:
            for s in self.db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    stop_codon = s.end + 1
                else:
                    stop_codon = s.start - 2
                    break
        return start_codon, stop_codon

    def isoform_is_coding(self, t):
        start_codon, stop_codon = self.get_codon_pair(t)
        return stop_codon is not None and start_codon is not None

    # assigns an ordered list of all known exons and introns to self.exons and self.introns
    # returns 2 maps, isoform id -> intron / exon list
    def set_introns_and_exons(self, keep_isoforms_without_codons):
        # dictionary: isoform id -> ordered list of intron coordinates
        all_isoforms_introns = {}
        # dictionary: isoform id -> ordered list of exon coordinates
        all_isoforms_exons = {}
        self.exon_to_geneid = {}

        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype='transcript', order_by='start'):
                if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                    continue

                all_isoforms_exons[t.id] = []
                for e in self.db.children(t, order_by='start'):
                    if e.featuretype == 'exon':
                        exon_block = (e.start, e.end)
                        all_isoforms_exons[t.id].append(exon_block)
                        if exon_block not in self.exon_to_geneid:
                            self.exon_to_geneid[exon_block] = set()
                        self.exon_to_geneid[exon_block].add(gene_db.id)

                all_isoforms_introns[t.id] = junctions_from_blocks(all_isoforms_exons[t.id])

        self.introns = set()
        self.exons = set()
        self.terminal_exons = set()
        for i in all_isoforms_introns.keys():
            self.introns.update(all_isoforms_introns[i])
        for i in all_isoforms_exons.keys():
            self.exons.update(all_isoforms_exons[i])
            self.terminal_exons.add(all_isoforms_exons[i][0])
            self.terminal_exons.add(all_isoforms_exons[i][-1])

        self.terminal_only_exons = copy.deepcopy(self.terminal_exons)
        for i in all_isoforms_exons.keys():
            non_terminal = set(all_isoforms_exons[i][1:-1])
            for terminal in self.terminal_exons:
                if terminal in non_terminal and terminal in self.terminal_only_exons:
                    self.terminal_only_exons.remove(terminal)


        self.introns = sorted(list(self.introns))
        self.exons = sorted(list(self.exons))
        print_debug(self.introns)
        print_debug(self.exons)

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
            print("Unequal number of starts and ends")

        if exon_blocks != sorted(exon_blocks):
            print("Somehow block are unsorted")

        return exon_blocks


    # calculate junction profiles for known isoforms
    def set_junction_profiles(self, profile_storage, i_introns, i_exons, keep_isoforms_without_codons):
        profile_storage.intron_profiles = {}
        profile_storage.exon_profiles = {}

        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype='transcript', order_by='start'):
                if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                    continue

                # setting up intron profiles for current isoform
                profile_storage.intron_profiles[t.id] = [-1 for i in range(0, len(self.introns) + 2)]
                pos = 0
                for intron in i_introns[t.id]:
                    while not equal_ranges(self.introns[pos], intron, 0):
                        pos += 1
                    profile_storage.intron_profiles[t.id][pos + 1] = 1

                # setting up exon profiles for current isoform
                profile_storage.exon_profiles[t.id] = [-1 for i in range(0, len(self.exons) + 2)]
                pos = 0
                for exon in i_exons[t.id]:
                    while not contains(exon, self.exons[pos]):
                        pos += 1
                    while pos < len(self.exons) and contains(exon, self.exons[pos]):
                        profile_storage.exon_profiles[t.id][pos + 1] = 1
                        pos += 1

                print_debug("Isoform " + t.id)
                print_debug(profile_storage.intron_profiles[t.id])
                print_debug(profile_storage.exon_profiles[t.id])

    # compute start-stop codon pair for known isoforms
    def get_codon_pairs(self, gene_db_list):
        codon_pairs = {}
        for gene_db in gene_db_list:
            for t in self.db.children(gene_db, featuretype='transcript', order_by='start'):
                start_codon, stop_codon = self.get_codon_pair(t)
                if stop_codon is not None and start_codon is not None:
                    codon_pairs[t.id] = (start_codon, stop_codon)
                else:
                    codon_pairs[t.id] = (None, None)
        return codon_pairs

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

# storage for all barcodes/reads mapped to a specific set of genes
class ReadProfilesInfo:
    gene_info = None
    # isoform id -> start-stop codon pair (pair of starting coordinates)
    codon_pairs = {}
    # read id -> read mapping info
    read_mapping_infos = {}
    args = None

    def __init__(self, gene_db_list, db, args, chr_bam_prefix = "", exon_count_mode = False):
        self.exon_count_mode = exon_count_mode
        self.gene_info = GeneInfo(gene_db_list, db, chr_bam_prefix, split_exons=not self.exon_count_mode)
        self.args = args
        self.codon_pairs = self.gene_info.get_codon_pairs(self.gene_info.gene_db_list)
        self.read_mapping_infos = {}

    # process alignment within a gene region
    def add_read(self, alignment, read_id):
        chr_name = alignment.reference_name.strip()
        gene_chr = self.gene_info.chr_id
        if not chr_name.startswith('chr') and self.gene_info.chr_id.startswith('chr'):
            gene_chr = self.gene_info.chr_id[3:]
        if chr_name != gene_chr:
            print("Unequal chromosomes")
            return

        blocks = alignment.get_blocks()
#        if read_id == 'e44d07f3-1abc-41f8-918a-8547c40d511a':
#            print(len(alignment.cigartuples))
#            print(alignment.cigartuples)
#            print(len(alignment.get_blocks()))
#            print(sorted(alignment.get_blocks()) == blocks)
#            print(alignment.get_blocks())

        if len(blocks) == 0:
            return

        read_start = blocks[0][0]
        read_end = blocks[-1][1]
        if not overlaps((read_start, read_end), (self.gene_info.start, self.gene_info.end)):
            return

        if read_id not in self.read_mapping_infos:
            self.read_mapping_infos[read_id] = \
                ReadMappingInfo(read_id, len(self.gene_info.introns) + 2, len(self.gene_info.exons) + 2,
                                check_flanking=CONSIDER_FLANKING_JUNCTIONS,
                                exon_counting_mode=self.exon_count_mode, delta=self.args.delta)
        self.read_mapping_infos[read_id].add_read(alignment, self.gene_info.introns, self.gene_info.exons)

    # match barcode/sequence junction profile to a known isoform junction profile, hint - potential candidates
    def find_matching_isofoms_by_profile(self, read_profile, profiles, hint = set()):
        read_sign_profile = map(sign, read_profile.profile)
        print_debug(read_sign_profile)

        matched_isoforms = set()
        for t in profiles.keys():
            if len(hint) > 0 and t not in hint:
                continue
            isoform_profile = profiles[t]
            if diff_only_present(isoform_profile, read_sign_profile) == 0:
                print_debug("Matched " + t)
                matched_isoforms.add(t)
        return matched_isoforms

    # returns true if alignment has no splice junctions and alignment outside of the gene region
    def is_empty_alignment(self, read_mapping_info):
        read_intron_profile = map(sign, read_mapping_info.junctions_counts.profile)
        read_exon_profile = map(sign, read_mapping_info.exons_counts.profile)
#        print_debug(read_mapping_info.read_id)
#        print_debug(read_intron_profile)
#        print_debug(read_exon_profile)

        if all(el != 1 for el in read_intron_profile) and all(el != 1 for el in read_exon_profile[1:-1]):
            return True
        return False

    # resolve assignment ambiguities caused by identical profiles
    def resolve_ambiguous(self, read_count_profile, isoform_profiles, matched_isoforms):
        read_profile = map(sign, read_count_profile)
        if all(el != 1 for el in read_profile[1:-1]):
            return matched_isoforms

        for t in matched_isoforms:
            matched_positions = find_matching_positions(isoform_profiles[t], read_profile)
            # print_debug(matched_positions)
            # print_debug(profile_storage.intron_profiles[t])

            all_positions_detected = True
            for i in range(len(matched_positions)):
                if matched_positions[i] == 0 and isoform_profiles[t][i] != -1:
                    all_positions_detected = False
                    break

            if all_positions_detected:
                # print_debug("Ambiguity resolved")
                # print_debug(profile_storage.intron_profiles[t])
                # print_debug(matched_positions)
                # print_debug(bacrode_jprofile)
                return set([t])

        return matched_isoforms

    # assign barcode/sequence alignment to a known isoform
    def assign_isoform(self, read_id, stat, coverage_cutoff):
        print_debug('=== ' + read_id + ' ===')

        read_mapping_info = self.read_mapping_infos[read_id]
        if read_mapping_info.total_reads < coverage_cutoff:
            stat.low_covered += 1
            return None, None          

        if self.is_empty_alignment(read_mapping_info):
            stat.empty += 1
            print_debug("Empty profile ")
            return None, None          

        # match read profile with isoform profiles
        intron_matched_isoforms = self.find_matching_isofoms_by_profile(read_mapping_info.junctions_counts,
                                                                        self.gene_info.all_rna_profiles.intron_profiles)
        exon_matched_isoforms = self.find_matching_isofoms_by_profile(read_mapping_info.exons_counts,
                                                                      self.gene_info.all_rna_profiles.exon_profiles)
        both_mathched_isoforms = intersection(intron_matched_isoforms, exon_matched_isoforms)

        transcript_id = None
        codon_pair = (None, None)
        if len(both_mathched_isoforms) == 1:
            stat.uniquely_assigned += 1
            print_debug("Unique match")
            transcript_id = list(both_mathched_isoforms)[0]
            codon_pair = self.codon_pairs[transcript_id]
            add_to_global_stats(read_id, both_mathched_isoforms)
        elif len(both_mathched_isoforms) == 0:
            if len(intron_matched_isoforms) == 1:
                print_debug("Extra exons, but unique intron profile " + read_id)
                transcript_id = list(intron_matched_isoforms)[0]
                codon_pair = self.codon_pairs[transcript_id]
                stat.unique_extra_exon += 1
                add_to_global_stats(read_id, intron_matched_isoforms)
            elif len(exon_matched_isoforms) == 1:
                print_debug("Extra intron, but unique exon profile " + read_id)
                transcript_id = list(exon_matched_isoforms)[0]
                codon_pair = self.codon_pairs[transcript_id]
                stat.unique_extra_intros += 1
                add_to_global_stats(read_id, exon_matched_isoforms)
            else:
                stat.contradictory += 1
                print_debug("Contradictory " + read_id)
                add_to_global_stats(read_id, [])
        else:
            if RESOLVE_AMBIGUOUS:
                intron_matched_isoforms = \
                    self.resolve_ambiguous(read_mapping_info.junctions_counts.profile,
                                           self.gene_info.all_rna_profiles.intron_profiles, intron_matched_isoforms)
                exon_matched_isoforms = \
                    self.resolve_ambiguous(read_mapping_info.exons_counts.profile,
                                           self.gene_info.all_rna_profiles.exon_profiles, exon_matched_isoforms)

                both_mathched_isoforms = intersection(intron_matched_isoforms, exon_matched_isoforms)

            if len(both_mathched_isoforms) == 1:
                stat.ambiguous_subisoform_assigned += 1
                print_debug("Unique match after resolution")
                transcript_id = list(both_mathched_isoforms)[0]
                codon_pair = self.codon_pairs[transcript_id]
                add_to_global_stats(read_id, both_mathched_isoforms)
            elif len(both_mathched_isoforms) == 0:
                stat.contradictory += 1
                print_debug("Contradictory " + read_id)
                add_to_global_stats(read_id, [])
            else:
                add_to_global_stats(read_id, both_mathched_isoforms)
                codons = set()
                if self.args.assign_codons_when_ambiguous:
                    for t in both_mathched_isoforms:
                        codons.add(self.codon_pairs[t])
                if len(codons) == 1:
                    codon_pair = list(codons)[0]
                    stat.ambiguous_codon_assigned += 1   
                    transcript_id = list(intron_matched_isoforms)[0]
                else:       
                    if any(mi in self.gene_info.all_rna_profiles.ambiguous for mi in both_mathched_isoforms):
                        stat.ambiguous_unassignable += 1
                        print_debug("Unassignable")
                    else:        
                        stat.ambiguous += 1
                        print_debug("Ambigous match")


        return transcript_id, codon_pair



# Class for processing entire bam file agains gene database
class GeneDBProcessor:
    bc_map = None
    db = None
    args = None
    bamfile_names = []
    output_prefix = ""
    stats = None
    chr_bam_prefix = ""
    
    def __init__(self, args, property_getter = None):
        self.args = args
        self.property_getter = property_getter
        self.bc_map = self.get_barcode_map(args.bam)
        self.bamfile_names = []
        for f in glob.glob(args.bam):
            self.bamfile_names.append(f)
        for bamfile_name in self.bamfile_names:
            if not os.path.isfile(bamfile_name):
                raise Exception("BAM file " + bamfile_name + " does not exist")
            bamfile_in = pysam.AlignmentFile(bamfile_name, "rb")
            if not bamfile_in.has_index:
                raise Exception("BAM file " + bamfile_name + " is not indexed, run samtools index")
            if args.change_chr_prefix and bamfile_in.references[0].startswith('chr'):
                print("Changing chomosome prefix")
                self.chr_bam_prefix = 'chr'
            bamfile_in.close()

        if not os.path.isfile(args.genedb):
            raise Exception("Gene database " + args.genedb + " does not exist")
        self.db = gffutils.FeatureDB(args.genedb, keep_order = True)
        
        self.stats = ReadAssignmentStats()
        self.output_prefix = args.output_prefix

    # read barcode map file generated along with a BAM file
    def get_barcode_map(self, bam_file_names):
        if len(bam_file_names) > 1:
            return None
        barcode_map = {}
        contigs_name, ext = os.path.splitext(bam_file_names[0])
        barcode_map_file = contigs_name + "_map.txt"
        if not os.path.isfile(barcode_map_file):
            print("Barcode file was not found")
            return None
        print("Reading barcode file from " + barcode_map_file)
        for line in open(barcode_map_file):
            tokens = line.strip().split("_barcodeIDs_")
            if len(tokens) != 2:    
                #print("Wrong format, _barcodeIDs_ was not found in " + line)
                continue
            barcode_map[tokens[0]] = filter(lambda x: x != "", tokens[1].replace("_", "-").split(','))
        return barcode_map

    # get barcode or sequence id depending on data type
    def get_sequence_id(self, query_name):
        if self.args.data_type == "10x":
            tokens = query_name.strip().split("___")
            if len(tokens) != 2:
                return ""
            return tokens[1]
        elif self.args.data_type == "contigs":
            return query_name.strip().split("_barcodeIDs_")[0]
        else:
            return query_name.strip().split(':')[-1]

    # add isoform stats whem mapping reference sequences
    def count_isoform_stats(self, isoform, read_id, gene_stats, gene_info):
        gene_isoform_ids = set(gene_info.coding_rna_profiles.intron_profiles.keys())
        gene_all_isoform_ids = set(gene_info.all_rna_profiles.intron_profiles.keys())

        original_isoform = read_id.split('_')[0]
        if isoform is not None:
            if original_isoform == isoform:
                gene_stats.correctly_assigned += 1
            elif original_isoform in gene_isoform_ids:
                gene_stats.incorrectly_assigned_same_gene += 1
                print("Incorrect assignment: isoform = " + isoform + " ; sequence = " + read_id + "\n")
            elif original_isoform in gene_all_isoform_ids:
                gene_stats.incorrectly_assigned_nc += 1
            else:
                gene_stats.incorrectly_assigned_other_gene += 1
                print("Alien assignment: isoform = " + isoform + " ; sequence = " + read_id + "\n")

        else:
            if original_isoform in gene_isoform_ids:
                if original_isoform in gene_info.all_rna_profiles.ambiguous:
                    gene_stats.unassignable += 1
                    #print_debug("ID " + read_id + " is unassignable")
                else:
                    gene_stats.unassigned += 1
                    #print_debug("ID " + read_id + " is unassigned")
            elif original_isoform in gene_all_isoform_ids:
                if original_isoform in gene_info.all_rna_profiles.ambiguous:
                    gene_stats.unassignable += 1
                    #print_debug("ID " + read_id + " is unassignable")
                else:
                    gene_stats.unassigned_nc += 1
                    #print_debug("ID " + read_id + " is unassigned")
            elif original_isoform in gene_info.all_rna_profiles.empty:
                gene_stats.empty_bc += 1
            else:
                gene_stats.mismapped += 1

    def count_unmapped_stats(self, read_profiles, gene_stats):
        gene_isoform_ids = set(read_profiles.gene_info.coding_rna_profiles.intron_profiles.keys())
        processed_ids = set()
        for t in read_profiles.read_mapping_infos.keys():
            if not self.bc_map or len(self.bc_map) == 0:
                processed_ids.add(t.split('_')[0])
            else:
                processed_ids.add(self.bc_map[t][0])
        for t in gene_isoform_ids:
            if t not in processed_ids:
                gene_stats.unmapped += 1

    # assign all reads/barcodes mapped to gene region
    def process_all_reads_from_file(self, read_profiles, bamfile_name):
        samfile_in = pysam.AlignmentFile(bamfile_name, "rb")
        gene_chr, gene_start, gene_end = read_profiles.gene_info.get_gene_region()

        # process all alignments
        # prefix is needed when bam file has chrXX chromosome names, but reference has XX names

        if not samfile_in.references[0].startswith('chr') and gene_chr.startswith('chr'):
            gene_chr = gene_chr[3:]

        for alignment in samfile_in.fetch(self.chr_bam_prefix + gene_chr, gene_start, gene_end):
            if alignment.reference_id == -1 or alignment.is_secondary:
                continue
            seq_id = self.get_sequence_id(alignment.query_name)
            if self.args.read_info_map is not None and seq_id in self.args.read_info_map:
                self.found_barcodes.add(self.args.read_info_map[seq_id])
            read_profiles.add_read(alignment, seq_id)
        samfile_in.close()

    def process_all_reads(self, read_profiles):
        for bamfile_name in self.bamfile_names:
            self.process_all_reads_from_file(read_profiles, bamfile_name)

    # assign all reads/barcodes mapped to gene region
    def assign_all_reads(self, read_profiles):
        self.process_all_reads(read_profiles)

        # read / barcode id -> (isoform, codon pair)
        assigned_reads = {}
        gene_stats = ReadAssignmentStats()

        #iterate over all barcodes / sequences and assign them to known isoforms
        for read_id in read_profiles.read_mapping_infos.keys():
            isoform, codons = read_profiles.assign_isoform(read_id, gene_stats, self.args.reads_cutoff)

            seq_id = read_id
            if self.bc_map is not None:
                seq_id = list(self.bc_map[read_id])[0]
                for bc in self.bc_map[read_id]:
                    assigned_reads[bc] = (isoform, codons)
            else:
                assigned_reads[read_id] = (isoform, codons)
            
            if self.args.count_isoform_stats:
                self.count_isoform_stats(isoform, seq_id, gene_stats, read_profiles.gene_info)
        
        if self.args.count_isoform_stats:
            self.count_unmapped_stats(read_profiles, gene_stats)

        print("Done. Read stats " + gene_stats.to_str())
        if self.args.count_isoform_stats:
            print("Done. Isoform stats " + gene_stats.isoform_stats())
        self.stats.merge(gene_stats)

        return assigned_reads

    # calculate exon counts
    def calculate_exon_counts(self, read_profiles):
        self.found_barcodes = set()
        self.process_all_reads(read_profiles)

        # cell group -> two list of counts for each exon inclusion/exclusion
        exon_counts = {}

        # iterate over all barcodes / sequences and assign them to known isoforms
        for read_id in read_profiles.read_mapping_infos.keys():
            group_id = self.property_getter.get_property(read_id)
            if group_id not in exon_counts:
                # first list for exon inclusion, second for exlusion
                exon_counts[group_id] = ([0 for i in range(len(read_profiles.gene_info.exons))],
                                         [0 for i in range(len(read_profiles.gene_info.exons))])

            exon_count_profile = read_profiles.read_mapping_infos[read_id].exons_counts.profile
            for i in range(len(exon_count_profile) - 1):
                if exon_count_profile[i] > 0:
                    exon_counts[group_id][0][i - 1] += exon_count_profile[i]
                elif exon_count_profile[i] < 0:
                    exon_counts[group_id][1][i - 1] -= exon_count_profile[i]
        #print(exon_counts)
        return exon_counts

    def gene_list_id_str(self, gene_db_list, delim = "_"):
        gene_names = [g.id for g in gene_db_list]
        return delim.join(gene_names)
    

    def write_gene_stats(self, gene_db_list, assignments):
        #writing TSV with read id -> isoform id
        outf = open(self.out_tsv, "a+")
        outf.write(self.gene_list_id_str(gene_db_list) + "\t" + str(len(assignments)) + "\n")
        for b in assignments.keys():
            if b != "" and assignments[b][0] is not None:
                outf.write(b + "\t" + assignments[b][0] + "\n")
        outf.close()

    def write_codon_tables(self, gene_db_list, gene_info, assigned_reads):
        if CODON_OUTPUT == "ignore_overlaps" and len(gene_db_list) == 1:
            self.write_codon_tables_for_genes(gene_db_list, gene_info, assigned_reads)
        elif CODON_OUTPUT == "merge":
            self.write_codon_tables_for_genes(gene_db_list, gene_info, assigned_reads)
        elif CODON_OUTPUT == "separate":
            for g in gene_db_list:
                self.write_codon_tables_for_genes([g], gene_info, assigned_reads)
        elif CODON_OUTPUT == "merge_exons":
            for genes in self.group_genes_with_overlapping_exons(gene_db_list):
                self.write_codon_tables_for_genes(genes, gene_info, assigned_reads)
        else:
            print("Unsupported codon output method")


    def group_genes_with_overlapping_exons(self, gene_db_list):
        sys.exit(-1)

    def write_codon_tables_for_genes(self, gene_db_list, gene_info, assigned_reads):
        #writing codon stats
        gene_codon_pairs =  set(gene_info.get_codon_pairs(gene_db_list).values())
        codon_count_table = {}
        start_codons = set()
        stop_codons = set()
        for read_id in assigned_reads.keys():
            if read_id != "" and assigned_reads[read_id][0] is None:
                continue

            codon_pair = assigned_reads[read_id][1]
            if codon_pair is None or codon_pair[0] is None or codon_pair[1] is None:
                continue
            if codon_pair[0] <= 0 or codon_pair[1] <= 0:
                print("Incorrect codon pair")
                continue

            if codon_pair not in gene_codon_pairs:
                continue

            start_codons.add(codon_pair[0])
            stop_codons.add(codon_pair[1])
            if codon_pair not in codon_count_table:
                codon_count_table[codon_pair] = 0
            codon_count_table[codon_pair] += 1

        if len(start_codons) >= MIN_CODON_COUNT and len(stop_codons) >= MIN_CODON_COUNT:
            #WRITE_CODON_COORDINATES
            outf = open(self.out_codon_stats, "a+")
            outf.write("====" + self.gene_list_id_str(gene_db_list) + "\n")
            outf.write(table_to_str(codon_count_table, WRITE_CODON_COORDINATES))
            outf.close()

    def write_exon_counts(self, exon_counts, read_profiles, chr_id, exon_to_genes):
        gene_coverage = 0 if self.args.distinct_barcodes == 0 else float(len(self.found_barcodes)) / self.args.distinct_barcodes

        for i in range(len(read_profiles.gene_info.exons)):
            total_counts = 0
            total_inclusion = 0
            exon = read_profiles.gene_info.exons[i]
#            if not self.args.keep_terminal and exon in read_profiles.gene_info.terminal_exons:
#                continue
            exon_id = chr_id + "_" + str(exon[0]) + "_" + str(exon[1])
            out_exons = open(self.out_exon_counts, "a+")

            exon_type = "."
            if len(exon_to_genes[exon]) > 1:
                exon_type += "M"
            if exon in read_profiles.gene_info.terminal_only_exons:
                exon_type += "X"
            elif exon in read_profiles.gene_info.terminal_exons:
                exon_type += "T"
            if exon_type != ".":
                exon_type = exon_type[1:]

            for group_id in exon_counts.keys():
                if group_id == 'unknown':
                    continue
                include_counts = exon_counts[group_id][0][i]
                exclude_counts = exon_counts[group_id][1][i]
                total_counts += exclude_counts + include_counts
                total_inclusion += include_counts
            #inclusion_rate = 0 if total_counts == 0 else float(total_inclusion) / float(total_counts)

            for group_id in exon_counts.keys():
                if group_id == 'unknown':
                    continue
                include_counts = exon_counts[group_id][0][i]
                exclude_counts = exon_counts[group_id][1][i]
                if exclude_counts == 0 and include_counts == 0:
                    continue
                out_exons.write(exon_id + "\t" + group_id + "\t" + str(include_counts) + "\t" + str(exclude_counts)
                                + "\t" + str(total_counts) + "\t" + str(total_inclusion) + "\t" + exon_type + "\t" +
                                ",".join(list(exon_to_genes[exon])) + "\t" + "{0:.4f}".format(gene_coverage) + "\n")
            out_exons.close()

            #out_exon_info = open(self.out_exon_genes, "a+")
            #out_exon_info.write(exon_id + "\t" + ",".join(list(exon_to_genes[exon])) + "\t" + str(gene_coverage) + "\t" + str(inclusion_rate) + "\n")
            #out_exon_info.close()

    # Process a set of genes given in gene_db_list
    def process_gene_list(self, gene_db_list):
        print("Processing " + str(len(gene_db_list)) + " gene(s): " + self.gene_list_id_str(gene_db_list, ", "))

        read_profiles = ReadProfilesInfo(gene_db_list, self.db, self.args, self.chr_bam_prefix,
                                         exon_count_mode=self.args.exon_count_mode)
        exon_to_genes = read_profiles.gene_info.exon_to_geneid
        if self.args.exon_count_mode:
            for gene_db in gene_db_list:
                single_gene_profiles = ReadProfilesInfo([gene_db], self.db, self.args, self.chr_bam_prefix,
                                         exon_count_mode=self.args.exon_count_mode)
                exon_counts = self.calculate_exon_counts(single_gene_profiles)
                self.write_exon_counts(exon_counts, single_gene_profiles, gene_db_list[0].seqid, exon_to_genes)
        else:
            assigned_reads = self.assign_all_reads(read_profiles)

            self.write_gene_stats(gene_db_list, assigned_reads)
            self.write_codon_tables(gene_db_list, read_profiles.gene_info, assigned_reads)

    #R un though all genes in db and count stats according to alignments given in bamfile_names
    def process_all_genes(self):
        if self.args.exon_count_mode:
            self.out_exon_counts = self.output_prefix + ".exon_counts.tsv"
            outf = open(self.out_exon_counts, "w")
            outf.write('#exon_id\tcell_type\tinclusion\texclusion\ttotal_count\ttotal_inclusion\texon_type\tgene_id\tgene_coverage\n')
            outf.close()
            #self.out_exon_genes = self.output_prefix + ".exon_to_geneid.tsv"
            #outf = open(self.out_exon_genes, "w")
            #outf.write('#exon_id\tgene_id\tgene_coverage\texon_inclusion\n')
            #outf.close()
        else:
            self.out_tsv = self.output_prefix + ".assigned_reads.tsv"
            outf = open(self.out_tsv, "w")
            outf.close()
            self.out_codon_stats = self.output_prefix + ".codon_stats.tsv"
            outf = open(self.out_codon_stats, "w")
            outf.close()
        
        gene_db_list = []
        current_chromosome = ""

        for g in self.db.features_of_type('gene', order_by=('seqid', 'start')):
            if current_chromosome != g.seqid:
                current_chromosome = g.seqid
                print("Processing chromosome " + current_chromosome)
            gene_name = g.id
            gene_db = self.db[gene_name]

            if len(gene_db_list) == 0 or any(genes_overlap(g, gene_db) for g in gene_db_list):
                gene_db_list.append(gene_db)
            else:
                self.process_gene_list(gene_db_list)
                gene_db_list = [gene_db]

        self.process_gene_list(gene_db_list)

        if self.args.exon_count_mode:
            print("Done")
        else:
            print("\nFinished. Total stats " + self.stats.to_str())
            if self.args.count_isoform_stats:
                print("Finished. Isoform stats " + self.stats.isoform_stats() + "\n")


#Print global stats for isoform assignment
def global_stats(bc_map):
    correct_unique = 0
    incorrect_unique = 0
    correct_amb = 0
    incorrect_amb = 0
    correct_unassignable = 0
    incorrect_unassignable = 0
    unassigned_unassignable = 0
    unassigned = 0

    for k in global_assignment_map.keys():
        b = k.split('_')[0]
        if bc_map is not None:
            b = bc_map[k][0]
        matched_isoforms = global_assignment_map[k]
        if len(matched_isoforms) == 0:
            if b in global_unassignable_set:
                unassigned_unassignable += 1
            else:
                unassigned += 1
        elif len(matched_isoforms) == 1:
            if b == list(matched_isoforms)[0]:
                correct_unique += 1
            else:
                incorrect_unique += 1
        else:
            if b in matched_isoforms:
                if b in global_unassignable_set:
                    correct_unassignable += 1
                else:
                    correct_amb += 1
            else:
                if b in global_unassignable_set:
                    incorrect_unassignable += 1
                else:
                    incorrect_amb += 1

    print("\nGlobal stats, total reads / barcodes processed " + str(len(global_assignment_map)))
    print("\t\tCorrect\tIncorrect")
    print("Unique\t\t" + str(correct_unique) + "\t" + str(incorrect_unique))
    print("Ambiguous\t" + str(correct_amb) + "\t" + str(incorrect_amb))
    print("Unassignable\t" + str(correct_unassignable) + "\t" + str(incorrect_unassignable))
    print("\nNot matched " + str(unassigned) + ", unassignable and not matched " + str(unassigned_unassignable))
        

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', metavar='BAM_FILE', nargs='+', type=str,  help='sorted and indexed BAM file(s)')
    parser.add_argument("--data_type", "-d", help="type of data to process, supported types are: contigs, 10x, long_reads, isoforms", type=str, default = "10x")
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    parser.add_argument("--change_chr_prefix", help="change prefix", type=bool, default=False)
    args = parser.parse_args()
    return args

#Tune algorithm params
def set_params(args):
    args.reads_cutoff = READS_CUTOFF if args.data_type == "10x" else 0
    args.assign_codons_when_ambiguous = ASSIGN_CODONS_WHEN_AMBIGUOUS and args.data_type != "isoforms"
    args.consider_flanking_junctions = CONSIDER_FLANKING_JUNCTIONS and args.data_type != "10x" 
    args.delta = LR_JUNCTION_DELTA  if args.data_type == "long_reads" else JUNCTION_DELTA
    args.count_isoform_stats = COUNT_ISOFORM_STATS and args.data_type == "isoforms"
    args.exon_count_mode = False


def main():
    args = parse_args()
    set_params(args)

    db_processor = GeneDBProcessor(args)
    db_processor.process_all_genes()

    if args.count_isoform_stats:
        global_stats(db_processor.bc_map)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

