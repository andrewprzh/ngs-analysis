############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import argparse
import pysam
from common import *
from traceback import print_exc

# global params, to be fixed
RESOLVE_AMBIGUOUS = True
DEDUCE_CODONS_FROM_CDS = True
READS_CUTOFF = 10
MIN_CODON_COUNT = 2
ASSIGN_CODONS_WHEN_AMBIGUOUS = True
CONSIDER_FLANKING_JUNCTIONS = True
JUNCTION_DELTA = 2
LR_JUNCTION_DELTA = 3
COUNT_ISOFORM_STATS = True
# merge --- merge overlaping genes and treat as one
# separate --- count start/stop codons independently for each gene
# ignore_overlaps --- do not output overlapping genes at all
# shared_exons --- merge ones with at least 2 shared exons
CODON_OUTPUT = "separate"
WRITE_CODON_COORDINATES = False

# global variables for carrying out the stats
global_barcode_map = {}
global_unassignable_set = set()

DEBUG = False
def print_debug(s):
    if DEBUG:
        print(s)

# class for storing support vectors for known features (junctions)
class FeatureVector:
    profile = []
    reads = 0
    check_flanking = True
    
    def __init__(self, num, check_flanking = True):
        self.profile = [0 for i in range(0, num)]  
        self.reads = 0
        self.check_flanking = check_flanking    

    # update vector using features from alignment
    def add_from_blocks(self, read_features, known_features): 
        read_pos = 0
        ref_pos = 0
        
        #print read_features
        #print known_features

        features_present = [0 for i in range(0, len(known_features) + 2)]

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

            if equal_ranges(known_features[ref_pos], read_features[read_pos]):
                features_present[ref_pos + 1] = 1
                ref_pos += 1
            elif overlaps(known_features[ref_pos], read_features[read_pos]):
                features_present[ref_pos + 1] = -1
                ref_pos += 1
            elif known_features[ref_pos] < read_features[read_pos]:
                ref_pos += 1
            else:
                read_pos +=1
               
        #print features_present
        self.fill_gaps(features_present)

        self.reads += 1
        for i in range(0, len(self.profile)):
            self.profile[i] += features_present[i]

    def fill_gaps(self, features_present):
        start = 0
        while start < len(features_present) and features_present[start] == 0:
            start += 1

        end = len(features_present) - 1
        while end > 0 and features_present[end] == 0:
            end -= 1

        for i in range(start, end):
            if features_present[i] == 0:
                features_present[i] = -1

# Feature vector + support information for a barcode
class BacrodeInfo:
    barcode = ""
    total_reads = 0
    junctions_counts = None

    def __init__(self, bc, junctions_num, check_flanking):
        self.barcode = bc
        self.total_reads = 0
        self.junctions_counts = FeatureVector(junctions_num, check_flanking)

    def add_read(self, alignment, known_junctions):
        self.total_reads += 1
        blocks = sorted(alignment.get_blocks())
        if len(blocks) >= 2:
            read_junctions = []
            for i in range(0, len(blocks) - 1):
                read_junctions.append((blocks[i][1], blocks[i+1][0]))
            self.junctions_counts.add_from_blocks(read_junctions, known_junctions)
            print_debug("Barcode " + self.barcode + "\n" + str(read_junctions))


# class for saving all the stats
class BarcodeAssignmentStats:
    low_covered = 0
    uniquely_assigned = 0
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

        total = self.correctly_assigned+self.unassigned+self.mismapped+self.unmapped+self.incorrectly_assigned_same_gene+self.incorrectly_assigned_other_gene +  self.empty_bc+ self.incorrectly_assigned_nc + self.unassigned_nc + self.unassignable
        s = "\nTotal\tcorrect\twrong_same\twrong_other\tunassigned\tmismapped\tunmapped\tempty\t\twrong_nc\tunassigned_nc\tunassignable\n"
        return s + "%d\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d" % \
            (total, self.correctly_assigned, self.incorrectly_assigned_same_gene, self.incorrectly_assigned_other_gene, self.unassigned, self.mismapped, self.unmapped, self.empty_bc, self.incorrectly_assigned_nc, self.unassigned_nc, self.unassignable)

    def to_str(self):
        total_bc = self.low_covered + self.uniquely_assigned + self.assigned_to_ncrna + self.contradictory + self.empty + self.ambiguous + self.ambiguous_codon_assigned + self.ambiguous_subisoform_assigned + self.ambiguous_unassignable
        s = "\nTotal\t\tlow_covered\tunique\t\tncrna\t\tcontradictory\tempty\t\tambiguous\tambiguous_codon\tambiguous_assigned\tunassignable\n"
        return s + "%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t\t%d" % \
            (total_bc, self.low_covered, self.uniquely_assigned, self.assigned_to_ncrna, self.contradictory, self.empty, self.ambiguous, self.ambiguous_codon_assigned, self.ambiguous_subisoform_assigned, self.ambiguous_unassignable)

# storage for feature profiles of all known isoforms of a gene or a set of overlapping genes
class ProfileStorage:
    isoform_profiles = {}
    isoform_exon_profiles = {}
    empty = set()
    ambiguous = set()

    def __init__(self):
        self.isoform_profiles = {}
        self.isoform_exon_profiles = {}
        self.empty = set()
        self.ambiguous = set()


    def detect_ambiguous(self):
        for t in self.isoform_profiles:
            isoform_profile = self.isoform_profiles[t]
            for t2 in self.isoform_profiles:
                if t == t2:
                    continue
                if is_subprofile(isoform_profile, self.isoform_profiles[t2]):
                    self.ambiguous.add(t)
                    #print("Unassignable " + t)
                    #print(isoform_profile)
                    #print(self.isoform_profiles[t2])
                    break
        for i in self.ambiguous:
            global_unassignable_set.add(i)
        #print(self.ambiguous)
         
# storage for all barcodes/sequences mapped to a specific set of genes
class GeneBarcodeInfo:
    gene_db_list = []
    chr_id = None
    start = 0
    end = 0
    db = None
    coding_rna_profiles = ProfileStorage()
    all_rna_profiles = ProfileStorage()
    junctions = []
    exons = []
    codon_pairs = {}
    barcodes = {}
    args = None

    def __init__(self, gene_db_list, db, args, chr_bam_prefix = ""):
        self.args = args
        self.codon_pairs = {}
        self.db = db
        self.gene_db_list = gene_db_list
        self.chr_id, self.start, self.end = self.get_gene_region()
        self.chr_id = chr_bam_prefix + self.chr_id
        self.barcodes = {}
        self.junctions = []
        self.exons = []
        self.coding_rna_profiles = ProfileStorage()
        self.all_rna_profiles = ProfileStorage()

        self.set_codon_pairs()
        i_junctions, i_exons = self.get_junctions_and_exons(True)
        self.set_junction_profiles(self.coding_rna_profiles, i_junctions, i_exons, False)
        self.set_junction_profiles(self.all_rna_profiles, i_junctions, i_exons, True)

        self.coding_rna_profiles.detect_ambiguous()
        self.all_rna_profiles.detect_ambiguous()
        #print("Gene has " + str(len(self.all_rna_profiles.ambiguous)) + " ambiguous isoforms")

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

    # get a set of all known exons and splice junctions in a set of genes
    def get_junctions_and_exons(self, keep_isoforms_without_codons):
        i_junctions = {}
        i_exons = {}
        
        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype = 'transcript', order_by='start'):
                if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                    continue

                cur_exon = None
                i_junctions[t.id] = set()
                i_exons[t.id] = set()
                for e in self.db.children(t, order_by='start'):
                    if e.featuretype == 'exon':
                        if cur_exon is None:
                            i_exons[t.id].add((e.start, e.end))
                            cur_exon = e
                            continue
                        i_junctions[t.id].add((cur_exon.end, e.start))
                        i_exons[t.id].add((e.start, e.end))
                        cur_exon = e

        self.junctions = set()
        self.exons = set()
        for i in i_junctions.keys():
            for j in i_junctions[i]:
                self.junctions.add(j)
        for i in i_exons.keys():
            for e in i_exons[i]:
                self.exons.add(e)
            
        self.junctions = sorted(list(self.junctions))
        self.exons = sorted(list(self.exons))
        print_debug(self.junctions)
        #print(self.exons)

        return i_junctions, i_exons

    # calculate junction profiles for known isoforms
    def set_junction_profiles(self, profile_storage, i_junctions, i_exons, keep_isoforms_without_codons):
        profile_storage.isoform_profiles = {}
        profile_storage.isoform_exon_profiles = {}

        for gene_db in self.gene_db_list:
            for t in self.db.children(gene_db, featuretype = 'transcript', order_by='start'):
                if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                    continue

                profile_storage.isoform_profiles[t.id] = [-1 for i in range(0, len(self.junctions) + 2)]
                profile_storage.isoform_exon_profiles[t.id] = [-1 for i in range(0, len(self.exons) + 2)]
                for j in i_junctions[t.id]:
                    pos = self.junctions.index(j)
                    profile_storage.isoform_profiles[t.id][pos + 1] = 1
                for e in i_exons[t.id]:
                    pos = self.exons.index(e)
                    profile_storage.isoform_exon_profiles[t.id][pos + 1] = 1
     
                print_debug("Isoform " + t.id)
                print_debug(profile_storage.isoform_profiles[t.id])
                #print_debug(profile_storage.isoform_exon_profiles[t.id])

                if all(x == -1 for x in profile_storage.isoform_profiles[t.id]):
                    del profile_storage.isoform_profiles[t.id]
                    del profile_storage.isoform_exon_profiles[t.id]
                    profile_storage.empty.add(t.id)

    # compute start-stop codon pair for known isoforms
    def get_codon_pairs(self, gene_db_list):
        codon_pairs = {}
        for gene_db in gene_db_list:
            for t in self.db.children(gene_db, featuretype = 'transcript', order_by='start'):
                start_codon, stop_codon = self.get_codon_pair(t)
                if stop_codon is not None and start_codon is not None:
                    codon_pairs[t.id] = (start_codon, stop_codon)
                else:
                    codon_pairs[t.id] = (None, None)
        return codon_pairs

    def set_codon_pairs(self):
        self.codon_pairs = self.get_codon_pairs(self.gene_db_list)

    # return region of overlapping gene set
    def get_gene_region(self):
        start = self.gene_db_list[0].start
        end = self.gene_db_list[-1].end
        chr_id =  self.gene_db_list[0].seqid
        
        for gene_db in self.gene_db_list:
            if start > gene_db.start:
                start = gene_db.start
            if end < gene_db.end:
                end = gene_db.end

        return chr_id, start, end       

    # process alignment within a gene region
    def add_read(self, alignment, barcode_id):
        chr_name = alignment.reference_name.strip()
        if chr_name != self.chr_id:
            return
        blocks = sorted(alignment.get_blocks())
        if len(blocks) == 0:
            return
        read_start = blocks[0][0]
        read_end = blocks[-1][1]

        if not overlaps((read_start, read_end), (self.start, self.end)):
            return

        if barcode_id not in self.barcodes:
            self.barcodes[barcode_id] = BacrodeInfo(barcode_id, len(self.junctions) + 2, self.args.consider_flanking_junctions)
        self.barcodes[barcode_id].add_read(alignment, self.junctions)

    # match barcode/sequence junction profile to a known isoform junction profile
    def find_matches(self, barcode_info, profile_storage):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
        print_debug(barcode_info.barcode)
        print_debug(bacrode_jprofile)

        matched_isoforms = set()
        for t in profile_storage.isoform_profiles.keys():
            isoform_jprofile = profile_storage.isoform_profiles[t]
            if diff_only_present(isoform_jprofile, bacrode_jprofile) == 0:
                print_debug("Matched " + t)
                matched_isoforms.add(t)
        return matched_isoforms

    # returns true if alignment has no splice junctions (e.g. single-block alignment)
    def is_empty_alignment(self, barcode_info):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
#        print_debug(barcode_info.barcode)
#        print_debug(bacrode_jprofile)

        if all(el != 1 for el in bacrode_jprofile):
            return True
        return False

    # assign barcode/sequence alignment to a known isoform
    def assign_isoform(self, barcode_id, stat, coverage_cutoff):
        print_debug('=== ' + barcode_id + ' ===')

        barcode_info = self.barcodes[barcode_id]
        if barcode_info.total_reads < coverage_cutoff:
            stat.low_covered += 1
            return None, None          

        if self.is_empty_alignment(barcode_info):
            stat.empty += 1
            print_debug("Empty profile ")
            return None, None          

        matched_isoforms = self.find_matches(barcode_info, self.all_rna_profiles)

        if barcode_id not in global_barcode_map:
            global_barcode_map[barcode_id] = set()
        for i in matched_isoforms:
            global_barcode_map[barcode_id].add(i)

        transcript_id = None
        codon_pair = (None, None)
        if len(matched_isoforms) == 0:
            stat.contradictory += 1
            print_debug("Contradictory")
        elif len(matched_isoforms) > 1:
            if RESOLVE_AMBIGUOUS:
                matched_isoforms = self.resolve_ambiguous(barcode_info, matched_isoforms, self.all_rna_profiles)
            
            if len(matched_isoforms) == 1:
                stat.ambiguous_subisoform_assigned += 1
                #print_debug("Unique match after resolution")
                transcript_id = list(matched_isoforms)[0]
                codon_pair = self.codon_pairs[transcript_id]
            else:
                codons = set()
                if self.args.assign_codons_when_ambiguous:
                    for t in matched_isoforms:
                        codons.add(self.codon_pairs[t])
                if len(codons) == 1:
                    codon_pair = list(codons)[0]
                    stat.ambiguous_codon_assigned += 1   
                    transcript_id = list(matched_isoforms)[0]
                else:       
                    if any(mi in self.all_rna_profiles.ambiguous for mi in matched_isoforms):
                        stat.ambiguous_unassignable += 1
                        print_debug("Unassignable")
                    else:        
                        stat.ambiguous += 1
                        print_debug("Ambigous match")
        else:
            stat.uniquely_assigned += 1
            print_debug("Unique match")
            transcript_id = list(matched_isoforms)[0] 
            codon_pair = self.codon_pairs[transcript_id]

        return transcript_id, codon_pair

    # resolve assignment ambiguities caused by identical profiles
    def resolve_ambiguous(self, barcode_info, matched_isoforms, profile_storage):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
        for t in matched_isoforms:
            matched_positions = find_matching_positions(profile_storage.isoform_profiles[t], bacrode_jprofile)
            #print_debug(matched_positions)
            #print_debug(profile_storage.isoform_profiles[t])
            
            all_junctions_detected = True
            for i in range(len(matched_positions)):
                if matched_positions[i] == 0 and profile_storage.isoform_profiles[t][i] != -1:
                    all_junctions_detected = False
                    break

            if all_junctions_detected:
                #print_debug("Ambiguity resolved")
                #print_debug(profile_storage.isoform_profiles[t])
                #print_debug(matched_positions)
                #print_debug(bacrode_jprofile)
                return set([t])

        return matched_isoforms


# Class for processing entire bam file agains gene database
class GeneDBProcessor:
    bc_map = None
    db = None
    args = None
    bamfile_name = ""
    output_prefix = ""
    stats = None
    chr_bam_prefix = ""
    
    def __init__(self, args):
        self.args = args
        self.bc_map = self.get_barcode_map(args.bam)
        self.bamfile_name = args.bam
        if not os.path.isfile(self.bamfile_name):
            raise Exception("BAM file " + self.samfile_nam + " does not exist")
        samfile_in = pysam.AlignmentFile(self.bamfile_name, "rb")
        if not samfile_in.has_index:
            raise Exception("BAM file " + self.samfile_nam + " is not indexed, run samtools index")
        if samfile_in.references[0].startswith('chr'):
            print("Changing chomosome prefix")
            self.chr_bam_prefix = 'chr'
        samfile_in.close()

        if not os.path.isfile(args.genedb):
            raise Exception("Gene database " + args.genedb + " does not exist")
        self.db = gffutils.FeatureDB(args.genedb, keep_order = True)
        
        self.stats = BarcodeAssignmentStats()
        self.output_prefix = args.output_prefix

    # read barcode map file generated along with a BAM file
    def get_barcode_map(self, sam_file_name):
        barcode_map = {}
        contigs_name, ext = os.path.splitext(sam_file_name)
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


    #get barcode or sequence id depending on data type
    def get_sequence_id(self, query_name):
        if self.args.data_type == "10x":
            tokens = query_name.strip().split("___")
            if len(tokens) != 2:
                return ""
            return tokens[1]
        elif self.args.data_type == "contigs":
            return query_name.strip().split("_barcodeIDs_")[0]
        else:
            return query_name.strip()

    #add isoform stats whem mapping reference sequences
    def count_isoform_stats(self, isoform, barcode_id, gene_stats, gene_info):
        gene_isoform_ids = set(gene_info.coding_rna_profiles.isoform_profiles.keys())
        gene_all_isoform_ids = set(gene_info.all_rna_profiles.isoform_profiles.keys())

        if isoform is not None:
            if barcode_id == isoform:
                gene_stats.correctly_assigned += 1
            elif barcode_id in gene_isoform_ids:
                gene_stats.incorrectly_assigned_same_gene += 1
                print_debug("Incorrect assignment: isoform = " + isoform + " ; sequence = " + barcode_id + "\n")
            elif barcode_id in gene_all_isoform_ids:
                gene_stats.incorrectly_assigned_nc += 1
            else:
                gene_stats.incorrectly_assigned_other_gene += 1
                print_debug("Alien assignment: isoform = " + isoform + " ; sequence = " + barcode_id + "\n")

        else:
            if barcode_id in gene_isoform_ids:
                if barcode_id in gene_info.all_rna_profiles.ambiguous:
                    gene_stats.unassignable += 1
                    #print_debug("Barcode " + barcode_id + " is unassignable")
                else:
                    gene_stats.unassigned += 1
                    #print_debug("Barcode " + barcode_id + " is unassigned")
            elif barcode_id in gene_all_isoform_ids:
                if barcode_id in gene_info.all_rna_profiles.ambiguous:
                    gene_stats.unassignable += 1
                    #print_debug("Barcode " + barcode_id + " is unassignable")
                else:
                    gene_stats.unassigned_nc += 1
                    #print_debug("Barcode " + barcode_id + " is unassigned")
            elif barcode_id in gene_info.all_rna_profiles.empty:
                gene_stats.empty_bc += 1
            else:
                gene_stats.mismapped += 1


    def count_unmapped_stats(self, gene_info,  gene_stats):
        gene_isoform_ids = set(gene_info.coding_rna_profiles.isoform_profiles.keys())
        processed_barcodes = set()
        for t in gene_info.barcodes.keys():
            processed_barcodes.add(self.bc_map[t][0])
        for t in gene_isoform_ids:
            if t not in processed_barcodes:
                gene_stats.unmapped += 1


    #assign all barcodes mapped to gene region
    def get_gene_barcodes(self, gene_info):
        samfile_in = pysam.AlignmentFile(self.bamfile_name, "rb")
        gene_chr, gene_start, gene_end = gene_info.get_gene_region()

        #process all alignments
        #prefix is needed when bam file has chrXX chromosome names, but reference has XX names
        for alignment in samfile_in.fetch(self.chr_bam_prefix + gene_chr, gene_start, gene_end):
            if alignment.reference_id == -1:
                continue
            seq_id = self.get_sequence_id(alignment.query_name)
            gene_info.add_read(alignment, seq_id)
        samfile_in.close()

        barcodes = {}
        gene_stats = BarcodeAssignmentStats()

        #iterate over all barcodes / sequences and assign them to known isoforms
        for b in gene_info.barcodes.keys():
            isoform, codons = gene_info.assign_isoform(b, gene_stats, self.args.reads_cutoff)

            barcode_id = b
            if self.bc_map is not None:
                barcode_id = list(self.bc_map[b])[0]
                for bc in self.bc_map[b]:
                    barcodes[bc] = (isoform, codons)
            else:
                barcodes[b] = (isoform, codons)
            
            if self.args.count_isoform_stats:
                self.count_isoform_stats(isoform, barcode_id, gene_stats, gene_info)
        
        if self.args.count_isoform_stats:
            self.count_unmapped_stats(gene_info,  gene_stats)

        print("Done. Barcodes stats " + gene_stats.to_str())
        if self.args.count_isoform_stats:
            print("Done. Isoform stats " + gene_stats.isoform_stats())
        self.stats.merge(gene_stats)

        return barcodes


    def gene_list_id_str(self, gene_db_list, delim = "_"):
        gene_names = [g.id for g in gene_db_list]
        return delim.join(gene_names)
    

    def write_gene_stats(self, gene_db_list, barcodes):
        #writing TSV with barcode -> isoform id
        outf = open(self.out_tsv, "a+")
        outf.write(self.gene_list_id_str(gene_db_list) + "\t" + str(len(barcodes)) + "\n")
        for b in barcodes.keys():
            if b != "" and barcodes[b][0] is not None:
                outf.write(b + "\t" + barcodes[b][0] + "\n")
        outf.close()

    def write_codon_tables(self, gene_db_list, gene_info, barcodes):
        if CODON_OUTPUT == "ignore_overlaps" and len(gene_db_list) == 1:
            self.write_codon_tables_for_genes(gene_db_list, gene_info, barcodes)
        elif CODON_OUTPUT == "merge":
            self.write_codon_tables_for_genes(gene_db_list, gene_info, barcodes)
        elif CODON_OUTPUT == "separate":
            for g in gene_db_list:
                self.write_codon_tables_for_genes([g], gene_info, barcodes)
        elif CODON_OUTPUT == "merge_exons":
            for genes in self.group_genes_with_overlapping_exons(gene_db_list):
                self.write_codon_tables_for_genes(genes, gene_info, barcodes)
        else:
            print("Unsupported codon output method")


    def group_genes_with_overlapping_exons(self, gene_db_list):
        sys.exit(-1)
        return []

    def write_codon_tables_for_genes(self, gene_db_list, gene_info, barcodes):
        #writing codon stats
        gene_codon_pairs =  set(gene_info.get_codon_pairs(gene_db_list).values())
        codon_count_table = {}
        start_codons = set()
        stop_codons = set()
        for b in barcodes.keys():
            if b != "" and barcodes[b][0] is None:
                continue

            codon_pair = barcodes[b][1]
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


    # Process a set of genes given in gene_db_list
    def process_gene_list(self, gene_db_list):
        print("Processing " + str(len(gene_db_list)) + " gene(s): " + self.gene_list_id_str(gene_db_list, ", "))
        
        gene_info = GeneBarcodeInfo(gene_db_list, self.db, self.args, self.chr_bam_prefix)
        barcodes = self.get_gene_barcodes(gene_info)

        self.write_gene_stats(gene_db_list, barcodes)
        self.write_codon_tables(gene_db_list, gene_info, barcodes)

    #Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_all_genes(self):
        self.out_tsv = self.output_prefix + ".barcodes.tsv"
        outf = open(self.out_tsv, "w")
        outf.close()
        self.out_codon_stats = self.output_prefix + ".codon_stats.tsv"
        outf = open(self.out_codon_stats, "w")
        outf.close()
        
        gene_db_list = []

        for g in self.db.features_of_type('gene', order_by=('seqid', 'start')):
            gene_name = g.id
            gene_db = self.db[gene_name]

            if len(gene_db_list) == 0 or any(genes_overlap(g, gene_db) for g in gene_db_list):
                gene_db_list.append(gene_db)
            else:
                self.process_gene_list(gene_db_list)
                gene_db_list = [gene_db]

        self.process_gene_list(gene_db_list)

        print("Finished. Total stats " + self.stats.to_str() + "\n")
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

    for k in global_barcode_map.keys():
        b = k
        if bc_map is not None:
            b = bc_map[k][0]
        matched_isoforms = global_barcode_map[k]
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

    print("\nGlobal stats, total barcodes processd " + str(len(global_barcode_map)))
    print("\t\tCorrect\tIncorrect")
    print("Unique\t\t" + str(correct_unique) + "\t" + str(incorrect_unique))
    print("Ambiguos\t" + str(correct_amb) + "\t" + str(incorrect_amb))
    print("Unassinable\t" + str(correct_unassignable) + "\t" + str(incorrect_unassignable))
    print("\nNot matched " + str(unassigned) + ", unassignable and not matched " + str(unassigned_unassignable))
        

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', metavar='BAM_FILE', type=str,  help='sorted and indexed BAM file')
    parser.add_argument("--data_type", "-d", help="type of data to process, supported types are: contigs, 10x, long_reads, isoforms", type=str, default = "10x")
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    args = parser.parse_args()
    return args

#Tune algorithm params
def set_params(args):
    args.reads_cutoff = READS_CUTOFF if args.data_type == "10x" else 0
    args.assign_codons_when_ambiguous = ASSIGN_CODONS_WHEN_AMBIGUOUS and args.data_type != "isoforms"
    args.consider_flanking_junctions = CONSIDER_FLANKING_JUNCTIONS and args.data_type != "10x" 
    args.junction_delta = LR_JUNCTION_DELTA  if args.data_type == "long_reads" else JUNCTION_DELTA
    args.count_isoform_stats = COUNT_ISOFORM_STATS and args.data_type == "isoforms"


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

