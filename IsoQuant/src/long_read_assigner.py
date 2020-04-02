############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from enum import *
from common import *
from gene_info import *
from long_read_profiles import *


logger = logging.getLogger('LRAssignment')


class AssignmentType(Enum):
    UNIQUE = 1


class ReadAssignment:
    def __init__(self, read_id):
        self.read_id = read_id
        self.main_assigned_feature = None
        self.main_assignment_type = None
        self.assigned_features = {}


class LongReadAssigner:
    def __init__(self, gene_info):
        self.gene_info = gene_info

    def assign_to_isoform(self, read_intron_profile, read_split_exon_profile):
        if any(el == -1 for el in read_intron_profile.read_profile):
            return "contradictory"
        elif all(el == 0 for el in read_split_exon_profile.read_profile):
            return "empty"
        elif read_split_exon_profile.read_profile[0] == 0 or read_split_exon_profile.read_profile[-1] == 0:
            return "extra"

    def assign_exons(self, read_exon_profile):
        pass

    # match read profiles to a known isoform junction profile, hint - potential candidates
    def find_matching_isofoms_by_profile(self, read_profile, profiles, hint=set()):
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



