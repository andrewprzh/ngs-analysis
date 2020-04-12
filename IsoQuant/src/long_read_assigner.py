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


class ReadAssignment:
    def __init__(self, read_id, assignment_type = [], assigned_features = []):
        self.read_id = read_id
        if isinstance(assignment_type, list):
            self.assignment_type = assignment_type
        else:
            self.assignment_type = [assignment_type]
        self.assigned_features = assigned_features

    def add_feature(self, feature, type):
        self.assigned_features.append(feature)
        self.assignment_type.append(type)


class LongReadAssigner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params

    # === Basic functions ===
    # match read profiles to a known isoform junction profile, hint - potential candidates
    def match_profile(self, read_gene_profile, isoform_profiles, hint = set()):
        isoforms = []
        for isoform_id in isoform_profiles.keys():
            if len(hint) > 0 and isoform_id not in hint:
                continue
            isoform_profile = isoform_profiles[isoform_id]
            diff = difference_in_present_features(isoform_profile, read_gene_profile)
            isoforms.append((isoform_id, diff))
        return sorted(isoforms, key=lambda x:x[1])

    # match read profiles to a known isoform junction profile, hint - potential candidates
    def find_matching_isofoms(self, read_gene_profile, isoform_profiles, hint = set()):
        isoforms = self.match_profile(read_gene_profile, isoform_profiles, hint)
        return set(filter(lambda x: x[1] == 0, isoforms))

    # === Isoforom matching function ===
    def assign_to_isoform(self, read_id, read_intron_profile, read_split_exon_profile):
        if all(el == 0 for el in read_split_exon_profile.read_profile):
            # none of the blocks matched
            return ReadAssignment(read_id, "empty")
        elif any(el == 0 for el in read_intron_profile.read_profile) or any(el == 0 for el in read_split_exon_profile.read_profile):
            # Read has extra flanking exons / introns
            return self.resolve_extra_flanking(read_id, read_intron_profile, read_split_exon_profile)
        elif any(el == -1 for el in read_intron_profile.read_profile) or any(el == -1 for el in read_split_exon_profile.read_profile):
            # Read has missing exons / introns
            return self.resolve_contradictory(read_id, read_intron_profile, read_split_exon_profile)
        else:
            assignment = self.match_non_contradictory(read_id, read_intron_profile, read_split_exon_profile)
            #check for extra flanking sequences
            if not self.params.ingnore_extra_flanking:
                isoform_exon_profile = self.gene_info.split_exon_profiles.profiles[assignment.assigned_features[0]]
                start_exon = self.gene_info.split_exon_profiles.features[isoform_exon_profile.index[1]]
                end_exon = self.gene_info.split_exon_profiles.features[isoform_exon_profile.rindex[1]]
                extra5 = read_split_exon_profile.read_features[0][0] + self.params.delta < start_exon[0]
                extra3 = read_split_exon_profile.read_features[-1][1] - self.params.delta > end_exon[0]
                if extra3 and extra5:
                    assignment.assignment_type[0] += "_extra_seq53"
                elif extra3:
                    assignment.assignment_type[0] += "_extra_seq3"
                elif extra5:
                    assignment.assignment_type[0] += "_extra_seq5"

    # match profile when all read features are assigned
    def match_non_contradictory(self, read_id, read_intron_profile, read_split_exon_profile):
        intron_matched_isoforms = self.find_matching_isofoms(read_intron_profile.gene_profile,
                                                                        self.gene_info.intron_profiles.profiles)
        exon_matched_isoforms = self.find_matching_isofoms(read_split_exon_profile.gene_profile,
                                                                      self.gene_info.split_exon_profiles.profiles)
        both_mathched_isoforms = intron_matched_isoforms.intersection(exon_matched_isoforms)

        read_assignment = None
        if len(both_mathched_isoforms) == 1:
            read_assignment = ReadAssignment(read_id, "unique", list(both_mathched_isoforms))

        elif len(both_mathched_isoforms) > 1:
            if not params.resolve_ambiguous or len(self.gene_info.ambiguous_isoforms.intersection(both_mathched_isoforms)) > 0:
                read_assignment = ReadAssignment(read_id, "ambiguous", list(both_mathched_isoforms))
            else:
                possible_isoforms = self.gene_info.ambiguous_isoforms.intersection(both_mathched_isoforms)
                intron_matched_isoforms = \
                    self.resolve_ambiguous(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles, possible_isoforms)
                exon_matched_isoforms = \
                    self.resolve_ambiguous(read_split_exon_profile.gene_profile, self.gene_info.split_exon_profiles.profile, possible_isoforms)

                both_mathched_isoforms = intron_matched_isoforms.intersection(exon_matched_isoforms)

                if len(both_mathched_isoforms) == 1:
                    read_assignment = ReadAssignment(read_id, "ambiguous_isoform", list(both_mathched_isoforms))
                else:
                    read_assignment = ReadAssignment(read_id, "ambiguous", list(both_mathched_isoforms))

        elif len(both_mathched_isoforms) == 0:
            if len(intron_matched_isoforms) > 0:
                if len(intron_matched_isoforms) == 1:
                    read_assignment = ReadAssignment(read_id, "unique_extra_exon", list(intron_matched_isoforms))
                else:
                    read_assignment = ReadAssignment(read_id, "ambiguous_extra_exon", list(intron_matched_isoforms))
            elif len(exon_matched_isoforms) > 0:
                if len(exon_matched_isoforms) == 1:
                    read_assignment = ReadAssignment(read_id, "unique_extra_intron_wtf", list(exon_matched_isoforms))
                else:
                    read_assignment = ReadAssignment(read_id, "ambiguous_extra_intron_wtf", list(exon_matched_isoforms))
            else:
                # alternative isoforms made of known introns/exons or intron retention
                return self.resolve_unmatched(read_id, read_intron_profile, read_split_exon_profile)
        return read_assignment

    # resolve assignment ambiguities caused by identical profiles
    def resolve_ambiguous(self, read_gene_profile, isoform_profiles, matched_isoforms):
        for t in matched_isoforms:
            matched_positions = find_matching_positions(isoform_profiles[t], read_gene_profile)

            all_positions_detected = True
            for i in range(len(matched_positions)):
                if matched_positions[i] == 0 and isoform_profiles[t][i] == 1:
                    all_positions_detected = False
                    break

            if all_positions_detected:
                return set([t])

        return matched_isoforms

    #resolve when there are 0s  at the ends of read profile
    def resolve_extra_flanking(self, read_id, read_intron_profile, read_split_exon_profile):
        if not self.params.ingnore_extra_flanking:
            assignment = ReadAssignment(read_id, "contradictory")
        else:
            assignment = self.match_non_contradictory(read_id, read_intron_profile, read_split_exon_profile)

        if read_intron_profile.read_profile[0] != 1 and read_intron_profile.read_profile[-1] != 1:
            assignment_type_suffix = "_extra_intron35"
        elif read_intron_profile.read_profile[0] != 1:
            assignment_type_suffix = "_extra_intron5"
        elif read_intron_profile.read_profile[-1] != 1:
            assignment_type_suffix = "_extra_intron3"
        elif read_intron_profile.read_profile[0] == 1 and read_intron_profile.read_profile[-1] == 1:
            assignment_type_suffix = "_extra_exon_wtf"
        else:
            assignment_type_suffix = "_extra_intron_wtf"
        assignment.assignment_type[0] += assignment_type_suffix
        return assignment

    #resolve when there are no exactly matching isoforms, but no -1s in read profiles
    def resolve_unmatched(self, read_id, read_intron_profile, read_split_exon_profile):
        #FIXME do not exec match profile for the second time
        intron_matching_isoforms = self.match_profile(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles)
        exon_matching_isoforms = self.match_profile(read_split_exon_profile.gene_profile, self.gene_info.split_exon_profiles.profiles)

        assignment = ReadAssignment(read_id)
        if intron_matching_isoforms[0][1] == 0:
            # difference is only in exons
            isoform_ids = get_first_best_from_sorted(exon_matching_isoforms)
            assignment_type = "ambiguous" if len(isoform_ids) > 1 else "unique"

            for matched_isoform in isoform_ids:
                isoform_exon_comparative_profile = mask_profile(self.gene_info.split_exon_profiles.profiles[isoform_id],
                                                                read_split_exon_profile.gene_profile)

                if isoform_exon_comparative_profile[0] != 1 and isoform_exon_comparative_profile[-1] != 1:
                    assignment_type_suffix = "_extra_known_exon35"
                elif isoform_exon_comparative_profile[0] != 1:
                    assignment_type_suffix = "_extra_known_exon5"
                elif isoform_exon_comparative_profile[-1] != 1:
                    assignment_type_suffix = "_extra_known_exon3"
                else:
                    assignment_type_suffix = "_extra_known_exon_wtf"
                assignment.add_feature(matched_isoform, assignment_type + assignment_type_suffix)

                return assignment

        # difference in introns
        return self.detect_differences(read_id, read_intron_profile, read_split_exon_profile, intron_matching_isoforms,
                                       exon_matching_isoforms)

    #resolve when there are -1s in read profile
    def resolve_contradictory(self, read_id, read_intron_profile, read_split_exon_profile):
        intron_matching_isoforms = self.match_profile(read_intron_profile.gene_profile, self.gene_info.intron_profiles.profiles)
        exon_matching_isoforms = self.match_profile(read_split_exon_profile.gene_profile,
                                           self.gene_info.split_exon_profiles.profiles)
        return self.detect_differences(read_id, read_intron_profile, read_split_exon_profile, intron_matching_isoforms, exon_matching_isoforms)

    # compare read to closest matching isoforms
    def detect_differences(self, read_id, read_intron_profile, read_split_exon_profile, intron_matching_isoforms, exon_matching_isoforms):
        # get isoforms that have closes intron and exon profiles
        best_intron_isoform_ids = get_first_best_from_sorted(intron_matching_isoforms)
        best_exon_isoforms = list(filter(lambda x: x[0] in best_intron_isoform_ids, exon_matching_isoforms))
        best_isoform_ids = get_first_best_from_sorted(best_exon_isoforms)

        assignment = ReadAssignment(read_id)
        for isoform_id in best_isoform_ids:
            # get intron coordinates
            isoform_introns = get_blocks_from_profile(self.gene_info.intron_profiles.features,
                                                      self.gene_info.intron_profiles.profiles[isoform_id])
            # read start-end coordinates
            read_region = (read_split_exon_profile.read_features[0][0], read_split_exon_profile.read_features[0][0])
            # isoform start-end
            isoform_exon_profile = self.gene_info.split_exon_profiles.profiles[isoform_id]
            isoform_start = self.gene_info.split_exon_profiles.features[isoform_exon_profile.index(1)]
            isoform_end = self.gene_info.split_exon_profiles.features[isoform_exon_profile.rindex(1)]
            isoform_region = (isoform_start, isoform_end)
            assignment_type = self.compare_junctions(read_intron_profile.read_features, read_region,
                                                     isoform_introns, isoform_region)
            assignment.add_feature(isoform_id, assignment_type)
        return assignment

    # compare read splice junctions against similar isoform
    def compare_junctions(self, read_junctions, read_region, isoform_junctions, isoform_region,):
        read_pos = 0
        isoform_pos = 0
        read_features_present = [0 for i in range(0, len(read_junctions))]
        isoform_features_present = [0 for i in range(0, len(isoform_junctions))]
        contradictory_region_pairs = []
        current_contradictory_region = (None, None)

        while read_pos < len(read_junctions) and isoform_pos < len(isoform_junctions):
            if equal_ranges(isoform_junctions[isoform_pos], read_junctions[read_pos], self.params.delta):
                # junctions are equal
                read_features_present[read_pos] = 1
                isoform_features_present[isoform_pos] = 1
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                read_pos += 1
                isoform_pos += 1

            elif overlaps(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # junctions overlap, but are unequal
                read_features_present[read_pos] = -1
                isoform_features_present[isoform_pos] = -1
                if current_contradictory_region == (None, None):
                    current_contradictory_region = ((read_pos, read_pos), (isoform_pos, isoform_pos))
                else:
                    current_contradictory_region = (
                    (current_contradictory_region[0][0], read_pos), (current_contradictory_region[1][0], isoform_pos))
                if (read_junctions[read_pos][1] < isoform_junctions[isoform_pos][1]):
                    read_pos += 1
                else:
                    isoform_pos += 1

            elif left_of(isoform_junctions[isoform_pos], read_junctions[read_pos]):
                # isoform junction is behing, move on
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                if read_pos > 0 or overlaps(read_region, isoform_junctions[isoform_pos]):
                    if (isoform_features_present[isoform_pos] != -1):
                        contradictory_region_pairs.append((None, (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
                isoform_pos += 1

            else:
                # read junction is behing, move on
                if (current_contradictory_region != (None, None)):
                    contradictory_region_pairs.append(current_contradictory_region)
                    current_contradictory_region = (None, None)
                if isoform_pos > 0 or overlaps(isoform_region, read_junctions[read_pos]):
                    if (read_features_present[read_pos] != -1):
                        contradictory_region_pairs.append(((read_pos, read_pos), None))
                    read_features_present[read_pos] = -1
                read_pos += 1

        if (current_contradictory_region != (None, None)):
            contradictory_region_pairs.append(current_contradictory_region)

        # check terminating regions
        while read_pos < len(read_junctions):
            if overlaps(isoform_region, read_junctions[read_pos]):
                if (read_features_present[read_pos] != -1):
                    contradictory_region_pairs.append(((read_pos, read_pos), None))
                    read_features_present[read_pos] = -1
            else:
                break
            read_pos += 1

        while isoform_pos < len(isoform_junctions):
            if overlaps(read_region, isoform_junctions[isoform_pos]):
                if (isoform_features_present[isoform_pos] != -1):
                    contradictory_region_pairs.append((None, (isoform_pos, isoform_pos)))
                    isoform_features_present[isoform_pos] = -1
            else:
                break
            isoform_pos += 1

        if any(el == -1 for el in read_features_present) or any(el == -1 for el in isoform_features_present):
            # classify contradictions
            return self.detect_contradiction_type(read_junctions, isoform_junctions, contradictory_region_pairs)
        else:
            logger.warn("No contradition detected")
            return "undefined"

    def detect_contradiction_type(self, read_junctions, isoform_junctions, contradictory_region_pairs):
        contradiction_events = []
        for pair in contradictory_region_pairs:
            # classify each contradictory area separately
            event = self.compare_overlapping_contradictional_regions(read_junctions, isoform_junctions, pair[0], pair[1])
            contradiction_events.append(event)

        contradiction_events_set = set(contradiction_events)
        if len(contradiction_events_set) == 1:
            return contradiction_events[0]
        else:
            return ",".join(list(contradiction_events_set))

    def compare_overlapping_contradictional_regions(self, read_junctions, isoform_junctions, read_cregion, isoform_cregion):
        if read_cregion is None:
            return "intron_retention"
        elif isoform_cregion is None:
            if self.are_known_introns(read_junctions, read_cregion):
                return "additional_known_intron"
            return "additional_unknown_intron"

        read_intron_total_len = sum(
            [read_junctions[i][1] - read_junctions[i][0] for i in range(read_cregion[0], read_cregion[1] + 1)])
        isoform_intron_total_len = sum(
            [isoform_junctions[i][1] - isoform_junctions[i][0] for i in range(isoform_cregion[0], isoform_cregion[1] + 1)])
        total_intron_len_diff = abs(read_intron_total_len - isoform_intron_total_len)

        read_introns_known = self.are_known_introns(read_junctions, read_cregion)

        if read_cregion[1] == read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if total_intron_len_diff < self.params.delta:
                if read_introns_known:
                    return "intron_migration"
                else:
                    return "intron_shift"
            else:
                if read_introns_known:
                    return "intron_altered"
                else:
                    return "intron_altered_novel"

        elif read_cregion[1] - read_cregion[0] == isoform_cregion[1] - isoform_cregion[0] and \
                total_intron_len_diff < self.params.delta:
            if read_introns_known:
                return "mutual_known_exons"
            else:
                return "mutual_novel_exons"

        elif read_cregion[1] == read_cregion[0] and isoform_cregion[1] > isoform_cregion[0]:
            if total_intron_len_diff < self.params.delta:
                return "small_exon_misalignment"
            else:
                if read_introns_known:
                    return "exon_skipping_known_intron"
                else:
                    return "exon_skipping_novel_intron"

        elif read_cregion[1] > read_cregion[0] and isoform_cregion[1] == isoform_cregion[0]:
            if read_introns_known:
                return "known_exon_gain"
            else:
                return "novel_exon_gain"

        else:
            if read_introns_known:
                return "alternative_known_exons"
            else:
                return "alternative_novel_exons"

    def are_known_introns(self, junctions, region):
        selected_junctions = []
        for i in range(region[0], region[1] + 1):
            selected_junctions.append(junctions[i])
        intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features,
                                                  (self.gene_info.start, self.gene_info.end),
                                                  comparator = partial(equal_ranges, delta = self.params.delta))
        selected_junctions_profile = intron_profile_constructor.construct_profile(selected_junctions)
        return all(el == 1 for el in selected_junctions_profile)

    # === Exon matching ==
    def assign_exons(self, read_exon_profile):
        pass