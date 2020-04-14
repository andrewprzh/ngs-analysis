############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial

from src.common import *

logger = logging.getLogger('IsoQuant')

# The following 2 classes are very similar, but lets keep them separately for now

class MappedReadProfile:
    def __init__(self, gene_profile, read_profile, read_features):
        self.gene_profile = gene_profile
        self.read_profile = read_profile
        self.read_features = read_features


class CombinedReadProfiles:
    def __init__(self, read_intron_profile, read_split_profile, read_split_exon_profile):
        self.read_intron_profile = read_intron_profile
        self.read_exon_profile = read_split_profile
        self.read_split_exon_profile = read_split_exon_profile


#accepts sorted gapless alignment blocks
class OverlappingFeaturesProfileConstructor:
    def __init__(self, known_introns, gene_region, comparator = partial(equal_ranges, delta = 0)):
        self.known_introns = known_introns
        self.gene_region = gene_region
        self.comparator = comparator

    def construct_profile(self, sorted_blocks):
        intron_profile = [0] * (len(self.known_introns))
        if len(sorted_blocks) < 2:
            return  MappedReadProfile(intron_profile, [], [])

        read_introns = junctions_from_blocks(sorted_blocks)
        read_profile = [0] * (len(read_introns))

        mapped_region = (sorted_blocks[0][0], sorted_blocks[-1][1])
        for i in range(len(intron_profile)):
            if overlaps(self.known_introns[i], mapped_region):
                intron_profile[i] = -1
        for i in range(len(read_profile)):
            if overlaps(read_introns[i], self.gene_region):
                read_profile[i] = -1

        gene_pos = 0
        read_pos = 0
        # TODO reduce excessive if statements
        while gene_pos < len(self.known_introns) and read_pos < len(read_introns):
            if self.comparator(read_introns[read_pos], self.known_introns[gene_pos]):
                intron_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                gene_pos += 1
            elif overlaps(read_introns[read_pos], self.known_introns[gene_pos]):
                intron_profile[gene_pos] = -1
                gene_pos += 1
            elif left_of(read_introns[read_pos], self.known_introns[gene_pos]):
                if read_profile[read_pos] == 0 and gene_pos > 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            else:
                if read_pos > 0:
                    intron_profile[gene_pos] = -1
                gene_pos += 1

        return MappedReadProfile(intron_profile, read_profile, read_introns)

#accepts sorted gapless alignment blocks
class NonOverlappingFeaturesProfileConstructor:
    def __init__(self, known_exons, comparator = overlaps):
        self.known_exons = known_exons
        self.comparator = comparator

    def construct_profile(self, sorted_blocks):
        exon_profile = [0] * (len(self.known_exons))
        read_profile = [0] * (len(sorted_blocks))
        read_exons = sorted_blocks
        gene_pos = 0
        read_pos = 0

        while gene_pos < len(self.known_exons) and read_pos < len(read_exons):
            if self.comparator(read_exons[read_pos], self.known_exons[gene_pos]):
                exon_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                gene_pos += 1
                read_pos += 1
            elif left_of(read_exons[read_pos], self.known_exons[gene_pos]):
                if gene_pos > 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            else:
                if read_pos > 0:
                    exon_profile[gene_pos] = -1
                gene_pos += 1

        return MappedReadProfile(exon_profile, read_profile, read_exons)
