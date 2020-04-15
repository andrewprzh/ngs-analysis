############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam

from src.long_read_assigner import *
from src.long_read_profiles import *

logger = logging.getLogger('IsoQuant')

# class for aggregating all assignment information
class LongReadAlginmentProcessor:
    def __init__(self, gene_info, bams, params, printer):
        self.gene_info = gene_info
        self.bams = bams
        self.printer = printer
        self.params = params

        gene_region = (gene_info.start, gene_info.end)
        self.assigner = LongReadAssigner(self.gene_info, self.params)
        self.intron_profile_construnctor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features, gene_region,
                                                  comparator = partial(equal_ranges, delta = self.params.delta))
        # TODO check for non split exons which do overlap
        self.exon_profile_construnctor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.exon_profiles.features, gene_region,
                                                  comparator = partial(equal_ranges, delta = self.params.delta))
        # TODO think whether overlaps should be changed to contains to avoid terminal partially covered exons
        self.split_exon_profile_construnctor = \
            NonOverlappingFeaturesProfileConstructor(self.gene_info.split_exon_profiles.features,
                                                     comparator=partial(overlaps_at_least, delta=self.params.delta))

    def process(self):
        for b in self.bams:
            self.process_single_file(b)

    def process_single_file(self, bam):
        bamfile_in = pysam.AlignmentFile(bam, "rb")

        for alignment in bamfile_in.fetch(self.gene_info.chr_id, self.gene_info.start, self.gene_info.end):
            if alignment.reference_id == -1:
                continue
            if self.params.skip_secondary and (alignment.is_secondary or alignment.is_supplementary):
                continue

            read_id = alignment.query_name
            logger.debug("=== Processing read " + read_id + " ===")
            concat_blocks = concat_gapless_blocks(sorted(alignment.get_blocks()), alignment.cigartuples)
            sorted_blocks = correct_bam_coords(concat_blocks)
            intron_profile = self.intron_profile_construnctor.construct_profile(sorted_blocks)
            split_exon_profile = self.split_exon_profile_construnctor.construct_profile(sorted_blocks)
            # FIXME
            combined_profile = CombinedReadProfiles(intron_profile, None, split_exon_profile)
            read_assignment = self.assigner.assign_to_isoform(read_id, combined_profile)
            logger.debug("=== Finished read " + read_id + " ===")
            self.printer.add_read_info(read_assignment, combined_profile)

        bamfile_in.close()



