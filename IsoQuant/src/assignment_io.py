############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from src.common import *

logger = logging.getLogger('IsoQuant')

class PrintAllFunctor:
    def check(sefl, assignment):
        return True


class PrintStartingWithSetFunctor:
    def __init__(self, *prefix_sets):
        self.prefix_set = set()
        self.update(*prefix_sets)

    def update(self, *prefix_sets):
        for s in prefix_sets:
            self.prefix_set.update(s)

    def check(self, assignment):
        return any(assignment.assignment_type[0].startswith(prefix) for prefix in self.prefix_set)


class AbstractAssignmentPrinter:
    def __init__(self, output_file_name, params, format, assignment_checker=PrintAllFunctor()):
        self.params = params
        self.format = format
        self.assignment_checker = assignment_checker
        self.output_file = open(output_file_name, "w")

    def __del__(self):
        self.output_file.close()

    def add_read_info(self, read_assignment, combined_read_profile = None):
        raise NotImplementedError()

    def flush(self):
        self.output_file.flush()


class ReadAssignmentCompositePrinter:
    def __init__(self, printers = []):
        self.pinters = printers

    def add_read_info(self, read_assignment, mapping_read_profile=None):
        for p in self.pinters:
            p.add_read_info(read_assignment, mapping_read_profile)

    def flush(self):
        for p in self.pinters:
            p.flush()


class BasicTSVAssignmentPrinter(AbstractAssignmentPrinter):
    def __init__(self, output_file_name, params, format = "TSV", assignment_checker=PrintAllFunctor()):
        AbstractAssignmentPrinter.__init__(self, output_file_name, params, format, assignment_checker)
        self.header = "#read_id\tisoform_id\tassignment_type"
        if self.params.print_additional_info:
            self.header += "\taligned_blocks\tintron_profile\tsplit_exon_profile"
        self.header += "\n"
        self.output_file.write(self.header)

    def add_read_info(self, read_assignment, combined_read_profile = None):
        if not self.assignment_checker.check(read_assignment):
            return

        line = read_assignment.read_id + "\t" + ",".join(read_assignment.assigned_features) + "\t" + ",".join(read_assignment.assignment_type)
        if self.params.print_additional_info:
            if combined_read_profile is None:
                line += "\t.\t.\t."
            else:
                line += "\t" + range_list_to_str(combined_read_profile.read_split_exon_profile.read_features) + "\t" + \
                    list_to_str(combined_read_profile.read_intron_profile.gene_profile) + "\t" + \
                        list_to_str(combined_read_profile.read_split_exon_profile.gene_profile)
        line += "\n"
        self.output_file.write(line)
