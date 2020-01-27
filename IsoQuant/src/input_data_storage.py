############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging

logger = logging.getLogger('InputData')


class InputDataStorage:
    input_type = ""
    samples = []
    labels = []

    def __init__(self, args):
        self.samples = []
        self.labels = []
        self.input_type = ""
        if args.fastq is not None:
            self.input_type = "fastq"
            for fq in args.fastq:
                samples.append([[fq]])
        elif args.bam is not None:
            self.input_type = "bam"
            for bam in args.bam:
                samples.append([[bam]])
        elif args.fastq_list is not None:
            self.input_type = "fastq"
            self.set_samples_from_file(args.fastq_list)
        elif args.bam_list is not None:
            self.input_type = "bam"
            self.set_samples_from_file(args.bam_list)
        else:
            logger.critical("Input data was not specified")
            exit(-1)

        if args.labels is not None:
            if len(args.labels) != len(self.samples):
                logger.critical("Number of labels is not equal to the numbe of samples")
                exit(-1)
            else:
                self.labels = args.labels
        else:
            self.set_labels()

    def set_samples_from_file(self, file_name):
        inf = open(file_name, "r")
        current_sample = []

        for l in inf.readline():
            if len(l.strip()) == 0 and len(current_sample) > 0:
                self.samples.append(current_sample)
                current_sample = []
            else:
                current_sample.append(l.strip().split())

        if len(current_sample) > 0:
            self.samples.append(current_sample)

    def set_labels(self):
        # TODO
        self.labels = []
        for i in range(self.samples):
            self.labels = str(i) + "_" + os.path.splitext(os.path.basename(self.samples[i][0][0]))[0]
