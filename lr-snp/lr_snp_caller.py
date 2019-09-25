############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
#import gffutils
import pysam
#import vcf
from Bio import SeqIO

NUCL_TO_INDEX_MAP = {'A': 0, 'C': 1, 'G':2, 'T':3}
INDEX_TO_NUCL = ['A', 'C', 'G', 'T']
GERMLINE_SNP = "germline"
SOMATIC_SNP = "somatic"
BOUND_CASE_SNP = "unidentified"


class Nucl:
    def __init__(self, ref_nucl):
        self.reference_nucl = ref_nucl
        self.counts = [0 for i in range(len(NUCL_TO_INDEX_MAP))]

    def set(self, ref_nucl):
        self.reference_nucl = ref_nucl
        for i in range(len(self.counts)):
            self.counts[i] = 0

    def clean(self):
        self.__init__(self.reference_nucl)

    def increment(self, nucl):
        if nucl in NUCL_TO_INDEX_MAP.keys():
            self.counts[NUCL_TO_INDEX_MAP[nucl]] += 1

    def total_cov(self):
        return sum(self.counts)

    def freq(self, nucl):
        return self.counts[nucl] / float(self.total_cov())

    def max_non_ref_freq(self):
        sum = self.total_cov()
        freqs = map(lambda x: float(x) / float(sum), self.counts)
        freqs[NUCL_TO_INDEX_MAP[self.reference_nucl]] = 0.0

        self.best_nucl = 0
        self.max_freq = 0
        for i in range(len(freqs)):
            if freqs[i] > self.max_freq:
                self.best_nucl = i
                self.max_freq = freqs[i]
        return (self.max_freq, self.best_nucl)

    def to_str(self):
        nucls = []
        for i in range(len(self.counts)):
            nucls.append((INDEX_TO_NUCL[i], self.counts[i]))
        return  ', '.join(map(lambda x: x[0] + ": " + str(x[1]), nucls))


class MultiSampleSNP:
    def __init__(self, snp_type, snp_list):
        self.per_sample_snps = snp_list
        self.snp_type = snp_type
        self.reference_nucl = snp_list[0].reference_nucl
        #print(snp_type + ": " + snp_list[0].to_str())

    def to_str(self, sample_names):
        res = "TYPE: " + self.snp_type + "; REF: " + self.reference_nucl
        for i in range(len(self.per_sample_snps)):
            snp = self.per_sample_snps[i]
            res += '\n' + sample_names[i] + ": " + snp.to_str() if snp is not None else ""
        return res


class SNPCaller:
    def __init__(self, args):
        self.bamfiles = args.bam_file
        self.reference_path = args.reference
        self.snp_map = {}
        self.args = args

    # find all different nucleotides in chromosome region [start, end]
    def process_region_for_bam(self, index, bam, chromosome_record, start, end):
        region_nucls = self.all_region_nucls[index]
        for i in range(start, end):
            region_nucls[i - start].set(chromosome_record.seq[i].upper())
        #sys.stderr.write("Processing alignments\n")
        for alignment in bam.fetch(chromosome_record.id, start, end):
            if alignment.reference_id == -1:
                continue
            reached_end = False
            read_seq = alignment.query_alignment_sequence
            read_pos = 0
            for block in sorted(alignment.get_blocks()):
                for pos in range(block[0], block[1]):
                    if pos < start:
                        read_pos += 1
                        continue
                    if pos > end:
                        reached_end = True
                        break
                    if read_pos >= len(read_seq):
                        print("ERROR: read position out of bounds")
                    region_nucls[pos - start].increment(read_seq[read_pos].upper())
                    read_pos += 1

                if reached_end:
                    break

    def nucl_list(self, all_region_nucls, pos):
        res = []
        for i in range(len(all_region_nucls)):
            res.append(all_region_nucls[i][pos])
        return res

    # compare nucleotide data from multiple samples
    def select_snp(self, all_region_nucls, pos):
        for i in range(len(all_region_nucls)):
            if all_region_nucls[i][pos].total_cov() < self.args.min_cov:
                return None

        nucls = []
        for i in range(len(all_region_nucls)):
            freq_and_nucl = all_region_nucls[i][pos].max_non_ref_freq()
            nucls.append(freq_and_nucl)

        min_freq = min(map(lambda x: x[0], nucls))
        max_freq = max(map(lambda x: x[0], nucls))
        if min_freq >= self.args.min_freq or \
                (max_freq >= self.args.min_freq and min_freq > 0 and max_freq / min_freq <= self.args.min_freq_factor):
            # all snps has good frequency
            detected_variatns = set(map(lambda x: x[1], nucls))
            #print(str(all_region_nucls[0][pos].counts))
            #print(str(nucls))
            #print("Good frequency, detected variants: " + str(len(detected_variatns)) + " " + str(detected_variatns))
            if len(detected_variatns) > 1:
                # more than one possible variant exist
                return MultiSampleSNP(SOMATIC_SNP, self.nucl_list(all_region_nucls, pos))
            else:
                if not self.args.keep_germline:
                    return None
                # one possible variant
                return MultiSampleSNP(GERMLINE_SNP, self.nucl_list(all_region_nucls, pos))
        elif max_freq < self.args.min_freq:
            # all positions have low frequency
            return None
        else:
            # some SNPs have frequency significantly different from others
            return MultiSampleSNP(SOMATIC_SNP, self.nucl_list(all_region_nucls, pos))

    # find all SNPs in all samples in chromosome region [start, end]
    def process_bams_in_region(self, chromosome_record, start, end):
        for i in range(len(self.bamfiles)):
            bam = pysam.AlignmentFile(self.bamfiles[i], "rb")
            self.process_region_for_bam(i, bam, chromosome_record, start, end)

        for pos in range(0, len(self.all_region_nucls[0])):
            multi_snp = self.select_snp(self.all_region_nucls, pos)
            if multi_snp is not None:
                #print("Position " + str(start + pos + 1))
                self.snp_map[chromosome_record.id][start + pos] = multi_snp

    # process all alignments
    def process(self, writer):
        self.all_region_nucls = []
        for bamfile_name in self.bamfiles:
            self.all_region_nucls.append([Nucl('A') for pos in range(0, self.args.window_lenth)])
        for record in SeqIO.parse(self.reference_path, "fasta"):
            region_start = 0#102000000
            sys.stderr.write("\nProcessing chormosome " + record.id + "\n")
            self.snp_map[record.id] = {}
            while region_start < len(record.seq): # and region_start < 103000000:
                region_end = min(len(record.seq) - 1, region_start + self.args.window_lenth - 1)
                self.process_bams_in_region(record, region_start, region_end)
                region_start += self.args.window_lenth
                sys.stderr.write("Processed " + str(region_start) + " bases\r")
                writer.dump_to_file(self.snp_map)
                self.snp_map[record.id] = {}

    # process all alignments gene by gene
    def process_with_genes(self, gene_db):
        pass

def print_snp_map(snp_map, sample_names):
    for chromosome in snp_map.keys():
        print("Chromosome " + chromosome)
        snps = snp_map[chromosome]
        for pos in sorted(snps.keys()):
            print("Position " + str(pos + 1))
            print(snps[pos].to_str(sample_names))

class SNPMapTSVWriter:
    main_header = ['CHR', 'POS', 'REF', 'ALT', 'TYPE']
    sample_header = ['TOTAL', 'READS', 'FREQ']
    delim = '\t'

    def __init__(self, out_prefix, sample_names, args):
        self.gemline_file_name = out_prefix + "germline_SNPs.tsv"
        self.somatic_file_name = out_prefix + "somatic_SNPs.tsv"
        self.sample_names = sample_names
        self.args = args
        print("Outputting results to " + self.somatic_file_name + " and " + self.gemline_file_name)
        somatic_file = open(self.somatic_file_name, 'w')
        germline_file = open(self.gemline_file_name, 'w')
        header = self.delim.join(self.main_header)
        for sample in self.sample_names:
            for h in self.sample_header:
                header += self.delim + sample + '_' + h
        somatic_file.write(header + '\n')
        germline_file.write(header + '\n')
        somatic_file.close()
        germline_file.close()

    def form_line(self, chromosome, pos, snp, nucl, type):
        l = self.delim.join([chromosome, str(pos), snp.reference_nucl, INDEX_TO_NUCL[nucl], type])
        for sample_snp in snp.per_sample_snps:
            l += self.delim + self.delim.join(map(str, [sample_snp.total_cov(), sample_snp.counts[nucl]]))
            l += self.delim + '{:.2f}'.format(sample_snp.freq(nucl))
        return l + '\n'

    def snp_is_detected(self, snp_freq, max_freq):
        return snp_freq >= self.args.min_freq or (snp_freq > 0 and max_freq >= self.args.min_freq
                                                  and max_freq / snp_freq <= self.args.min_freq_factor)

    def detect_type(self, snp, nucl):
        freqs = map(lambda x: x.freq(nucl), snp.per_sample_snps)
        min_freq = min(freqs)
        max_freq = max(freqs)
        if min_freq > 0 and max_freq / min_freq <= self.args.min_freq_factor:
            return GERMLINE_SNP
        return SOMATIC_SNP


    def dump_to_file(self, snp_map):
        somatic_file = open(self.somatic_file_name, 'a+')
        germline_file = open(self.gemline_file_name, 'a+')

        for chromosome in snp_map.keys():
            #print("Processing chromosome " + chromosome)
            snps = snp_map[chromosome]
            for pos in sorted(snps.keys()):
                snp = snps[pos]
                detected_variants = set()
                overall_max_freq = max(map(lambda x:x.max_freq, snp.per_sample_snps))
                for sample_snp in snp.per_sample_snps:
                    if self.snp_is_detected(sample_snp.max_freq, overall_max_freq):
                        detected_variants.add(sample_snp.best_nucl)
                for nucl in detected_variants:
                    nucl_type = self.detect_type(snp, nucl)
                    if nucl_type == SOMATIC_SNP:
                        somatic_file.write(self.form_line(chromosome, pos, snp, nucl, nucl_type))
                    else:
                        germline_file.write(self.form_line(chromosome, pos, snp, nucl, nucl_type))
        somatic_file.close()
        germline_file.close()
        #print("Outputting done")
