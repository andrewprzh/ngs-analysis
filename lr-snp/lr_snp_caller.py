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

NUCL_COUNT = 4
NUCL_TO_INDEX_MAP = {'A': 0, 'C': 1, 'G':2, 'T':3}
NUCL_TO_INDEX_MAP_FULL = {'A': 0, 'C': 1, 'G':2, 'T':3, 'a': 0, 'c': 1, 'g':2, 't':3}
INDEX_TO_NUCL = ['A', 'C', 'G', 'T']
NUCL_TO_INDEX = [-1 if chr(x) not in NUCL_TO_INDEX_MAP_FULL.keys() else NUCL_TO_INDEX_MAP_FULL[chr(x)] for x in range(0, 256)]

GERMLINE_SNP = "germline"
SOMATIC_SNP = "somatic"
BOUND_CASE_SNP = "unidentified"


class NuclStorage:
    def __init__(self, region_size):
        self.nucl_counts = []
        for i in range(NUCL_COUNT):
            self.nucl_counts.append([0 for j in range(region_size)])

    def clean(self, region_size):
        if region_size > len(self.nucl_counts[0]):
            self.__init__(region_size)
        else:
            for i in range(NUCL_COUNT):
                for j in range(region_size):
                    self.nucl_counts[i][j] = 0

    def increment(self, pos, nucl):
        self.nucl_counts[NUCL_TO_INDEX[ord(nucl)]][pos] += 1

    def total_cov(self, pos):
        return sum([self.nucl_counts[i][pos] for i in range(NUCL_COUNT)])

    def count(self, pos, nucl):
        return self.nucl_counts[NUCL_TO_INDEX[ord(nucl)]][pos]

    def freq(self, pos, nucl):
        return float(self.nucl_counts[NUCL_TO_INDEX[ord(nucl)]][pos]) / self.total_cov(pos)

    def non_ref_freqs(self, pos, ref_nucl):
        return [self.freq(pos, nucl) if nucl != ref_nucl else 0.0 for nucl in INDEX_TO_NUCL]

    def get_counts(self, pos):
        return [self.nucl_counts[i][pos] for i in range(NUCL_COUNT)]


class SelectedSNP:
    def __init__(self, ref_nucl, alt_nucl, snp_type, sample_cov, sample_counts):
        self.reference_nucl = ref_nucl
        self.alternative_nucl = alt_nucl
        self.snp_type = snp_type
        self.sample_coverage = sample_cov
        self.sample_counts = sample_counts


class SNPCaller:
    def __init__(self, args):
        self.bamfiles = args.bam_file
        self.reference_path = args.reference
        self.snp_map = {}
        self.count_storages = [NuclStorage(1) for i in range(len(self.bamfiles))]
        self.args = args

    # find all different nucleotides in chromosome region [start, end]
    def process_region_for_bam(self, index, bam, chromosome_record, start, end):
        count_storage = self.count_storages[index]
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
                    count_storage.increment(pos - start, read_seq[read_pos])
                    read_pos += 1
                if reached_end:
                    break

    def nucl_list(self, all_region_nucls, pos):
        res = []
        for i in range(len(all_region_nucls)):
            res.append(all_region_nucls[i][pos])
        return res

    # compare nucleotide data from multiple samples
    def select_snp(self, ref_nucl, pos):
        for i in range(len(self.count_storages)):
            if self.count_storages[i].total_cov(pos) < self.args.min_cov:
                return []

        all_frequences = [count_storage.non_ref_freqs(pos, ref_nucl) for count_storage in self.count_storages]
        resulting_snps = []
        for nucl in INDEX_TO_NUCL:
            if nucl == ref_nucl:
                continue

            nucl_index = NUCL_TO_INDEX[ord(nucl)]
            nucl_per_sample_freqs = [freq[nucl_index] for freq in all_frequences]
            max_freq = max(nucl_per_sample_freqs)
            if max_freq < self.args.min_freq:
                continue

            min_freq = min(nucl_per_sample_freqs)
            if min_freq >= self.args.min_freq or \
                    (max_freq >= self.args.min_freq and min_freq > 0 and max_freq / min_freq <= self.args.min_freq_factor):
                # all snps has good frequency
                if not self.args.keep_germline:
                    continue
                resulting_snps.append(SelectedSNP(ref_nucl, nucl, GERMLINE_SNP,
                                                  [count_storage.total_cov(pos) for count_storage in self.count_storages],
                                                  [count_storage.count(pos, nucl) for count_storage in self.count_storages]))
            else:
                # some SNPs have frequency significantly different from others
                resulting_snps.append(SelectedSNP(ref_nucl, nucl, SOMATIC_SNP,
                                                  [count_storage.total_cov(pos) for count_storage in self.count_storages],
                                                  [count_storage.count(pos, nucl) for count_storage in self.count_storages]))
        return resulting_snps

    # find all SNPs in all samples in chromosome region [start, end]
    def process_bams_in_region(self, chromosome_record, start, end):
        for i in range(len(self.bamfiles)):
            bam = pysam.AlignmentFile(self.bamfiles[i], "rb")
            self.count_storages[i].clean(end - start + 1)
            self.process_region_for_bam(i, bam, chromosome_record, start, end)

        for pos in range(start, end + 1):
            multi_snp = self.select_snp(chromosome_record.seq[pos], pos - start)
            if multi_snp is not None and len(multi_snp) > 0:
                self.snp_map[chromosome_record.id][pos] = multi_snp

    # process all alignments
    def process(self, writer):
        for record in SeqIO.parse(self.reference_path, "fasta"):
            region_start = 100000000
            sys.stderr.write("\nProcessing chormosome " + record.id + "\n")
            self.snp_map[record.id] = {}
            while region_start < len(record.seq) and region_start < 105000000:
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

    def form_line(self, chromosome, pos, snp):
        l = self.delim.join([chromosome, str(pos), snp.reference_nucl, snp.alternative_nucl, snp.snp_type])
        for i in range(len(snp.sample_counts)):
            l += self.delim + self.delim.join(map(str, [snp.sample_coverage[i], snp.sample_counts[i]]))
            l += self.delim + '{:.2f}'.format(float(snp.sample_counts[i]) / float(snp.sample_coverage[i]))
        return l + '\n'

    def dump_to_file(self, snp_map):
        somatic_file = open(self.somatic_file_name, 'a+')
        germline_file = open(self.gemline_file_name, 'a+')

        for chromosome in snp_map.keys():
            #print("Processing chromosome " + chromosome)
            snps = snp_map[chromosome]
            for pos in sorted(snps.keys()):
                multi_snp = snps[pos]
                for snp in multi_snp:
                    l = self.form_line(chromosome, pos, snp)
                    if snp.snp_type == SOMATIC_SNP:
                        somatic_file.write(l)
                    else:
                        germline_file.write(l)
        somatic_file.close()
        germline_file.close()

