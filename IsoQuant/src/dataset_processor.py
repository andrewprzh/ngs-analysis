############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging
import argparse
from traceback import print_exc

import gffutils
import pysam
from Bio import SeqIO


# Class for processing all samples against gene database
class DatasetProcessor:
    db = None
    args = None
    input_dataset = None
    output_dir = ""
    output_prefix = ""

    def __init__(self, args):
        self.args = args
        self.bc_map = self.get_barcode_map(args.bam)
        self.bamfile_name = args.bam
        if not os.path.isfile(self.bamfile_name):
            raise Exception("BAM file " + self.samfile_nam + " does not exist")
        samfile_in = pysam.AlignmentFile(self.bamfile_name, "rb")
        if not samfile_in.has_index:
            raise Exception("BAM file " + self.bamfile_name + " is not indexed, run samtools index")
        if args.change_chr_prefix and samfile_in.references[0].startswith('chr'):
            print("Changing chomosome prefix")
            self.chr_bam_prefix = 'chr'
        samfile_in.close()

        if not os.path.isfile(args.genedb):
            raise Exception("Gene database " + args.genedb + " does not exist")
        self.db = gffutils.FeatureDB(args.genedb, keep_order=True)

        self.output_prefix = args.output_prefix

    def process_all_samples(self):
        pass

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
            return query_name.strip()

    # assign all reads/barcodes mapped to gene region
    def assign_all_reads(self, read_profiles):
        samfile_in = pysam.AlignmentFile(self.bamfile_name, "rb")
        gene_chr, gene_start, gene_end = read_profiles.gene_info.get_gene_region()

        # process all alignments
        # prefix is needed when bam file has chrXX chromosome names, but reference has XX names
        for alignment in samfile_in.fetch(self.chr_bam_prefix + gene_chr, gene_start, gene_end):
            if alignment.reference_id == -1 or alignment.is_secondary:
                continue
            seq_id = self.get_sequence_id(alignment.query_name)
            read_profiles.add_read(alignment, seq_id)
        samfile_in.close()

        # read / barcode id -> (isoform, codon pair)
        assigned_reads = {}
        gene_stats = ReadAssignmentStats()

        # iterate over all barcodes / sequences and assign them to known isoforms
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

    # Process a set of genes given in gene_db_list
    def process_gene_list(self, gene_db_list):
        print("Processing " + str(len(gene_db_list)) + " gene(s): " + self.gene_list_id_str(gene_db_list, ", "))

        read_profiles = ReadProfilesInfo(gene_db_list, self.db, self.args, self.chr_bam_prefix)
        assigned_reads = self.assign_all_reads(read_profiles)

        self.write_gene_stats(gene_db_list, assigned_reads)
        self.write_codon_tables(gene_db_list, read_profiles.gene_info, assigned_reads)

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_all_genes(self):
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

        print("\nFinished. Total stats " + self.stats.to_str())
        if self.args.count_isoform_stats:
            print("Finished. Isoform stats " + self.stats.isoform_stats() + "\n")