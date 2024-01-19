#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gffutils
from traceback import print_exc
from collections import defaultdict
import glob


class GeneInfo:
    def __init__(self, gene_name, gene_len, spliced):
        self.gene_name = gene_name
        self.gene_len = gene_len
        self.spliced = spliced
        self.tpms = []


def collect_genes(db):
    gene_dict = {}
    for gene in db.features_of_type('gene'):
        transcript_lengths = []
        spliced = False
        for t in db.children(gene, featuretype=('transcript', 'mRNA')):
            tlen = 0
            exon_count = 0
            for e in db.children(t, order_by='start'):
                if e.featuretype == 'exon':
                    tlen += e.end - e.start + 1
                    exon_count += 1
            transcript_lengths.append(tlen)
            if exon_count > 1: spliced = True
        gene_dict[gene.id] = GeneInfo(gene.id, max(transcript_lengths), spliced)
    return gene_dict


def add_tpms_from_tpm_file(tpm_file, sample_num, gene_dict):
    for l in open(tpm_file):
        if l.startswith("#") or l.startswith("_"):
            continue
        v = l.strip().split("\t")
        gene_id = v[0]
        tpm = float(v[1])
        assert sample_num == len(gene_dict[gene_id].tpms)
        gene_dict[gene_id].tpms.append(tpm)


def add_tpms_from_assignment_file(assignment_file, gene_dict):
    gene_counts = defaultdict(float)
    for l in open(assignment_file):
        if l.startswith("#"):
            continue
        v = l.strip().split("\t")
        gene_id = v[4]
        assignment = v[5]
        if assignment.startswith("unique"):
            gene_counts[gene_id] += 1.0

    scale_factor = 1000000.0 / sum(gene_counts.values())
    for gene_id in gene_counts.keys():
        gene_dict[gene_id].tpms.append(gene_counts[gene_id] * scale_factor)


def get_output_dir_and_sample_name(dir, sample_name):
    dirs = []
    all_files = glob.glob(os.path.join(dir, "*"))
    for f in all_files:
        if os.path.isdir(f):
            dirs.append(os.path.basename(f))

    if len(dirs) == 0:
        print("No output was found for %s, check you output" % sample_name)

    out_dir = ""
    if len(dirs) == 1:
        out_dir = dirs[0]
    else:
        for d in dirs:
            if d.startswith("OUT"):
                out_dir = d
                break
            elif d.lower().startswith("curio"):
                out_dir = d

    if out_dir == "":
        print("No proper output dir was found for %s, check you output" % sample_name)
        out_dir = dirs[0]

    if out_dir.startswith("OUT"):
        return out_dir, sample_name
    return out_dir, out_dir


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix", required=True)
    parser.add_argument("--genedb", "-g", type=str, help="gffutils genedb used in IsoQuant", required=True)
    parser.add_argument("--isoquant_folder", "-i", type=str, help="IsoQuant folder with samples", required=True)
    parser.add_argument("--sample_list", "-s", type=str, nargs="+", help="space separated list of samples", required=True)
    parser.add_argument("--use_assignments", help="use raw read2gene assignments instead of IsoQuant TPMs", default=False, action='store_true')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gffutils_db = gffutils.FeatureDB(args.genedb)
    print("Loading gene db " + args.genedb)
    gene_dict = collect_genes(gffutils_db)

    sample_names = []
    for i, sample in enumerate(args.sample_list):
        print("Processing sample " + sample)
        out_dir = str(os.path.join(args.isoquant_folder, sample))
        sample_id, sample_name = get_output_dir_and_sample_name(out_dir, sample)
        sample_names.append(sample_name)
        sample_dir = os.path.join(out_dir, sample_id)

        if args.use_assignments:
            assignment_file = os.path.join(sample_dir, sample_id + ".read_assignments.tsv")
            add_tpms_from_assignment_file(assignment_file, gene_dict)
        else:
            tpm_file = os.path.join(sample_dir, sample_id + ".gene_tpm.tsv")
            add_tpms_from_tpm_file(tpm_file, i, gene_dict)

    print("Dumping output to " + args.output)
    all_tpm_outf = open(args.output + ".all_genes.tsv", "w")
    all_tpm_outf.write("#geneid\tlength\t" + "\t".join(sample_names) + "\n")
    spliced_outf = open(args.output + ".spliced_genes.tsv", "w")
    spliced_outf.write("#geneid\tlength\t" + "\t".join(sample_names) + "\n")

    max_bin = 5
    max_bin_cutoff = max_bin * 1000
    bin_outf_list = []
    for i in range(max_bin + 1):
        bin_outf_list.append(open(args.output + ".ceil_len_%dkb.tsv" % i, "w"))
        bin_outf_list[-1].write("#geneid\tlength\t" + "\t".join(sample_names) + "\n")

    for geneid in gene_dict.keys():
        gene_info = gene_dict[geneid]
        if gene_info.gene_len >= max_bin_cutoff:
            bin = max_bin
        else:
            bin = gene_info.gene_len // 1000
        gene_str = "%s\t%d\t" % (geneid, gene_info.gene_len) + "\t".join(map(lambda x: "%.6f" % x, gene_info.tpms)) + "\n"
        all_tpm_outf.write(gene_str)
        if gene_info.spliced:
            spliced_outf.write(gene_str)
        bin_outf_list[bin].write(gene_str)

    all_tpm_outf.close()
    spliced_outf.close()
    for f in bin_outf_list:
        f.close()
    print("Done")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
