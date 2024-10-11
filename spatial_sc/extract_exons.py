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
from traceback import print_exc
import glob
import gzip
from Bio import SeqIO, Seq, SeqRecord


REGION_LEN = 300


base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverse_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


def get_regions(fasta_record, start, end, strand):
    exon_region = fasta_record[start:end+1]
    if strand == "+":
        return exon_region
    else:
        return reverse_complement(exon_region)


def process_exon_file(exon_file, output_dir, fasta_dict):
    name = os.path.basename(exon_file)
    seq_buffer = []

    for l in open(exon_file):
        v = l.strip().split("\t")
        exon_id = v[0]
        exon_info = exon_id.split("_")
        assert len(exon_info) == 4
        chr_id, start, end, strand = exon_info[0], int(exon_info[1]), int(exon_info[2]), exon_info[3]

        exon_seq = get_regions(fasta_dict[chr_id].seq, start, end, strand)
        exon_str = "_".join(v)
        seq_buffer.append(SeqRecord.SeqRecord(seq=Seq.Seq(exon_seq), id="%s_exon" % exon_str, description=""))

    out_fasta = os.path.join(output_dir, name + ".exons.fasta")
    SeqIO.write(seq_buffer, out_fasta, 'fasta')


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", default="./")
    parser.add_argument("--reference", "-r", type=str, help="reference genome", required=True)
    parser.add_argument("--input", "-i", type=str, nargs='+',
                        help="one or more files with input exon lists (allows wildcards)", required=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print("Loading reference genome from %s" % args.reference)
    _, outer_ext = os.path.splitext(args.reference)
    if outer_ext.lower() in ['.gz', '.gzip']:
        with gzip.open(args.reference, "rt") as handle:
            reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        reference_record_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    for inf in args.input:
        if os.path.isdir(inf):
            for f in glob.glob(inf):
                if os.path.isfile(f):
                    print("Processing %s" % f)
                    process_exon_file(f, args.output, reference_record_dict)
        elif os.path.isfile(inf):
            print("Processing %s" % inf)
            process_exon_file(inf, args.output, reference_record_dict)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
