############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import math
import os
import sys
import argparse
from Bio import SeqIO
import random
from common import *
from traceback import print_exc
import numpy


POLYA_LEN = 30
LINKER = "TCTTCAGCGTTCCCGAGA"
UPS_PRIMER = "AAGCAGTGGTATCAACGCAGAGT"
PCR_PRIMER = "CTACACGACGCTCTTCCGATCT"
UPS_PRIMER_REV = reverese_complement(UPS_PRIMER)
LEFT_BC_LENGTH = 8
RIGHT_BC_LENGTH = 6
BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
UMI_LENGTH = 9

TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
R1 = "CTACACGACGCTCTTCCGATCT" # R1
BARCODE_LEN_10X = 16
UMI_LEN_10X = 12
UMI_LEN_10X_V2 = 10


STEREO_LINKER = "TTGTCTTCCTAAGAC"
TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
STEREO_BC_LEN = 25
STEREO_UMI_LEN = 10

# defaults for --strand_switch_prob and --concat_probabilities (P of 1, 2, 3, 4 molecules)
CONCAT_STRAND_SWITCH_PROB = 0.4
CONCAT_PROBABILITIES = [0.5, 0.35, 0.1, 0.05]

NUCLS = ['A', 'C', 'G', 'T']
VNUCLS =  ['A', 'C', 'G']


def get_random_seq(length, non_t_tail=0):
    seq = ""
    for i in range(length-non_t_tail):
        index = random.randint(0, len(NUCLS) - 1)
        seq += NUCLS[index]
    for i in range(non_t_tail):
        index = random.randint(0, len(VNUCLS) - 1)
        seq += VNUCLS[index]
    return seq


def create_template_spatial(sequence, barcode, umi):
    assert len(barcode) == BC_LENGTH
    barcoded_part = PCR_PRIMER + barcode[:LEFT_BC_LENGTH] + LINKER + barcode[LEFT_BC_LENGTH:] + umi
    template_seq = barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + UPS_PRIMER_REV
    return reverese_complement(template_seq)


def create_template_stereo(sequence, barcode, umi):
    barcoded_part = PC1_PRIMER + barcode + STEREO_LINKER + umi
    template_seq = barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + TSO5
    return reverese_complement(template_seq)


VISIUM_HD_3PRIMER = "CTACACGACGCTCTTCCGATCT"
VISIUM_HD_5PRIMER = "ATGTACTCTGCGTTGATACCACTGCTT"
VISIUM_HD_BC_LEN = 30
VISIUM_HD_UMI_LEN = 9


def create_template_visiumhd(sequence: str, barcode: str, umi:str):
    barcode = barcode.replace("|", "")
    barcoded_part = VISIUM_HD_3PRIMER + umi + barcode
    template_seq = barcoded_part + get_random_seq(2, 2) + "T" * POLYA_LEN + reverese_complement(sequence) + VISIUM_HD_5PRIMER
    return reverese_complement(template_seq)


#AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT—--DNA------CCCATGTACTCTGCGTTGATACCACTGCTT

# p5="AATGATACGGCGACCACCGAGATCTACAC"
# i5="NNNNNNNNNN"
# truseqread1="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# 10xbc="NNNNNNNNNNNNNNNN"
# umi="NNNNNNNNNNNN"
# polyt="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
# tso="CCCATGTACTCTGCGTTGATACCACTGCTT"
def create_template_10x(sequence, barcode, umi):
    assert len(barcode) == BARCODE_LEN_10X
    barcoded_part = R1 + barcode + umi
    template_seq =  barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + TSO
    return reverese_complement(template_seq)


def load_counts(inf, total_count):
    count_dict = {}
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split("\t")
        count_dict[v[0]] = math.ceil(float(v[2]) * total_count / 1000000)
    return count_dict


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--transcriptome", "-t", help="reference transcriptome", required=True, type=str)
    parser.add_argument("--counts", "-c", help="expression table", type=str, required=True)
    parser.add_argument("--template_count", help="total number of templates to generate", type=int, default=100000)
    parser.add_argument("--barcodes", "-b", help="barcode list (random if not set)", type=str)
    parser.add_argument("--umis", "-u", help="UMI list (random if not set)", type=str)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--mode", help="[spatial | 10x | 10x_v2 | stereo | visium_hd]", type=str, default="stereo")
    parser.add_argument("--concatenate_templates", action="store_true", help="concatenate different cDNA templates together", default=False)
    parser.add_argument("--concat_probabilities", type=str, default=",".join(map(str, CONCAT_PROBABILITIES)),
                        help="comma-separated probabilities of 1, 2, 3, ... molecules per template "
                             "(with --concatenate_templates), default %(default)s")
    parser.add_argument("--strand_switch_prob", type=float, default=CONCAT_STRAND_SWITCH_PROB,
                        help="probability of each concatenated molecule to be reverse-complemented "
                             "(with --concatenate_templates), default %(default)s")
    parser.add_argument("--seed", type=int, help="random seed")

    args = parser.parse_args()
    args.concat_probabilities = [float(p) for p in args.concat_probabilities.split(",")]
    if any(p < 0 for p in args.concat_probabilities) or abs(sum(args.concat_probabilities) - 1.0) > 1e-6:
        print("--concat_probabilities must be non-negative and sum to 1, got %s" % args.concat_probabilities)
        exit(-1)
    if not 0 <= args.strand_switch_prob <= 1:
        print("--strand_switch_prob must be within [0, 1], got %f" % args.strand_switch_prob)
        exit(-1)
    return args


def main():
    args = parse_args()
    if args.seed is not None:
        random.seed(args.seed)
        numpy.random.seed(args.seed)
    isoforms = SeqIO.to_dict(SeqIO.parse(args.transcriptome, "fasta"))
    print("Loaded %d reference sequences from %s" % (len(isoforms), args.transcriptome))
    barcodes = load_barcodes(args.barcodes) if args.barcodes else []
    print("Using %d reference barcodes from %s" % (len(barcodes), args.barcodes))
    umis = load_barcodes(args.umis) if args.umis else []
    print("Using %d reference UMIs from %s" % (len(umis), args.umis))
    count_dict = load_counts(args.counts, args.template_count)
    print("Loaded %d counts from %s" % (len(count_dict), args.counts))

    if args.mode == "10x":
        bc_len = BARCODE_LEN_10X
        umi_len = UMI_LEN_10X
        template_func = create_template_10x
    elif args.mode == "10x_v2":
        bc_len = BARCODE_LEN_10X
        umi_len = UMI_LEN_10X_V2
        template_func = create_template_10x
    elif args.mode == "spatial":
        bc_len = BC_LENGTH
        umi_len = UMI_LENGTH
        template_func = create_template_spatial
    elif args.mode == "stereo":
        bc_len = STEREO_BC_LEN
        umi_len = STEREO_UMI_LEN
        template_func = create_template_stereo
    elif args.mode == "visium_hd":
        bc_len = VISIUM_HD_BC_LEN
        umi_len = VISIUM_HD_UMI_LEN
        template_func = create_template_visiumhd
    else:
        print("Unknown mode %s" % args.mode)
        exit(-1)

    scale_factor = 1000000.0 / sum(count_dict.values())
    counts_file = args.output + ".counts.tsv"
    with open(counts_file, "w") as expr_table:
        expr_table.write("#transcript_id\tcounts\ttpm\n")
        for t_id in count_dict:
            expr_table.write("%s\t%d\t%.4f\n" % (t_id, count_dict[t_id], scale_factor * count_dict[t_id]))
    print("Saved templates counts to %s" % counts_file)

    template_fname = args.output + ".templates.fasta"
    output_fasta = open(template_fname, "w")
    output_counts = open(args.output + ".template_counts.tsv", "w")
    # NanoSim treats the first line of the expression file as a header and skips it
    output_counts.write("#template_id\tcounts\ttpm\n")
    concat_counts = list(range(1, len(args.concat_probabilities) + 1))
    for t_id in count_dict:
        transcript_seq = isoforms[t_id]
        template_list = []
        for i in range(count_dict[t_id]):
            seq_id = "READ_%d_%s" % (i, t_id)
            if args.concatenate_templates:
                template_count = numpy.random.choice(concat_counts, p=args.concat_probabilities)
            else:
                template_count = 1

            template = ""
            for j in range(template_count):
                bc = barcodes[random.randint(0, len(barcodes) - 1)] if barcodes else get_random_seq(bc_len)
                umi = umis[random.randint(0, len(umis) - 1)] if umis else get_random_seq(umi_len, 0)
                sub_template = template_func(transcript_seq, bc, umi)
                strand = "+"
                if args.concatenate_templates and numpy.random.random() < args.strand_switch_prob:
                    strand = "-"
                    sub_template = reverese_complement(sub_template)

                seq_id += "_%s_%s_%s" % (bc, umi, strand)
                template += sub_template
            output_counts.write("%s\t%d\t%.4f\n" % (seq_id, 1, scale_factor))
            template_list.append(SeqIO.SeqRecord(seq=Seq.Seq(template), id=seq_id, description="", name=""))

        SeqIO.write(template_list, output_fasta, "fasta")
    output_fasta.close()
    output_counts.close()
    print("Saved templates counts to %s" % template_fname)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
