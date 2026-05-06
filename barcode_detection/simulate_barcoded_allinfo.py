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


STEREO_LINKER = "TTGTCTTCCTAAGAC"
TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
STEREO_BC_LEN = 25
STEREO_UMI_LEN = 10

CONCAT_STRAND_SWITCH_PROB = 0.4
CONCAT_PROBABILITIES = [0.5, 0.35, 0.1, 0.05]
CONCAT_COUNT = [1,2,3,4]

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
VISIUM_HD_UMI_LEN = 10


def create_template_visiumhd(sequence: str, barcode: str, umi:str):
    barcode = barcode.replace("|", "")
    barcoded_part = VISIUM_HD_3PRIMER + umi + barcode
    template_seq = barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + VISIUM_HD_5PRIMER
    return reverese_complement(template_seq)


def create_template_10x(sequence, barcode, umi):
    assert len(barcode) == BARCODE_LEN_10X
    barcoded_part = R1 + barcode + umi
    template_seq =  barcoded_part + "T" * POLYA_LEN + reverese_complement(sequence) + TSO
    return reverese_complement(template_seq)


def load_cell_type_assignments(inf):
    """
    Load cell type assignment TSV file with columns:
    gene_id, cell_type, barcode, umi, isoform_id
    """
    assignments = []
    with open(inf) as f:
        #header = f.readline().strip().split('\t')
        #print("Header columns:", header)
        
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            if parts[4] in {'novel', 'None'}: continue
            if parts[3] == 'None':
                parts[3] = get_random_seq(10, 2)

            if parts[1][-1].isdigit():
                parts[1] = parts[1][:-1]
                
            assignment = {
                'gene_id': parts[0],
                'cell_type': parts[1],
                'barcode': parts[2],
                'umi': parts[3],
                'isoform_id': parts[4]
            }
            assignments.append(assignment)
    
    return assignments


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--transcriptome", "-t", help="reference transcriptome", required=True, type=str)
    parser.add_argument("--cell_type_assignments", "-a", help="cell type assignment TSV file", type=str, required=True)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--mode", help="[spatial | 10x | stereo | visium_hd]", type=str, default="10x")
    parser.add_argument("--concatenate_templates", action="store_true", help="concatenate different cDNA templates together", default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    # Load transcriptome
    isoforms = SeqIO.to_dict(SeqIO.parse(args.transcriptome, "fasta"))
    #isoforms = {k.split('.')[0] : raw_isoforms[k] for k in raw_isoforms.keys()}
    print(list(isoforms.keys())[:10])
    print("Loaded %d reference sequences from %s" % (len(isoforms), args.transcriptome))
    
    # Load cell type assignments
    assignments = load_cell_type_assignments(args.cell_type_assignments)
    print("Loaded %d cell type assignments from %s" % (len(assignments), args.cell_type_assignments))

    # Set mode-specific parameters
    if args.mode == "10x":
        bc_len = BARCODE_LEN_10X
        umi_len = UMI_LEN_10X
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

    # Open output files
    template_fname = args.output + ".templates.fasta"
    counts_fname = args.output + ".counts.tsv"
    template_counts_fname = args.output + ".template_counts.tsv"
    
    output_fasta = open(template_fname, "w")
    output_counts = open(template_counts_fname, "w")
    
    # Count occurrences for expression table
    isoform_counts = {}
    missed_isoforms = set()
    
    # Process each assignment
    for idx, assignment in enumerate(assignments):
        gene_id = assignment['gene_id']
        cell_type = assignment['cell_type']
        barcode = assignment['barcode']
        umi = assignment['umi']
        isoform_id = assignment['isoform_id']
        
        # Check if isoform exists in transcriptome
        if isoform_id not in isoforms:
            missed_isoforms.add(isoform_id)
            #print("Warning: isoform %s not found in transcriptome, skipping" % isoform_id)
            continue
        
        transcript_seq = isoforms[isoform_id]
        
        # Create sequence ID with gene_id and cell_type
        seq_id = "READ_%d_%s_%s_%s" % (idx, gene_id, cell_type, isoform_id)
        
        # Handle template concatenation if enabled
        if args.concatenate_templates:
            template_count = numpy.random.choice(CONCAT_COUNT, p=CONCAT_PROBABILITIES)
        else:
            template_count = 1
        
        template = ""
        for j in range(template_count):
            sub_template = template_func(transcript_seq, barcode, umi)
            strand = "+"
            if args.concatenate_templates and numpy.random.random() < CONCAT_STRAND_SWITCH_PROB:
                strand = "-"
                sub_template = reverese_complement(sub_template)
            
            seq_id += "_%s_%s_%s" % (barcode, umi, strand)
            template += sub_template
        
        # Write template to FASTA
        template_record = SeqIO.SeqRecord(seq=Seq.Seq(template), id=seq_id, description="", name="")
        SeqIO.write([template_record], output_fasta, "fasta")
        
        # Write to template counts (count=1, TPM=1.0 for each template)
        output_counts.write("%s\t%d\t%.4f\n" % (seq_id, 1, 1.0))
        
        # Count for expression table
        isoform_ct = (isoform_id, cell_type)
        if isoform_ct not in isoform_counts:
            isoform_counts[isoform_ct] = 0
        isoform_counts[isoform_ct] += 1
    
    output_fasta.close()
    output_counts.close()
    print("Isoforms missed: %d" % len(missed_isoforms))
    
    print("Saved templates to %s" % template_fname)
    print("Saved template counts to %s" % template_counts_fname)
    
    # Write expression table with actual counts and TPMs
    total_count = sum(isoform_counts.values())
    scale_factor = 1000000.0 / total_count if total_count > 0 else 0
    
    with open(counts_fname, "w") as expr_table:
        expr_table.write("#transcript_id\tcounts\ttpm\n")
        for isoform_ct in sorted(isoform_counts.keys()):
            count = isoform_counts[isoform_ct]
            tpm = scale_factor * count
            expr_table.write("%s\t%s\t%d\t%.4f\n" % (isoform_ct[0], isoform_ct[1], count, tpm))
    
    print("Saved expression counts to %s" % counts_fname)
    print("Processed %d templates total" % len(assignments))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
