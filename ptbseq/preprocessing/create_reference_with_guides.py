############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from traceback import print_exc


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--guide_tsv", "-g", help="TSV with guide sequences", type=str)
    parser.add_argument("--seq_column", help="column with guide sequence in TSV", type=int, default=2)
    parser.add_argument("--head_fasta", help="CROP-seq head sequence", type=str)
    parser.add_argument("--tail_fasta", help="CROP-seq tail sequence", type=str)
    parser.add_argument("--output", "-o", help="output FASTA", type=str)

    args = parser.parse_args()

    return args


def read_tsv_guide(inf):
    result = []
    for l in open(inf, "r"):
        if l.startswith("#"):
            continue
        result.append(l.strip().split())

    return result


def generate_sequence(tsv_entry, column, head_seq, tail_seq):
    seq = head_seq + tsv_entry[column] + tail_seq
    header = "_".join(tsv_entry[:column] + tsv_entry[column+1:])
    gene_id = "GENE_" + header
    t_id = "TRANSCRIPT_" + header
    start = len(head_seq) - 10
    end = len(head_seq) + len(tsv_entry[column]) + 10
    gtf_gene = "\t".join([header, "GuideRNA", "gene", str(start), str(end), ".", "+", ".", 'gene_id "' + gene_id + '";'])
    gtf_transcript = "\t".join([header, "GuideRNA", "transcript", str(start), str(end), ".", "+", ".",
                                'gene_id "' + gene_id + '"; transcript_id "' + t_id + '";'])
    gtf_exon = "\t".join([header, "GuideRNA", "exon", str(start), str(end), ".", "+", ".",
                          'gene_id "' + gene_id + '"; transcript_id "' + t_id + '";'])
    return SeqRecord.SeqRecord(seq=Seq.Seq(seq), id=header, description=""), "\n".join([gtf_gene, gtf_transcript, gtf_exon])


def main():
    args = parse_args()
    head_seq = str(list(SeqIO.parse(args.head_fasta, "fasta"))[0].seq).upper()
    tail_seq = str(list(SeqIO.parse(args.tail_fasta, "fasta"))[0].seq).upper()
    guide_sequences = read_tsv_guide(args.guide_tsv)
    output_records = []
    with open(args.output + ".gtf", "w") as output_gtf:
        for guide_seq in guide_sequences:
            entry = generate_sequence(guide_seq, args.seq_column, head_seq, tail_seq)
            output_records.append(entry[0])
            output_gtf.write(entry[1] + "\n")
    with open(args.output, "w") as output_handle:
        SeqIO.write(output_records, output_handle, "fasta")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
