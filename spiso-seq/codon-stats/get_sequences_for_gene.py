import os
import sys
from Bio import SeqIO


def find_features(gene_id, feature_type, infile_name):
    infile = open(infile_name)
    found_gene = False
    l = infile.readline()
    while l and not found_gene:
        if l.startswith("gene\t" + gene_id):
            found_gene = True
            break
        l = infile.readline()

    if not found_gene:
        print("Gene was not detected in one in " + infile_name)
        return []

    feature_list = []

    l = infile.readline()
    while l and not l.startswith("gene"):
        tokens = l.strip().split('\t')
        if len(tokens) < 3:
            print("Wrong format")
            break
        if tokens[0] == feature_type:
            if len(tokens) < 4:
                print("Wrong format")
                break
            feature_list.append(tokens[-1])
        l = infile.readline()

    return feature_list


def extract_sequences_to_file(in_fasta_file, id_list, out_fasta_file):
    in_fasta_records = SeqIO.index(in_fasta_file, "fasta")
    selected_records = []
    for seq_id in id_list:
        selected_records.append(in_fasta_records[seq_id])
    SeqIO.write(selected_records, out_fasta_file, "fasta")


def process_gene(gene_id, gene_coord_file, contig_coord_file, cdna_fasta_file, contig_fasta_file, out_contigs, out_cdna):
    extract_sequences_to_file(cdna_fasta_file, find_features(gene_id, "transcript", gene_coord_file), out_cdna)
    extract_sequences_to_file(contig_fasta_file, find_features(gene_id, "contig", contig_coord_file), out_contigs)


def main():
    if len(sys.argv) < 6:
        print("Usage: " + sys.argv[0] + " <gene id> <gene coord file> <contig coord file> <cdna.fa> <contigs.fa>  [output folder = ./]")
        exit(0)

    gene_id = sys.argv[1]
    gene_coord_file = sys.argv[2]
    contig_coord_file = sys.argv[3]
    cdna_fasta_file = sys.argv[4]
    contig_fasta_file = sys.argv[5]
    out_dir = sys.argv[6] if len(sys.argv) > 6 else "./"
    out_contigs = os.path.join(out_dir, gene_id + "_contigs.fa")
    out_cdna = os.path.join(out_dir, gene_id + "_cdna.fa")
    process_gene(gene_id, gene_coord_file, contig_coord_file, cdna_fasta_file, contig_fasta_file, out_contigs, out_cdna)

if __name__ == "__main__":
   main()
