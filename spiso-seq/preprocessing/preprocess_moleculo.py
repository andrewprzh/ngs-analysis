import os
import sys

from common import *

def get_moleculo_barcode(contig_id):
    delimeter = "-Barcode="
    tokens = contig_id.split(delimeter)
    if len(tokens) != 2:
        print("Wrong fromat " + contig_id)

    return tokens[0], tokens[1].strip().split(':')[0] if len(tokens) > 1 else ""


if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <FASTA file> <GMAP index> <output dir>")
    exit(0)

contigs_file = sys.argv[1]
gmap_index = sys.argv[2]
output_dir = sys.argv[3]

contigs_name, short_id_contigs_name = convert_fasta_with_barcodes(contigs_file, output_dir, get_moleculo_barcode)

align_fasta(output_dir, contigs_name, short_id_contigs_name, gmap_index, "star")

