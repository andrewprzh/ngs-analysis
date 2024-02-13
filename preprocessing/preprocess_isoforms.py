import os
import sys

from common import *

def get_transcript_name(contig_id):
    tokens = contig_id.split(" ")
    return tokens[0], tokens[0]

if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <FASTA file> <GMAP index> <output dir>")
    exit(0)

contigs_file = sys.argv[1]
gmap_index = sys.argv[2]
output_dir = sys.argv[3]

contigs_name, short_id_contigs_name = convert_fasta_with_barcodes(contigs_file, output_dir, get_transcript_name)

align_fasta(output_dir, contigs_name, short_id_contigs_name, gmap_index, "star")

