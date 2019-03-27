import os
import sys
from plot_common import *


if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <gene id> <gene coord file> <contig coord file> [output figure name = gene_id]")
    exit(0)

gene_id = sys.argv[1]
gene_coord_file = open(sys.argv[2])
contig_coord_file = open(sys.argv[3])
out_name = sys.argv[4] if len(sys.argv) > 4 else gene_id

found_gene = False
lg = gene_coord_file.readline()
while lg and not found_gene:
    if lg.startswith("gene\t" + gene_id):
        found_gene = True
        break
    lg = gene_coord_file.readline()

found_contigs = False
lc = contig_coord_file.readline()
while lc and not found_contigs:
    if lc.startswith("gene\t" + gene_id):
        found_contigs = True
        break
    lc = contig_coord_file.readline()


if not found_gene and not found_contigs:
    print("Gene was not detected in both files")
    exit(0)

total_contigs = 0
total_transcripts = 0
gene_coords = (0, 0)
chr_id = ""
if found_contigs:
    tokens = lc.strip().split('\t')
    if len(tokens) != 7:
        print("Wrong format")
        exit(0)
    gene_coords = (int(tokens[4]), int(tokens[5]))
    total_contigs = int(tokens[6])
    chr_id = tokens[2]

if found_gene:
    tokens = lg.strip().split('\t')
    if len(tokens) != 8:
        print("Wrong format")
        exit(0)
    gene_coords = (int(tokens[4]), int(tokens[5]))
    total_transcripts = int(tokens[6])
    chr_id = tokens[2]


drawer = GeneDrawer(total_transcripts, total_contigs, gene_coords, out_name)
drawer.draw_gene(gene_coords)
drawer.label_genome(gene_coords, gene_id, chr_id)

lg = gene_coord_file.readline()
transcript_num = 0
while lg and not lg.startswith("gene"):
    tokens = lg.strip().split('\t')
    if len(tokens) < 3:
        print("Wrong format")
        break
    coords = (int(tokens[1]), int(tokens[2]))
    if tokens[0] == "transcript":
        transcript_num += 1
        drawer.draw_transcript(coords, transcript_num)
        drawer.label_alignment(coords, transcript_num, tokens[-1])
    elif tokens[0] == "CDS":
        drawer.draw_exon(coords, transcript_num)        
    elif tokens[0] == "start_codon":
        drawer.draw_start_codon(coords, transcript_num, -1)
        drawer.draw_start_codon(coords)                
    elif tokens[0] == "stop_codon":
        drawer.draw_stop_codon(coords, transcript_num, -1)        
        drawer.draw_stop_codon(coords)   
    lg = gene_coord_file.readline()

lc = contig_coord_file.readline()
contig_num = 0
cov = 1
while lc and not lc.startswith("gene"):
    tokens = lc.strip().split('\t')
    if len(tokens) < 3:
        print("Wrong format")
        break
    coords = (int(tokens[1]), int(tokens[2]))
    if tokens[0] == "contig":
        contig_num += 1
        cov = int(tokens[3])
        drawer.draw_contig(coords, contig_num)
        drawer.label_alignment(coords, contig_num, tokens[-1])
    elif tokens[0] == "block":
        drawer.draw_block(coords, contig_num, cov)        
    elif tokens[0] == "start_codon":
        if contig_num == 0:
            drawer.draw_start_codon(coords)                
        else:
            drawer.draw_start_codon(coords, contig_num, 1)
    elif tokens[0] == "stop_codon":
        if contig_num == 0:
            drawer.draw_stop_codon(coords)   
        else:
            drawer.draw_stop_codon(coords, contig_num, 1)        
    lc = contig_coord_file.readline()
    

