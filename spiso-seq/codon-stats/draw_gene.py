import os
import sys
from plot_common import *


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <gene coord file> <gene id> [output figure name = gene_id]")
    exit(0)

inf = open(sys.argv[1])
gene_id = sys.argv[2]
out_name = sys.argv[3] if len(sys.argv) > 3 else gene_id

found_gene = False
l = inf.readline()
while l and not found_gene:
    if not l.startswith("gene\t" + gene_id):
        l = inf.readline()
        continue

    tokens = l.strip().split('\t')
    if len(tokens) != 8:
        print("Wrong format")
        break

    found_gene = True
    gene_coords = (int(tokens[4]), int(tokens[5]))
    total_transcripts = int(tokens[6])
    drawer = GeneDrawer(total_transcripts, gene_coords, out_name)
    drawer.draw_gene(gene_coords)

    l = inf.readline()
    transcript_num = 0
    while l and not l.startswith("gene"):
        tokens = l.strip().split('\t')
        if len(tokens) != 3:
            print("Wrong format")
            break
        coords = (int(tokens[1]), int(tokens[2]))
        if tokens[0] == "transcript":
            transcript_num += 1
            drawer.draw_transcript(coords, transcript_num)
        elif tokens[0] == "CDS":
            drawer.draw_exon(coords, transcript_num)        
        elif tokens[0] == "start_codon":
            drawer.draw_start_codon(coords, transcript_num)
            drawer.draw_start_codon(coords)                
        elif tokens[0] == "stop_codon":
            drawer.draw_stop_codon(coords, transcript_num)        
            drawer.draw_stop_codon(coords)   
        l = inf.readline()
    

