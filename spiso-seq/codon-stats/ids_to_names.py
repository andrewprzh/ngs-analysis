import os
import sys

if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <GTF> <gene ID list>")
    exit(0)

gtf = open(sys.argv[1])
inf = open(sys.argv[2])

gene_dict = {}

for l in gtf:
    tokens = l.strip().split()
    if len(tokens) < 3:
        continue
    if tokens[2] == "gene":
        gene_id = ""
        gene_name = ""
        for i in range(3, len(tokens)):
            if tokens[i] == "gene_id":
                gene_id = tokens[i+1][1:-2]
            if tokens[i] == "gene_name":
                gene_name = tokens[i+1][1:-2]

        if gene_id != "" and gene_name != "":
            gene_dict[gene_id] = gene_name

for l in inf:
    print(gene_dict[l.strip()])

