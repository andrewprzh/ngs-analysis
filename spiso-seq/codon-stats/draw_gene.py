import os
import sys
import plot_common


if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <gene coord file> <gene id> <output figure name>")
    exit(0)

inf = open(sys.argv[1])
gene_id = sys.argv[2]
out_name = sys.argv[3]

l = f.readline()
while l:
    if not l.startswith("gene\t" + gene_id):
        l = f.readline()
        continue

    tokens = l.strip().split('\t')
    if len(tokens) != 8:
        print("Wrong format")
        break

    init(int(tokens[6]), (int(tokens[4]), int(tokens[5])))
    context = init_context(out_name)

    draw_gene()
    while l and not l.startswith("gene"):
        l = f.readline()
        
    

