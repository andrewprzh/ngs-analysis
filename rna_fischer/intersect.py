import sys

def read_genes(inf):
    d = set()
    for l in open(inf):
        if l.strip() != "":
            d.add(l.strip().split()[0].split('.')[0])
    return d

genes1 = read_genes(sys.argv[1])
genes2 = read_genes(sys.argv[2])

for g in genes1:
    if g in genes2:
        print(g)
