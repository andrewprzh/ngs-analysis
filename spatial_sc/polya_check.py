import sys
from collections import defaultdict

gene_dict = defaultdict(lambda: defaultdict(int))
for l in open(sys.argv[1]):
    v = l.strip().split()
    c = int(v[0])
    g = v[1]
    p = v[2]
    if v == 'NoPolyA': continue

    gene_dict[g][p] = c


for l in open(sys.argv[2]):
    v = l.strip().split('\t')
    g = v[0]
    p1 = v[6]
    p2 = v[7]

    if g not in gene_dict: continue
    other_counts = []
    used_counts = []

    for p in gene_dict[g].keys():
        if p not in [p1, p2]:
            other_counts.append(gene_dict[g][p])
        else:
            used_counts.append(gene_dict[g][p])

    print(p1, p2)
    print(used_counts)
    print(other_counts)
    print(min(used_counts) / max(other_counts))

