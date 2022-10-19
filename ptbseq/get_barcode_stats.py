import sys
from collections import defaultdict


gene_dict = defaultdict(list)
for l in open(sys.argv[1]):
    v = l.strip().split()
    gene_dict[v[1]].append((v[2], v[0]))


gene_counts = defaultdict(int)
gene_counts_filt = defaultdict(int)

for g in gene_dict.keys():
    gene_set = set([x[0] for x in gene_dict[g]]) # all genes for this batcode
    gene_counts[len(gene_set)] += 1
    if len(gene_set) > 1:
        gene2umi = defaultdict(list)
        for x in gene_dict[g]:
            gene2umi[x[0]].append(x[1])
        max_umis = 0
        supported_by_2 = 0
        for bc in gene2umi.keys():
            max_umis = max(max_umis, len(gene2umi[bc]))
            if len(gene2umi[bc]) > 1:
                supported_by_2 += 1
        if max_umis > 1:
            gene_counts_filt[supported_by_2] += 1
        else:
            gene_counts_filt[len(gene_set)] += 1
    else:
        gene_counts_filt[len(gene_set)] += 1

print(gene_counts)
print(gene_counts_filt)

print(sum(gene_counts.values()))
print(sum(gene_counts_filt.values()))



