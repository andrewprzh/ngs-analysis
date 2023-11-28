import sys
from collections import defaultdict
import numpy

exon_counts = defaultdict(int)
total_exons = 0
total_exons_in_spliced = 0
gene_count_dict = defaultdict(int)
spliced_exon_counts = []


for l in open(sys.argv[1]):
    v = l.strip().split("\t")
    exon_count = v[8].count('%')
    total_exons += exon_count
    exon_counts[exon_count] += 1
    if exon_count > 1:
        spliced_exon_counts.append(exon_count)
        total_exons_in_spliced += exon_count
        gene_id = v[1]
        gene_count_dict[gene_id] += 1

read_count = sum(exon_counts.values())
spliced_reads = sum([v if k > 1 else 0 for k, v in exon_counts.items()])

print("Total reads: %d, total exons: %d, mean:%.2f" % (read_count, total_exons, total_exons / read_count))
print("Total spliced reads: %d, total exons: %d, mean:%.2f" % (spliced_reads, total_exons_in_spliced, total_exons_in_spliced / spliced_reads))

print("Exon count histogram:")
for k in sorted(exon_counts.keys()):
    print("%d\t%d" % (k, exon_counts[k]))

gene_counts = list(gene_count_dict.values())
for cutoff in [20, 50, 100]:
    print("Genes with >= %d spliced reads: %d" %(cutoff, len(list(filter(lambda x: x>=cutoff, gene_counts)))))

bins = [i * 10 for i in range(50)] + [100000]
h, b = numpy.histogram(gene_counts, bins)
print(h)


