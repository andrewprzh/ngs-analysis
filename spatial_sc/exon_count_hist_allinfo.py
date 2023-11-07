import sys
from collections import defaultdict

exon_counts = defaultdict(int)
total_exons = 0
total_exons_in_spliced = 0

for l in open(sys.argv[1]):
    v = l.strip().split("\t")
    exon_count = v[8].count('%')
    total_exons += exon_count
    exon_counts[exon_count] += 1
    if exon_count > 1:
        total_exons_in_spliced += exon_count

read_count = sum(exon_counts.values())
spliced_reads = sum([v if k > 1 else 0 for k, v in exon_counts.items()])

print("Total reads: %d, total exons: %d, mean:%.2f" % (read_count, total_exons, total_exons / read_count))
print("Total spliced reads: %d, total exons: %d, mean:%.2f" % (spliced_reads, total_exons_in_spliced, total_exons_in_spliced / spliced_reads))

print("Exon count histogram:")
for k in sorted(exon_counts.keys()):
    print("%d\t%d" % (k, exon_counts[k]))

