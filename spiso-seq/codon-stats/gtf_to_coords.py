import os
import sys
import extract_exons


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <GTF> <output coord file>")
    exit(0)

gtf = open(sys.argv[1])
out_name = sys.argv[2]
out_f = open(out_name, 'w')

gene_info_list = []
current_gene = None
count = 0

for l in gtf:
    count += 1
    if count % 10000 == 0:
        sys.stderr.write("\r   " + str(count) + " lines processed")

    if l.startswith("#"):
        continue
    feature = Feature(l)

    if feature.feature_type == "gene":
        if current_gene is not None:
            gene_info_list.append(GeneInfo(current_gene))
        current_gene = Gene(feature)
    else:
        current_gene.add_feature(feature)

    if len(codon_info_list) > 1000:
        for c in codon_info_list:
            out_f.write(c.to_string())
        del codon_info_list[:]

codon_info_list.append(CodonInfo(current_gene))
stats.add_stats(codon_info_list)
for c in codon_info_list:
    out_f.write(c.to_string())
out_f.close()


