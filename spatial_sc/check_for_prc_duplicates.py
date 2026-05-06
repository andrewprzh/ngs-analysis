import sys
from collections import defaultdict

gene_barcode_pairs = defaultdict(list)
umis = defaultdict(int)
c = 0

for l in open(sys.argv[1]):
#m84039_250130_084836_s3/168629832/ccs/8052_8677 ENSG00000290825.2       None    s_002um_02478_01877     GGGGGTACA       ;%;     NoTSS   NoPolyA ;%;chr1_14400_14696_-   known_ambiguous 0       ENST00000832823.1       lncRNA
    v = l.strip().split()
    gene_barcode_pairs[(v[1],v[3])].append(v[4])
    umis[(v[1],v[3],v[4], v[8][3:].split('_')[0])] += 1
    c += 1

#    if c > 100000: break

print(len(umis))
for t in umis.keys():
    if umis[t] > 1:
        print(t)
#for k in gene_barcode_pairs.keys():
#    if len(gene_barcode_pairs[k]) > 1:
#        print(k, gene_barcode_pairs[k])
