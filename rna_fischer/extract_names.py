import sys
import os
import gffutils


gene_set = set()
for l in open(sys.argv[1]):
    gene_set.add(l.strip().upper())

genes_found = set()
name_set = set()
db = gffutils.FeatureDB(sys.argv[2], keep_order = True)
for g in db.features_of_type('gene', order_by=('seqid', 'start')):
    gene_db = db[g.id]
    gene_name = gene_db.attributes['gene_name'][0].upper()
    gene_id = g.id.upper()
    if gene_id in gene_set:
        genes_found.add(gene_id)
        name_set.add(gene_name)

for gene_name in name_set:
    print(gene_name)

for gene_id in gene_set:
    if gene_id not in genes_found:
        print("Unable to find gene: " + gene_id)





