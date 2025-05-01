import sys
import gffutils



t_len = {}
db = gffutils.FeatureDB(sys.argv[1])
for gene in db.features_of_type("gene"):
    for t in db.children(gene, featuretype=('transcript', 'mRNA')):
        exons = []
        for e in db.children(t):
            if e.featuretype == 'exon':
                exons.append((e.start, e.end))

        t_len[t.id] = sum(map(lambda x: x[1] - x[0] + 1, exons))

for t_id in sorted(t_len.keys()):
    print("%s\t%d" % (t_id, t_len[t_id]))
