import sys

# usage: <read-barcode-umi> <read_assignment.tsv>

read_bc_umi = {}
for l in open(sys.argv[1]):
    v = l.strip().split('\t')
    read_bc_umi[v[0][1:]] = (v[1], v[2])

bc_umi_gene = set()
for l in open(sys.argv[2]):
    v = l.strip().split('\t')
    if v[3] == '.':
        continue
    read_id = v[0]
    if read_id in read_bc_umi:
        print("%s\t%s\t%s\t%s" % (read_id, v[3], read_bc_umi[read_id][0], read_bc_umi[read_id][1]))
        bc_umi_gene.add((v[3], read_bc_umi[read_id][0], read_bc_umi[read_id][1]))

#for t in bc_umi_gene:
#    print("%s\t%s\t%s" % t)