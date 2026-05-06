import sys
import gzip

gb = set()
gbu =set()
spliced = 0
total_reads = 0

handle = gzip.open(sys.argv[1], "rt") if sys.argv[1].endswith('.gz') else open(sys.argv[1])

for l in handle:
    v = l.strip().split()
    total_reads += 1
    if v[5] != ";%;": spliced += 1
    gb.add((v[1], v[3]))
    gbu.add((v[1], v[3], v[4]))

print(total_reads, spliced, len(gb), len(gbu))
