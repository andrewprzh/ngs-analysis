import os
import sys


if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " <oligs in CSV> > <out FASTA>")
    exit(0)

header = True
seqid = 1
for l in open(sys.argv[1]):
    if header:
        header = False
        continue

    print('>Olig_' + str(seqid))
    print(l.split(',')[0])
    seqid += 1




