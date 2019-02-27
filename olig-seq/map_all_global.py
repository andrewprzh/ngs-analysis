import os
import sys
from pyfasta import Fasta
import ssw
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <oligs in FASTA> <processed reads in FASTA> > <output.txt>")
    exit(0)

oligs = Fasta(sys.argv[1])
reads = Fasta(sys.argv[2])

aligner = ssw.Aligner()
for r in reads.keys():
    if len(str(reads[r][:])) != 169:
        continue
    max_score = 0
    best_al = None
    best_seqs = None
    for o in oligs.keys():
        score = pairwise2.align.globalxx(str(oligs[o][:]), str(reads[r][:]), score_only = True)
        if max_score < score:
            max_score = score
            best_seqs = (r, o)

    best_al = pairwise2.align.globalxx(str(oligs[best_seqs[1]][:]), str(reads[best_seqs[0]][:]), score_only = False)[0]
    print(best_seqs[0] + "\t" + best_seqs[1])
    print(format_alignment(*best_al))


