import numpy
import sys
from Bio import SeqIO

transcript_dict = {}
for r in SeqIO.parse(sys.argv[1], "fasta"):
    transcript_dict[r.id] = int(numpy.random.gamma(0.7, 200, 1))

scale_factor = 1000000.0 / sum(transcript_dict.values())

print("#target_id\test_counts\ttpm")
for t in transcript_dict.keys():
    print("%s\t%d\t%.8f\t" % (t, transcript_dict[t], scale_factor * transcript_dict[t]))


