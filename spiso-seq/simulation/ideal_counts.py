import sys
from Bio import SeqIO

transcript_list = []
for r in SeqIO.parse(sys.argv[1], "fasta"):
    transcript_list.append(r.id)

count = int(sys.argv[2]) if len(sys.argv) > 2 else 20
tpm = 1000000.0 / float(len(transcript_list))

print("target_id\test_counts\ttpm")
for t in transcript_list:
    print("%s\t%d\t%.8f\t" % (t, count, tpm))


