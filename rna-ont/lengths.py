#!/usr/bin/env python


from Bio import SeqIO
import sys

def get_len_distr(in_fasta):
    ld = {}
    for record in SeqIO.parse(in_fasta, "fasta"):
        length = len(record.seq)
        if length not in ld:
            ld[length] = 0
        ld[length] += 1
    return ld


ld1 = get_len_distr(sys.argv[1])
ld2 = get_len_distr(sys.argv[2])

all_lengths = set()
for length in ld1.keys():
    all_lengths.add(length)
for length in ld2.keys():
    all_lengths.add(length)

for length in range(1, max(all_lengths)):
    l1 = 0 if length not in ld1 else ld1[length]
    l2 = 0 if length not in ld2 else ld2[length]
    print(str(length) + "\t" + str(l1) + "\t" + str(l2))
