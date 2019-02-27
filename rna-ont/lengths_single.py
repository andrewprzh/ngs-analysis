#!/usr/bin/env python


from Bio import SeqIO
import sys

def get_len_distr(in_fasta):
    f = SeqIO.index(in_fasta, "fasta")
    print(len(f['NODE_1_length_25851_cov_13.576154_NODE_g0_i0'].seq))

    for record in f:
        print(record)
        


get_len_distr(sys.argv[1])


