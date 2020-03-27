############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


inf = open(sys.argv[1])
l = "1"

while l:
    reads = inf.readline().strip().replace('/', '-')
    pos = inf.readline().strip()

    seq1 = inf.readline().strip()
    seq2 = inf.readline().strip()
    seq1_file = reads.split()[0] + ".fasta"
    seq2_file = reads.split()[1] + ".fasta"

    record = SeqRecord(Seq(seq1, IUPACAmbiguousDNA()), name=reads.split()[0])
    SeqIO.write([record], seq1_file, "fasta")
    record = SeqRecord(Seq(seq2, IUPACAmbiguousDNA()), name=reads.split()[1])
    SeqIO.write([record], seq2_file, "fasta")

    os.system("/Bmo/prjbel/tools/minimap2-2.8_x64-linux/minimap2 "
              "-t 1 -x ava-ont " + seq2_file + " " + seq1_file + " > " + seq1_file + "_" + seq2_file + ".sam")
    break

