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


inf = open(sys.argv[1])
l = inf.readline()

while l:
    reads = l.strip().replace('/', '-')
    pos = inf.readline().strip()

    seq1 = inf.readline().strip()
    seq2 = inf.readline().strip()
    read_id1 = reads.split()[0]
    read_id2 = reads.split()[1]
    seq1_file = read_id1 + ".fasta"
    seq2_file = read_id2 + ".fasta"

    record = SeqRecord(Seq(seq1), id=read_id1, description=pos.split()[0])
    SeqIO.write([record], seq1_file, "fasta")
    record = SeqRecord(Seq(seq2), id=read_id2, description=pos.split()[1])
    SeqIO.write([record], seq2_file, "fasta")

    sam_out = "sam/" + read_id1 +  "_" + read_id2 + ".sam"
    os.system("/Bmo/prjbel/tools/minimap2-2.8_x64-linux/minimap2 "
              "-t 1 -x ava-ont " + seq2_file + " " + seq1_file + " > " + sam_out)

    if os.stat(sam_out).st_size == 0:
        os.remove(sam_out)
    l = inf.readline()

