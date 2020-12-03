############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
from Bio import SeqIO


def append_polya(input_fasta, output_fasta, poly_len):
    new_fasta = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        record.seq += "A" * poly_len
        new_fasta.append(record)

    SeqIO.write(new_fasta, output_fasta, "fasta")


append_polya(sys.argv[1], sys.argv[2], 80)
