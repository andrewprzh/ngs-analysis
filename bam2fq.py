import sys
import pysam
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

open(sys.argv[2], "w").close()
handle = open(sys.argv[2], "a")
buffer = []
processed_ids = set()

for a in pysam.AlignmentFile(sys.argv[1], "r", check_sq=False):
    seq = a.get_forward_sequence()
    if not seq:
        continue
    read_id = a.query_name
    if read_id in processed_ids:
        continue
    processed_ids.add(read_id)
    qual = a.get_forward_qualities()

    buffer.append(SeqRecord.SeqRecord(seq=Seq.Seq(seq), id=read_id, description="", letter_annotations={'phred_quality':qual}))

    if len(buffer) > 100000:
        SeqIO.write(buffer, handle, "fastq")
        buffer = []

SeqIO.write(buffer, handle, "fastq")
handle.close()
