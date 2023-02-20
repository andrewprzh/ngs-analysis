import os
import sys
import pysam
import numpy


def read_reads(inf):
    read_set = set()
    for l in open(inf):
        read_set.add(l.strip().split('\t')[0])
    return read_set


def filter_reads(in_file_name, out_file_name, read_set):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    outf = pysam.AlignmentFile(out_file_name, "wb", template=inf)

    count = 0
    passed = 0
    lengths = []

    for read in inf:
        if read.reference_id == -1 or read.is_secondary:
            continue

        count += 1
        if count % 10000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.query_name in read_set:
            outf.write(read)
            passed += 1

    print("Processed " + str(count) + " reads, written " + str(passed))
    inf.close()
    outf.close()
    pysam.index(out_file_name)


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <input.bam> <read list>")
    exit(0)
filter_reads(sys.argv[1], os.path.splitext(sys.argv[1])[0] + ".filtered.bam", read_reads(sys.argv[2]))
