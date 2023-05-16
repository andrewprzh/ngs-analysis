import os
import sys
import pysam
import numpy


def read_reads(inf):
    read_set = set()
    for l in open(inf):
        read_set.add(l.strip().split('\t')[0])
    return read_set


def read_called_barcodes(inf):
    read_set = set()
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) >= 7 and v[6] == "*":
            continue
        read_set.add(v[0])
    print("Loaded %d read ids" % len(read_set))
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
    print("Usage: " + sys.argv[0] + " <input.bam> <read list> [suffix]")
    exit(0)

if len(sys.argv) == 4:
    suffix = "." + sys.argv[3]
else:
    suffix = ""
filter_reads(sys.argv[1], os.path.splitext(sys.argv[1])[0] + ".filtered" + suffix + ".bam", read_called_barcodes(sys.argv[2]))
