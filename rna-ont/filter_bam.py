import os
import sys
import pysam
import numpy


def filter_reads(in_file_name, out_file_name, min_aligned_len):
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

        blocks = sorted(read.get_blocks())
        total_blocks_len = 0
        for block in blocks:
            total_blocks_len += block[1] - block[0]

        lengths.append(total_blocks_len)
        if total_blocks_len >= min_aligned_len:
            outf.write(read)
            passed += 1

    print("Processed " + str(count) + " reads, written " + str(passed))
    inf.close()
    outf.close()
    v, k = numpy.histogram(lengths, bins = [500 * i for i in range(10)] + [max(6000, max(lengths))])
    for i in range(len(v)):
        print(str(k[i]) + '\t' + str(v[i]))

    pysam.index(out_file_name)


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <input.bam> <aligned length>")
    exit(0)
filter_reads(sys.argv[1], os.path.splitext(sys.argv[1])[0] + ".filtered.bam", int(sys.argv[2]))
