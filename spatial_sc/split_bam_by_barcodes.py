import os
import sys
import pysam
from collections import defaultdict


def read_called_barcodes(inf):
    barcoded_reads = {}
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) < 7:
            continue
        bc = v[6]
        if bc == "*":
            continue
        barcoded_reads[v[0]] = v[6]
    print("Loaded %d read ids" % len(barcoded_reads))
    return barcoded_reads


def filter_reads(in_file_name, out_file_prefix, barcoded_reads, min_reads = 100):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    barcode_to_reads = defaultdict(list)
    count = 0
    passed = 0
    unaligned = 0
    unbarcoded = 0

    for read in inf:
        if read.reference_id == -1:
            unaligned += 1
            continue

        count += 1
        if count % 10000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.query_name not in barcoded_reads:
            unbarcoded += 1
            continue

        barcode_to_reads[barcoded_reads[read.query_name]].append(read)
        passed += 1

    print("Processed %d reads, unaligned %d" % (count, unaligned))
    print("Among aligned reads, %d are not barcoded, %d were collected" % (unbarcoded, passed))

    bc_count = 0
    for bc in barcode_to_reads.keys():
        if len(barcode_to_reads[bc]) < min_reads:
            continue
        out_file_name = out_file_prefix + bc + ".bam"
        outf = pysam.AlignmentFile(out_file_name, "wb", template=inf)
        bc_count += 1
        for read in barcode_to_reads[bc]:
            outf.write(read)

        outf.close()
        pysam.index(out_file_name)
    print("Created " + str(bc_count) + " files")
    inf.close()


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <input.bam> <read list> [prefix]")
    exit(0)

if len(sys.argv) == 4:
    prefix = sys.argv[3]
else:
    prefix = ""

filter_reads(sys.argv[1], prefix, read_called_barcodes(sys.argv[2]))
