import os
import sys
import pysam
from collections import defaultdict


def filter_reads(in_file_name, out_file_name_prefix, read_dict, max_iter):
    inf = pysam.AlignmentFile(in_file_name, "rb")
    out_files = []
    out_file_names = []
    for i in range(max_iter):
        out_file_name = out_file_name_prefix + ".%d.bam" % i
        out_file_names.append(out_file_name)
        out_files.append(pysam.AlignmentFile(out_file_name, "wb", template=inf))

    count = 0
    counts = defaultdict(int)

    print("Sorting reads into %d files" % len(out_files))
    for read in inf:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write("Processed " + str(count) + " reads\r")

        if read.query_name in read_dict:
            for index in read_dict[read.query_name]:
                out_files[index].write(read)
                counts[index] += 1

    print("Processed %d reads, total written %d, min %d, max %d" % (count, sum(counts.values()),
                                                                    min(counts.values()), max (counts.values())))
    inf.close()
    for outf in out_files:
        outf.close()
    print("Done. Indexing now.")
    for out_file_name in out_file_names:
        pysam.index(out_file_name)


def load_lists(indir, max_iter):
    print("Loading read lists")
    read_dict = defaultdict(list)
    for i in range(max_iter):
        for l in open(os.path.join(indir, "Run%d/sampled_AllInfo.tsv" % i)):
            read_dict[l.strip().split("\t")[0]].append(i)
    print("Loaded %d reads" % len(read_dict))
    return read_dict



if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <input.bam> <runs folder> <output folder> [max=100]")
    exit(0)

if len(sys.argv) == 5:
    max_iter = int(sys.argv[4])
else:
    max_iter = 100

if not os.path.exists(sys.argv[3]):
    os.makedirs(sys.argv[3])
filter_reads(sys.argv[1], os.path.join(sys.argv[3], os.path.splitext(sys.argv[1])[0]), load_lists(sys.argv[2], max_iter), max_iter)
