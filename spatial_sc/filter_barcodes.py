import sys


def filter_barcodes(inf):
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) < 10 or v[6] == "*": # or v[9] != "True"
            continue
        sys.stdout.write(l)


filter_barcodes(sys.argv[1])