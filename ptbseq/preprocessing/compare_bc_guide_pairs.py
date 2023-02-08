import sys


def load_pairs(inf):
    pair_set = set()
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 2:
            continue
        pair_set.add((v[0], v[1].replace("GENE_", "")))
    return pair_set


set1 = load_pairs(sys.argv[1])
set2 = load_pairs(sys.argv[2])
print("Set1: %d, Set2: %d" % (len(set1), len(set2)))
print("Common: %d, only 1: %d, only 2: %d" % (len(set1.intersection(set2)), len(set1.difference(set2)), len(set2.difference(set1))))