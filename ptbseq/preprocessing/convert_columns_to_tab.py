import sys

for l in open(sys.argv[1]):
    v = l.strip().split('\t')
    print("\t".join(v[:4] + v[4].split('::')))
