import sys

is_gene = False
for l in open(sys.argv[1]):
    if l.startswith('=') or l.startswith('E'):
        is_gene = True
        print(l.strip())
        continue

    if is_gene:
        is_gene = False
        continue

    tokens = l.strip().split()
    print('\t'.join(tokens[1:]))
