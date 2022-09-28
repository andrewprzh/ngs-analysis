import sys

id_set = set()
for l in open(sys.argv[1]):
    if l.startswith("target"):
        continue

    v = l.strip().split()
    id_set.add(v[0])


for l in open(sys.argv[2]):
    if l.startswith("#") or l.startswith("Gene"):
        continue
    v = l.strip().split()
    if v[0] in id_set:
        print(l.strip())