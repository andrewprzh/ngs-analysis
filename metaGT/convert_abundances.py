import os
import sys


name_dict = {}
for l in open(sys.argv[1]):
    v = l.strip().split('\t')
    if v[0] == v[1]:
        continue

    if v[0].startswith("NODE") and v[1].startswith("MED"):
        name_dict[v[0]] = v[1]
    elif v[1].startswith("NODE") and v[0].startswith("MED"):
        name_dict[v[1]] = v[0]


for l in open(sys.argv[2]):
    v = l.strip().split('\t')
    if v[0].startswith("MED"):
        print(l.strip())
    elif v[0] in name_dict:
        print(name_dict[v[0]] + "\t" + "\t".join(v[1:]))
    else:
        print(l.strip())
