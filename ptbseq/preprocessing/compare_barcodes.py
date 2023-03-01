import sys
def load_dict(inf):
    bc_dict = {}
    for l in open(inf):
        v = l.strip().split('\t')
        if len(v) != 2:
            continue
        bc_dict[v[0]] = v[1]
    return bc_dict


bc_dict1 = load_dict(sys.argv[1])
bc_dict2 = load_dict(sys.argv[2])

same = 0
diff = 0
only1 = 0
only2 = 0

for k in bc_dict1.keys():
    if k not in bc_dict2:
        only1 += 1
        continue

    if bc_dict2[k] == bc_dict1[k]:
        same += 1
    else:
        diff += 1

for k in bc_dict2.keys():
    if k not in bc_dict1:
        only2 += 1

print("Set1: %d, Set2: %d" % (len(bc_dict1), len(bc_dict2)))
print("Common: %d, different: %d, only 1: %d, only 2: %d" % (same, diff, only1, only2))