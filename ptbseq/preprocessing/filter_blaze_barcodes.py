import sys


inf = sys.argv[1]
threshold = float(sys.argv[2])
whitelist = sys.argv[3]

whitelist_set = set()
for l in open(whitelist):
    whitelist_set.add(l.split()[0])

for l in open(inf):
    if l.startswith("read_id,"):
        continue
    v = l.strip().split(',')
    if len(v) < 3 or not v[2]:
        continue
    q = float(v[2])
    if v[1] not in whitelist_set or q < threshold:
        continue

    print("%s\t%s" % (v[0], v[1]))