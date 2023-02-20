import sys


inf = sys.argv[1]
threshold = float(sys.argv[2])

for l in open(inf):
    if l.startswith("read_id,"):
        continue
    v = l.strip().split(',')
    if len(v) < 3 or not v[2]:
        continue
    q = float(v[2])
    if q < threshold:
        continue

    print("%s\t%s" % (v[0], v[1]))