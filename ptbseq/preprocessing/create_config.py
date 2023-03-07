import sys

guides = set()
negative_control = set()
intergenic_control = set()
for l in open(sys.argv[1]):
    v = l.strip().split("\t")
    c = v[1].split("::")
    if c[1] == "Negative_control":
        negative_control.add(c[2])
    elif c[1] == "Splicing_factor":
        guides.add(c[2])
    elif c[1] == "Intergenic_control":
        intergenic_control.add(c[2])

for g in guides:
    print("%s\t%s" % (g, ",".join(list(sorted(negative_control)))))
for g in guides:
    print("%s\t%s" % (g, ",".join(list(sorted(intergenic_control)))))

