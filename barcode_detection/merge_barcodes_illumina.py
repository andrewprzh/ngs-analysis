import sys


inf1 = open(sys.argv[1])
inf2 = open(sys.argv[2])
while True:
    l1 = inf1.readline()
    l2 = inf2.readline()
    v1 = l1.strip().split('\t')
    v2 = l2.strip().split('\t')
    read1 = v1[0]
    read2 = v2[0]
    assert read2 == read1
    bc1 = v1[6]
    umi1 = v1[9]
    bc2 = v2[6]
    umi2 = v2[9]

    if bc2 == bc1 == "*":
        continue
    elif bc1 == "*":
        sys.stdout.write(l2)
    elif bc2 == "*":
        sys.stdout.write(l1)
    else:
        sys.stderr.write("Both reads have the barcode")
        continue







