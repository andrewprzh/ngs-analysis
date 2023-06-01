import sys


inf1 = open(sys.argv[1])
inf2 = open(sys.argv[2])
while True:
    l1 = inf1.readline()
    l2 = inf2.readline()
    if not l1 or not l2:
        break
    v1 = l1.strip().split('\t')
    v2 = l2.strip().split('\t')
    if len(v1) < 10 or len(v2) < 10:
        continue
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
        if abs(len(umi1) - 9) < abs(len(umi2) - 9):
            sys.stdout.write(l1)
        else:
            sys.stdout.write(l2)







