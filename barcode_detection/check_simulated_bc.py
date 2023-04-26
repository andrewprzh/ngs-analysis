import sys
# READ_1_ENST00000294189.11_TGACAATCTCCGTA_TACAATAGG_0_aligned_25_R_65_816_74     128     86      95      112     -       TGACAATCTCCGTA  TACAATAGG       14      True

count = 0
barcoded = 0
correct = 0
for l in open(sys.argv[1]):
    count += 1
    v = l.split("\t")
    bc = v[6]
    if bc == "*": continue
    barcoded += 1

    true_bc = v[0].split("_")[3]
    if true_bc == bc:
        correct += 1

print("Total\t%d\nBarcoded\t%d\nCorrect\t%d" % (count, barcoded, correct))
    