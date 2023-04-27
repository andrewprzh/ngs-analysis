import sys
import editdistance
# READ_1_ENST00000294189.11_TGACAATCTCCGTA_TACAATAGG_0_aligned_25_R_65_816_74     128     86      95      112     -       TGACAATCTCCGTA  TACAATAGG       14      True

count = 0
barcoded = 0
correct = 0
umis_dists = [0] * 15
true_umis_dists = [0] * 15

for l in open(sys.argv[1]):
    count += 1
    v = l.split("\t")
    bc = v[6]
    umi = v[7]
    if bc == "*": continue
    barcoded += 1

    readv = v[0].split("_")
    true_bc = readv[3]
    true_umi = readv[4]
    if true_bc == bc:
        correct += 1
    if len(umi) > 12:
        umi = umi[:13]
    umi_ed = editdistance.eval(true_umi, umi)
    umis_dists[umi_ed] += 1
    if v[9].startswith("T"):
        true_umis_dists[umi_ed] += 1

print("Total\t%d\nBarcode\t%d\nCorrect\t%d" % (count, barcoded, correct))
print("UMIs within 1\t%d\nUMIs within 2\t%d\nUMIs within 3\t%d" % (sum(umis_dists[:2]), sum(umis_dists[:3]), sum(umis_dists[:4])))
print(true_umis_dists)
print("UMIs within 1\t%d\nUMIs within 2\t%d\nUMIs within 3\t%d" % (sum(true_umis_dists[:2]), sum(true_umis_dists[:3]), sum(true_umis_dists[:4])))
print(true_umis_dists)
print()
