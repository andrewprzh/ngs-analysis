import sys
import editdistance
from collections import defaultdict
# READ_1_ENST00000294189.11_TGACAATCTCCGTA_TACAATAGG_0_aligned_25_R_65_816_74     128     86      95      112     -       TGACAATCTCCGTA  TACAATAGG       14      True

count = 0
barcoded = 0
correct = 0
umis_dists = [0] * 15
true_umis_dists = [0] * 15


def print_dict(d):
    l = [(x, d[x]) for x in d.keys()]
    l = sorted(l, key=lambda x:x[1], reverse=True)
    for k,v in l:
        print("%s\t%d" % (k, v))


incorrectly_called = defaultdict(int)
incorrectly_detected = defaultdict(int)

for l in open(sys.argv[1]):
    if l.startswith("#"): continue
    count += 1
    v = l.split("\t")
    bc = v[1]
    umi = v[2]
    if bc == "*": continue
    barcoded += 1

    readv = v[0].split("_")
    true_bc = readv[3]
    true_umi = readv[4]
    if true_bc == bc:
        correct += 1
    else:
        incorrectly_called[bc] += 1
        incorrectly_detected[true_bc] += 1
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
print("Incorrect calls: %d" % sum(incorrectly_called.values()))
print()
print_dict(incorrectly_called)
print()
print_dict(incorrectly_detected)
print()
