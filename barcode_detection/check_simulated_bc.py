import sys
import editdistance
from collections import defaultdict
from numpy import histogram
# READ_1_ENST00000294189.11_TGACAATCTCCGTA_TACAATAGG_0_aligned_25_R_65_816_74     128     86      95      112     -       TGACAATCTCCGTA  TACAATAGG       14      True

count = 0
barcoded = 0
correct = 0
umis_dists = [0] * 15
true_umis_dists = [0] * 15

score_barcoded = defaultdict(int)
score_correct = defaultdict(int)

SCORES = [11, 12, 13, 14]
MIN_SCORE = 13

barcode_barcoded = defaultdict(int)
barcode_correct = defaultdict(int)

true_barcodes = set()
all_barcodes = set()


def print_dict(d):
    l = [(x, d[x]) for x in d.keys()]
    l = sorted(l, key=lambda x:x[1], reverse=True)
    for k,v in l:
        print("%s\t%d" % (k, v))


incorrectly_called = defaultdict(int)
incorrectly_detected = defaultdict(int)
basic_stats = defaultdict(int)

# read_id        barcode UMI     BC_score        valid_UMI       strand  polyT_start     primer_end      linker_start    linker_end
for l in open(sys.argv[1]):
    if l.startswith("#"): continue

    v = l.split("\t")
    bc = v[1]
    umi = v[2]
    score = int(v[3])

    count += 1
    if v[6] != '-1':
        basic_stats["polyA"] += 1
    if v[7] != '-1':
        basic_stats["primer"] += 1
    if v[8] != '-1':
        basic_stats["linker"] += 1
    if bc == "*":
        continue
    else:
        basic_stats["barcode"] += 1
    barcoded += 1
    all_barcodes.add(bc)
    score_barcoded[score] += 1
    if score >= MIN_SCORE:
        barcode_barcoded[bc] += 1

    readv = v[0].split("_")
    true_bc = readv[3]
    true_umi = readv[4]
    true_barcodes.add(true_bc)

    if v[6] != '-1':
        basic_stats["polyA"] += 1
    if v[7] != '-1':
        basic_stats["primer"] += 1
    if v[8] != '-1':
        basic_stats["linker"] += 1

    if true_bc == bc:
        correct += 1
        score_correct[score] += 1
        if score >= MIN_SCORE:
            barcode_correct[bc] += 1
    else:
        incorrectly_called[bc] += 1
        incorrectly_detected[true_bc] += 1
    if len(umi) > 12:
        umi = umi[:13]
    umi_ed = editdistance.eval(true_umi, umi)
    umis_dists[umi_ed] += 1
    if v[9].startswith("T"):
        true_umis_dists[umi_ed] += 1

print("All barcodes %d, true barcodes %d" % (len(all_barcodes), len(true_barcodes)))

print("Total\t%d\nBarcode\t%d\nCorrect\t%d" % (count, barcoded, correct))

print("MinScore\tPrecision\tRecall")
for score in SCORES:
    total_barcoded = sum(score_barcoded[s] for s in range(score, 15))
    total_correct = sum(score_correct[s] for s in range(score, 15))
    print("%d\t%.2f\t%.2f" % (score, (100 * total_correct / total_barcoded), (100 * total_correct / count)))

print("UMIs within 1\t%d\nUMIs within 2\t%d\nUMIs within 3\t%d" % (sum(umis_dists[:2]), sum(umis_dists[:3]), sum(umis_dists[:4])))
print(true_umis_dists)
print("UMIs within 1\t%d\nUMIs within 2\t%d\nUMIs within 3\t%d" % (sum(true_umis_dists[:2]), sum(true_umis_dists[:3]), sum(true_umis_dists[:4])))
print(true_umis_dists)
print("Incorrect calls: %d" % sum(incorrectly_called.values()))
print()
print("Basic stats")
print_dict(basic_stats)

per_barcode_precision_values = []
for bc in barcode_barcoded:
    if bc not in true_barcodes:
        continue
    correct = barcode_correct[bc] if bc in barcode_correct else 0
    precision = 100 * correct / barcode_barcoded[bc]
    per_barcode_precision_values.append((bc, precision))

precision_arr = [x[1] for x in per_barcode_precision_values]
h = histogram(precision_arr, [5 * x for x in range(21)])
print(h)

low_prec_barcodes = list(sorted(filter(lambda x: x[1] <= 80, per_barcode_precision_values), key=lambda x: x[1]))

print(low_prec_barcodes)

#print_dict(incorrectly_called)
print()
#print_dict(incorrectly_detected)
print()
