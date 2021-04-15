import sys
import glob
from collections import  defaultdict
import os

def print_stats(consistent_dict, inconsistent_dict):
    qrange = (10, 20)
    print("shift\ttotal\t" + "\t".join([str(x) for x in range(qrange[0], qrange[1] + 1)]))
    for shift in range(-10, 11):
        total_cons = sum(consistent_dict[shift].values())
        total_inc = sum(inconsistent_dict[shift].values())
        total = total_inc + total_cons
        total_rate = 0 if total == 0 else 100 * total_inc / total
        rates = [total_rate]
        for q in range(qrange[0], qrange[1] + 1):
            cons = consistent_dict[shift][q]
            inc = inconsistent_dict[shift][q]
            total = cons + inc
            rates.append(0 if total == 0 else 100 * inc / cons)
        print("%d\t%s" % (shift, "\t".join([str(x) for x in rates])))


sys.stderr.write("Loading qualities\n")
current_id = ""
read_qualities = defaultdict(float)
for l in open(sys.argv[2]):
    if l.startswith(">"):
        current_id = l.strip()[1:]
    else:
        read_qualities[current_id] = float(l.strip())

consistent_acc_dict = defaultdict(lambda: defaultdict(int))
inconsistent_acc_dict = defaultdict(lambda: defaultdict(int))
consistent_donor_dict = defaultdict(lambda: defaultdict(int))
inconsistent_donor_dict = defaultdict(lambda: defaultdict(int))

for f in glob.glob(sys.argv[1] + "*.tsv"):
    sys.stderr.write("Processing %s\n" % f)
    for l in open(f):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        read_id = t[0]
        read_type = t[5]
        read_q = int(round(read_qualities[read_id]))
        if read_type == "consistent":
            donor_up = -int(t[6])
            donor_down = int(t[7])
            if donor_up != 0 and donor_down != 0:
                consistent_donor_dict[donor_up][read_q] += 1
                consistent_donor_dict[donor_down][read_q] += 1
            elif donor_up != 0:
                consistent_donor_dict[donor_up][read_q] += 2
            elif donor_down != 0:
                consistent_donor_dict[donor_down][read_q] += 2
            acc_up = -int(t[8])
            acc_down = int(t[9])
            if acc_up != 0 and acc_down != 0:
                consistent_acc_dict[acc_up][read_q] += 1
                consistent_acc_dict[acc_down][read_q] += 1
            elif acc_up != 0:
                consistent_acc_dict[acc_up][read_q] += 2
            elif acc_down != 0:
                consistent_acc_dict[acc_down][read_q] += 2
        elif read_type == "incosistent":
            donor_diff = int(t[10])
            acc_diff = int(t[11])
            if acc_diff != 0:
                inconsistent_acc_dict[acc_diff][read_q] += 2
            if donor_diff != 0:
                inconsistent_donor_dict[donor_diff][read_q] += 2

sys.stderr.write("Outputting stats\n")

print("donor stats")
print_stats(consistent_donor_dict, inconsistent_donor_dict)
print("acceptor stats")
print_stats(consistent_acc_dict, inconsistent_acc_dict)
