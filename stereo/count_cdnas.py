import sys
from collections import defaultdict


count_dict = defaultdict(int)
total_count = 0
current_readid = ""
current_count = 1
total_reads = 0

for l in open(sys.argv[1]):
    v = l.strip().split('\t')
    total_count += 1

    if v[-1] == '-1' and v[1] == '*':
        count_dict[0] += 1
        total_reads += 1
        continue

    if "/ccs/" in v[0]:
        rp = v[0].split('ccs')
        read_id = rp[0] + 'ccs' + rp[1].split('_')[0]
    else:
        read_id = v[0].split('_')[0]
    if read_id == current_readid:
        if v[1] != '*':
            current_count += 1
        continue

    count_dict[current_count] += 1
    total_reads += 1
    current_readid = read_id
    current_count = 1

print(total_count, total_reads)
for k in sorted(count_dict.keys()):
    print("%d\t%d" % (k, count_dict[k]))

