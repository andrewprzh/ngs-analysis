import sys
from collections import defaultdict

subs_index = {'A':0, 'C':1, 'G':2, 'T':3}

sub_dict = defaultdict(int)
sub_events = defaultdict(int)
for l in open(sys.argv[1]):
    v = l.strip().split('\t')
    if v[7] == 'AllSubs': continue
    gsubs = set(v[12].split(' ')) if v[12] != '-' else set()
    for s in v[7].split(' '):
        if s == '-': continue
        if s in gsubs: continue
        c = int(v[6][1:-1].split(',')[subs_index[s[1]]])
        sub_dict[s] += c
        sub_events[s] += 1

print(", ".join(("%s: %d" % (k, sub_dict[k])) for k in sorted(sub_dict.keys())))
print(", ".join(("%s: %d" % (k, sub_events[k])) for k in sorted(sub_events.keys())))
