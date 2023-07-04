import random
import editdistance
from numpy import histogram
from datetime import datetime



N = 10000
L = 9
ALPHABET = ["A", "C", "G", "T"]


def get_seq(length=L):
    return "".join([random.choice(ALPHABET) for _ in range(length)])


random.seed(datetime.now().timestamp())
dist = []
for i in range(N):
    s1 = get_seq()
    s2 = get_seq()
    dist.append(editdistance.eval(s1, s2))


bins = [i for i in range(L + 3)]
print(histogram(dist, bins))