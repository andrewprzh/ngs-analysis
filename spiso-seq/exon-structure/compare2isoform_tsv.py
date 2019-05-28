import os
import sys

m1 = sys.argv[1]
m2 = sys.argv[2]

d1 = {}
d2 = {}
with open(m1, "r") as infile:
    for line in infile.readlines():
        d1[line.strip().split("\t")[0].split("-")[0]] = line.strip().split("\t")[1]
with open(m2, "r") as infile:
    for line in infile.readlines():
        d2[line.strip().split("\t")[0].split("-")[0]] = line.strip().split("\t")[1]


print("Total barcodes in  " + m1 + " : " + str(len(d1)))
print("Total barcodes in  " + m2 + " : " + str(len(d2)))

common = []
for i in d2.keys():
    if i in d1.keys():
        common.append(i)
print("Common: ", len(common))

equal = []
for i in common:
    if d1[i] == d2[i]:
        equal.append(i)

print("Equal: ", len(equal))

