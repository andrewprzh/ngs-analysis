import os
import sys

m1 = sys.argv[1]
m2 = sys.argv[2]

exons_number =3
if m1.upper().startswith("B"):
    exons_number= 6

d1 = {}
d2 = {}
with open(m1, "r") as hagen_file:
    for line in hagen_file.readlines():
        d1[line.split("\t")[0]] = [line.split("\t")[i].strip() for i in range(1,exons_number + 1)]
with open(m2, "r") as hagen_file:
    for line in hagen_file.readlines():
        d2[line.split("\t")[0]] = [line.split("\t")[i].strip() for i in range(1,exons_number + 1)]


print("Total barcodes in  " + m1 + " : " + str(len(d1)))
print("Total barcodes in  " + m2 + " : " + str(len(d2)))

common = []
for i in d1.keys():
    if i in d2.keys():
        common.append(i)

print("Common: ", len(common))

equal = []
for i in common:
    if d1[i] == d2[i]:
        equal.append(i)

print("Equal: ", len(equal))

