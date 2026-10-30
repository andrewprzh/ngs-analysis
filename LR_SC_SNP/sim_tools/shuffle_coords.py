"""Permute the (x, y) pairs among cells of a cell_coords TSV (header kept). Usage: shuffle_coords.py IN OUT SEED"""
import sys
import numpy

inp, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
lines = open(inp).read().splitlines()
header = [l for l in lines if l.startswith("#")]
rows = [l.split("\t") for l in lines if not l.startswith("#")]
perm = numpy.random.default_rng(seed).permutation(len(rows))
with open(out, "w") as f:
    f.write("\n".join(header) + "\n")
    for row, j in zip(rows, perm):
        f.write("%s\t%s\t%s\n" % (row[0], rows[j][1], rows[j][2]))
