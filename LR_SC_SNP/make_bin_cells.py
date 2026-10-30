#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Group square grid barcodes (e.g. Visium HD 2 um bins) into larger square bins used as pseudo-cells.

Input: TSV with barcode, x column, y row (integer grid). A barcode goes to bin (x // size, y // size);
e.g. --size 4 turns 2 um bins into 8 um bins. Use when no cell segmentation is available (simulated data).

Outputs PREFIX.barcode_to_cell.tsv (barcode, cell) and PREFIX.cell_coords.tsv (cell, x, y of the bin centre in
input grid units), the inputs of spatial_filter.py.
"""

import argparse
import gzip


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--barcode_coords", required=True, help="TSV with barcode, x, y (integer grid)")
    parser.add_argument("--size", type=int, default=4, help="bin size in grid units [%(default)s]")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    cells = set()
    with open(args.output + ".barcode_to_cell.tsv", "w") as outf:
        outf.write("#barcode\tcell\n")
        for line in open_text(args.barcode_coords):
            if line.startswith("#"):
                continue
            barcode, x, y = line.rstrip("\n").split("\t")[:3]
            bx, by = int(x) // args.size, int(y) // args.size
            cells.add((bx, by))
            outf.write("%s\tbin%d_%d_%d\n" % (barcode, args.size, bx, by))
    centre = (args.size - 1) / 2
    with open(args.output + ".cell_coords.tsv", "w") as outf:
        outf.write("#cell\tx\ty\n")
        for bx, by in sorted(cells):
            outf.write("bin%d_%d_%d\t%g\t%g\n" % (args.size, bx, by, bx * args.size + centre, by * args.size + centre))
    print("%d bins of %dx%d" % (len(cells), args.size, args.size))


if __name__ == "__main__":
    main()
