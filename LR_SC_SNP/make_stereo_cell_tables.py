#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Barcode (DNB) to cell tables for one Stereo-seq sample, from the Spl-IsoFind cell segmentation and annotation.

Inputs: the barcode to cell map of the slide (SampleN_barcodeToPos.CellID.tsv[.gz], columns barcode and
CellID-original), the cell annotation (labels.csv with sample, CellID-original, first_type, spot_class, region,
subregion, x, y; or the obs of SampleN_cellbin_adjusted.h5ad saved as CSV) and the sample's barcoded reads.
Only annotated cells of --sample are used, and only barcodes that occur in the barcoded reads are written.
A cell has a type (first_type) only if it is a singlet; other annotated cells get NA.

Outputs:
  PREFIX.barcode2spot.tsv     barcode, cell_id, cell_type, region, subregion (IsoQuant --barcode2spot FILE:0:1,2,3,4);
  PREFIX.barcode_types.tsv    barcode, cell type, barcodes of typed (singlet) cells only (celltype_filter.py);
  PREFIX.barcode_to_cell.tsv  barcode, cell_id (spatial_filter.py);
  PREFIX.cell_coords.tsv      cell_id, x, y (spatial_filter.py);
  PREFIX.stats.txt            counts.
"""

import argparse
import csv
import gzip
from collections import Counter

NO_TYPE = "NA"


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


def cell_id(value):
    """CellID-original as an integer string: the barcode map writes 91689.0, the annotation 91689."""
    return str(int(float(value)))


def load_cells(fname, sample):
    cells = {}
    with open_text(fname) as f:
        for row in csv.DictReader(f):
            if row["sample"] != sample:
                continue
            cell_type = row["first_type"] if row["spot_class"] == "singlet" else NO_TYPE
            cells[cell_id(row["CellID-original"])] = (cell_type, row["region"], row["subregion"], row["x"], row["y"])
    return cells


def load_read_barcodes(fnames):
    barcodes = set()
    reads = 0
    for fname in fnames:
        for line in open_text(fname):
            if line.startswith("#"):
                continue
            reads += 1
            barcodes.add(line.split("\t", 2)[1])
    barcodes.discard("*")
    return barcodes, reads


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--barcode_to_cell", required=True, help="slide barcode to cell map (barcodeToPos.CellID.tsv)")
    parser.add_argument("--labels", required=True, help="cell annotation CSV (labels.csv)")
    parser.add_argument("--sample", required=True, help="slide name in the sample column of --labels, e.g. S1")
    parser.add_argument("--barcoded_reads", required=True, nargs="+", help="IsoQuant barcoded reads TSV(s)")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    cells = load_cells(args.labels, args.sample)
    read_barcodes, reads = load_read_barcodes(args.barcoded_reads)
    stats = Counter()
    cells_with_barcodes = set()
    with open_text(args.barcode_to_cell) as inf, open(args.output + ".barcode2spot.tsv", "w") as spot_out, \
            open(args.output + ".barcode_types.tsv", "w") as type_out, \
            open(args.output + ".barcode_to_cell.tsv", "w") as cell_out:
        spot_out.write("#barcode\tcell_id\tcell_type\tregion\tsubregion\n")
        type_out.write("#barcode\tcell_type\n")
        cell_out.write("#barcode\tcell_id\n")
        header = inf.readline().rstrip("\n").split("\t")
        barcode_col, cell_col = header.index("barcode"), header.index("CellID-original")
        for line in inf:
            v = line.rstrip("\n").split("\t")
            stats["map barcodes"] += 1
            cell = cell_id(v[cell_col])
            if cell not in cells:
                continue
            stats["map barcodes of annotated cells"] += 1
            barcode = v[barcode_col]
            if barcode not in read_barcodes:
                continue
            cell_type, region, subregion = cells[cell][:3]
            stats["written barcodes"] += 1
            cells_with_barcodes.add(cell)
            spot_out.write("%s\t%s\t%s\t%s\t%s\n" % (barcode, cell, cell_type, region, subregion))
            cell_out.write("%s\t%s\n" % (barcode, cell))
            if cell_type != NO_TYPE:
                stats["written typed barcodes"] += 1
                type_out.write("%s\t%s\n" % (barcode, cell_type))
    with open(args.output + ".cell_coords.tsv", "w") as outf:
        outf.write("#cell_id\tx\ty\n")
        for cell in sorted(cells, key=int):
            outf.write("%s\t%s\t%s\n" % (cell, cells[cell][3], cells[cell][4]))

    stats["annotated cells"] = len(cells)
    stats["typed cells"] = sum(1 for c in cells.values() if c[0] != NO_TYPE)
    stats["annotated cells with read barcodes"] = len(cells_with_barcodes)
    stats["barcoded reads"] = reads
    stats["read barcodes"] = len(read_barcodes)
    with open(args.output + ".stats.txt", "w") as outf:
        for key in ["barcoded reads", "read barcodes", "map barcodes", "map barcodes of annotated cells",
                    "written barcodes", "written typed barcodes", "annotated cells", "typed cells",
                    "annotated cells with read barcodes"]:
            outf.write("%s\t%d\n" % (key, stats[key]))
        for t, n in sorted(Counter(c[0] for c in cells.values()).items()):
            outf.write("cells of type %s\t%d\n" % (t, n))
    print(open(args.output + ".stats.txt").read(), end="")


if __name__ == "__main__":
    main()
