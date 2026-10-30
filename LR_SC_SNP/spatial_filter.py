#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Phase 2 spatial filter: is the alt fraction of a variant spatially autocorrelated?

Barcodes (spots / bins) are aggregated to cells with --barcode_to_cell: alt and depth UMIs are summed per cell, and
every cell covering a variant gets the value alt / depth; cells that do not cover it are left out.
Moran's I and its permutation p-value are computed by Spl-IsoFind (SplIsoFind.sv.moransI_sparse, all cells together)
with a k-nearest-neighbour graph over the covering cells; q-values are Benjamini-Yekutieli adjusted by Spl-IsoFind.
Spl-IsoFind skips a variant (reports NaN) when fewer than --mincells cells cover it, when its values do not vary,
or when fewer than --mincellspergroup cells (or a fraction below --imb) have alt fraction > 0.5 or <= 0.5.
Spl-IsoFind parameters that are not given keep the package defaults.

Cell coordinates come from --cell_coords (cell, x, y) or are the mean of barcode coordinates (--barcode_coords).

Output: PREFIX.spatial_sites.tsv, one line per variant.
Run in an environment with SplIsoFind (pip install SplIsoFind).
"""

import argparse
import gzip
import sys
from collections import defaultdict

import numpy
import pandas
from scipy.sparse import csr_matrix

import SplIsoFind


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


def load_two_columns(fname):
    result = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if len(v) >= 2:
            result[v[0]] = v[1]
    return result


def load_coords(fname):
    coords = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if len(v) >= 3:
            coords[v[0]] = (float(v[1]), float(v[2]))
    return coords


def cell_coordinates(args, barcode_to_cell):
    if args.cell_coords:
        return load_coords(args.cell_coords)
    sums = defaultdict(lambda: [0.0, 0.0, 0])
    for line in open_text(args.barcode_coords):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        cell = barcode_to_cell.get(v[0])
        if cell is None:
            continue
        s = sums[cell]
        s[0] += float(v[1])
        s[1] += float(v[2])
        s[2] += 1
    return {cell: (x / n, y / n) for cell, (x, y, n) in sums.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="PREFIX.sites.tsv from call_variants.py")
    parser.add_argument("--barcode_counts", required=True, help="PREFIX.barcode_counts.tsv.gz from call_variants.py")
    parser.add_argument("--barcode_to_cell", required=True, help="TSV with barcode and cell id")
    coords = parser.add_mutually_exclusive_group(required=True)
    coords.add_argument("--cell_coords", help="TSV with cell id, x, y")
    coords.add_argument("--barcode_coords", help="TSV with barcode, x, y; a cell is placed at the mean of its barcodes")
    parser.add_argument("--fdr", type=float, default=0.05, help="q-value cut-off for the pass column [%(default)s]")
    parser.add_argument("--k", type=int, help="nearest neighbours in the spatial graph (Spl-IsoFind default)")
    parser.add_argument("--nperm", type=int, help="permutations (Spl-IsoFind default)")
    parser.add_argument("--mincells", type=int, help="minimal number of covering cells (Spl-IsoFind default)")
    parser.add_argument("--mincellspergroup", type=int,
                        help="minimal number of cells with alt fraction > 0.5 and <= 0.5 (Spl-IsoFind default)")
    parser.add_argument("--imb", type=float, help="minimal fraction of the minority group (Spl-IsoFind default)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    barcode_to_cell = load_two_columns(args.barcode_to_cell)
    coords = cell_coordinates(args, barcode_to_cell)
    cells = sorted(coords)
    cell_index = {c: i for i, c in enumerate(cells)}
    print("Cells with coordinates: %d" % len(cells), file=sys.stderr)

    variant_ids = [line.rstrip("\n").split("\t")[2] for line in open(args.sites) if not line.startswith("#")]
    variant_index = {v: i for i, v in enumerate(variant_ids)}
    # (variant, cell) -> [alt, depth]
    counts = defaultdict(lambda: [0, 0])
    unmapped = 0
    for line in open_text(args.barcode_counts):
        if line.startswith("#"):
            continue
        vid, barcode, alt, depth = line.rstrip("\n").split("\t")
        cell = barcode_to_cell.get(barcode)
        if cell is None or cell not in cell_index or vid not in variant_index:
            unmapped += int(depth)
            continue
        c = counts[(variant_index[vid], cell_index[cell])]
        c[0] += int(alt)
        c[1] += int(depth)
    print("UMIs in barcodes without a cell: %d" % unmapped, file=sys.stderr)

    rows = numpy.fromiter((cell for (_, cell) in counts), dtype=numpy.int64, count=len(counts))
    cols = numpy.fromiter((v for (v, _) in counts), dtype=numpy.int64, count=len(counts))
    values = numpy.fromiter((a / d for a, d in counts.values()), dtype=float, count=len(counts))
    # explicit zeros mark covering cells without alt UMIs; Spl-IsoFind takes every stored value
    x_sparse = csr_matrix((values, (rows, cols)), shape=(len(cells), len(variant_ids)))
    covered = numpy.bincount(cols, minlength=len(variant_ids))
    alt_cells = numpy.bincount(cols[values > 0.5], minlength=len(variant_ids))

    labels = pandas.DataFrame({"x": [coords[c][0] for c in cells], "y": [coords[c][1] for c in cells]},
                              index=cells)
    ids = pandas.DataFrame({"Transcript ID": variant_ids})
    params = {name: getattr(args, name) for name in ["k", "nperm", "mincells", "mincellspergroup", "imb"]
              if getattr(args, name) is not None}
    m_i, pval, qval = SplIsoFind.sv.moransI_sparse(x_sparse, labels, ids, celltypes=["All"], n_jobs=args.threads,
                                                   seed=args.seed, **params)

    passed = 0
    with open(args.output + ".spatial_sites.tsv", "w") as outf:
        outf.write("#variant_id\tcovered_cells\tcells_alt_fraction_gt_0.5\tmoransI\tpvalue\tqvalue\tpass\n")
        for i, vid in enumerate(variant_ids):
            mi, p, q = m_i["All"].iloc[i], pval["All"].iloc[i], qval["All"].iloc[i]
            is_pass = int(not numpy.isnan(q) and q < args.fdr)
            passed += is_pass
            outf.write("%s\t%d\t%d\t%s\t%s\t%s\t%d\n" % (vid, covered[i], alt_cells[i],
                                                         "NA" if numpy.isnan(mi) else "%.4f" % mi,
                                                         "NA" if numpy.isnan(p) else "%.3g" % p,
                                                         "NA" if numpy.isnan(q) else "%.3g" % q, is_pass))
    tested = int(pval["All"].notna().sum())
    print("Variants: %d, tested: %d, passing FDR %g: %d" % (len(variant_ids), tested, args.fdr, passed))


if __name__ == "__main__":
    main()
