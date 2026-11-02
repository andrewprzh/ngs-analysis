#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Phase 2 spatial filter: is the alt fraction of a variant spatially autocorrelated?

Barcodes (spots / bins) are aggregated to cells with --barcode_to_cell: alt and depth UMIs are summed per cell, and
every cell covering a variant gets the value alt / depth; cells that do not cover it are left out.
A variant is tested when it passes the filters of Spl-IsoFind (SplIsoFind.sv.moransI_sparse): at least --mincells
covering cells, values that vary and, when cells lie on both sides of 0.5 (alt fraction > 0.5 or <= 0.5), at least
--mincellspergroup cells and a fraction --imb on each side. Spl-IsoFind parameters that are not given keep the
package defaults.

Moran's I is computed by esda over Spl-IsoFind's k-nearest-neighbour graph of the covering cells, with an analytic
one-sided p-value: the z-score of I under randomisation (values shuffled over the cells), so there is no
permutation floor. P-values are Benjamini-Hochberg adjusted over the tested variants. Variants passing --fdr are
checked by Spl-IsoFind with --check_perm permutations: their permutation p-values replace the analytic ones, BH is
recomputed, and a variant passes when it was checked and its final q-value is below --fdr.
With --restrict_to only the variants also called in another call set are tested.

Cell coordinates come from --cell_coords (cell, x, y) or are the mean of barcode coordinates (--barcode_coords).
The barcode lines of each variant must be consecutive in --barcode_counts, as call_variants.py writes them.

Output: PREFIX.spatial_sites.tsv, one line per variant.
Run in an environment with SplIsoFind (pip install SplIsoFind).
"""

import argparse
import inspect
import sys
from collections import defaultdict
from multiprocessing import Pool

import esda
import numpy
import pandas
from scipy.sparse import csc_matrix
from scipy.stats import norm

import SplIsoFind
from SplIsoFind.spatially_variable import _calculate_weight_matrix_sklearn

from celltype_filter import add_restrict_argument, benjamini_hochberg, load_variant_ids, open_text

# Spl-IsoFind draws a (notebook) progress bar for every call, one variant per call here
SplIsoFind.spatially_variable.tqdm = lambda iterable, **kwargs: iterable

SPLISOFIND_PARAMS = ["k", "mincells", "mincellspergroup", "imb"]
CHUNK_LINES = 5000000

# cell coordinates as Spl-IsoFind labels and as an array, set before the worker pool forks
LABELS = None
COORDS = None


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


class CellValues:
    """Per-variant cell values from barcode counts, kept only for the variants passing Spl-IsoFind's filters."""
    def __init__(self, tested, params):
        self.tested = tested  # bool per variant of the sites file
        self.params = params
        self.covered = numpy.zeros(len(tested), dtype=numpy.int64)
        self.alt_cells = numpy.zeros(len(tested), dtype=numpy.int64)
        self.seen = numpy.zeros(len(tested), dtype=bool)
        self.eligible = {}  # variant -> (cell indices, values)
        self.unmapped_umis = 0

    def add(self, v, c, alt, depth):
        """Variant, cell (-1: none), alt and depth per barcode line, holding every line of each variant in v."""
        blocks = v[numpy.r_[True, v[1:] != v[:-1]]]
        blocks = blocks[blocks >= 0]
        if len(numpy.unique(blocks)) != len(blocks) or self.seen[blocks].any():
            sys.exit("Barcode lines of a variant are not consecutive in the barcode counts file")
        self.seen[blocks] = True
        self.unmapped_umis += int(depth[(v >= 0) & (c < 0)].sum())
        keep = (v >= 0) & (c >= 0)
        keep[keep] = self.tested[v[keep]]
        v, c, alt, depth = v[keep], c[keep], alt[keep], depth[keep]
        if not len(v):
            return
        # barcodes of the same cell are summed; cells are in ascending order within a variant
        order = numpy.lexsort((c, v))
        v, c, alt, depth = v[order], c[order], alt[order], depth[order]
        starts = numpy.flatnonzero(numpy.r_[True, (v[1:] != v[:-1]) | (c[1:] != c[:-1])])
        v, c = v[starts], c[starts]
        values = numpy.add.reduceat(alt, starts) / numpy.add.reduceat(depth, starts)
        bounds = numpy.flatnonzero(numpy.r_[True, v[1:] != v[:-1]])
        ends = numpy.r_[bounds[1:], len(v)]
        variants = v[bounds]
        n = ends - bounds
        high = numpy.add.reduceat((values > 0.5).astype(numpy.int64), bounds)
        varies = numpy.maximum.reduceat(values, bounds) > numpy.minimum.reduceat(values, bounds)
        self.covered[variants] = n
        self.alt_cells[variants] = high
        # Spl-IsoFind checks the groups x > 0.5 and x <= 0.5 only when both are present
        smaller = numpy.where((high > 0) & (high < n), numpy.minimum(high, n - high), n)
        ok = (n >= self.params["mincells"]) & varies & \
            (smaller >= self.params["mincellspergroup"]) & (smaller / n >= self.params["imb"])
        for i in numpy.flatnonzero(ok).tolist():
            self.eligible[int(variants[i])] = (c[bounds[i]:ends[i]].astype(numpy.int32), values[bounds[i]:ends[i]])


def read_cell_values(fname, barcode_index, barcode_cell, variant_index, cell_values):
    """Stream the barcode counts in chunks; the lines of the last variant of a chunk wait for the next chunk."""
    carry = None
    reader = pandas.read_csv(fname, sep="\t", skiprows=1, header=None, names=["v", "b", "a", "d"],
                             dtype={"v": str, "b": str, "a": numpy.int64, "d": numpy.int64}, chunksize=CHUNK_LINES)
    for chunk in reader:
        b = barcode_index.get_indexer(chunk["b"])
        lines = (variant_index.get_indexer(chunk["v"]), numpy.where(b >= 0, barcode_cell[b], -1),
                 chunk["a"].to_numpy(), chunk["d"].to_numpy())
        if carry is not None:
            lines = tuple(numpy.concatenate([x, y]) for x, y in zip(carry, lines))
        v = lines[0]
        other = numpy.flatnonzero(v != v[-1])
        last = other[-1] + 1 if len(other) else 0
        carry = tuple(x[last:] for x in lines)
        cell_values.add(*(x[:last] for x in lines))
    if carry is not None:
        cell_values.add(*carry)


def analytic_moran(task):
    """Moran's I over Spl-IsoFind's kNN weights, its randomisation z-score and one-sided p-value."""
    v, cells, values, k = task
    w = _calculate_weight_matrix_sklearn(COORDS[cells], k)
    res = esda.Moran(values, w, transformation="b", permutations=0)
    return v, res.I, res.z_rand, float(norm.sf(res.z_rand))


def permutation_moran(task):
    """Spl-IsoFind's permutation p-value of one variant, seeded as in a single run over all variants."""
    v, vid, cells, values, seed, params = task
    x = csc_matrix((values, (numpy.arange(len(cells)), numpy.zeros(len(cells), dtype=numpy.int64))),
                   shape=(len(cells), 1))
    _, pval, _ = SplIsoFind.sv.moransI_sparse(x, LABELS.iloc[cells], pandas.DataFrame({"Transcript ID": [vid]}),
                                              celltypes=["All"], n_jobs=1, seed=seed, **params)
    return v, float(pval["All"].iloc[0])


def fmt(x, form):
    return "NA" if numpy.isnan(x) else form % x


def main():
    global LABELS, COORDS
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="PREFIX.sites.tsv from call_variants.py")
    parser.add_argument("--barcode_counts", required=True, help="PREFIX.barcode_counts.tsv.gz from call_variants.py")
    parser.add_argument("--barcode_to_cell", required=True, help="TSV with barcode and cell id")
    coords = parser.add_mutually_exclusive_group(required=True)
    coords.add_argument("--cell_coords", help="TSV with cell id, x, y")
    coords.add_argument("--barcode_coords", help="TSV with barcode, x, y; a cell is placed at the mean of its barcodes")
    add_restrict_argument(parser)
    parser.add_argument("--fdr", type=float, default=0.05, help="BH-adjusted p-value cut-off [%(default)s]")
    parser.add_argument("--check_perm", type=int, default=1000000,
                        help="permutations checking the variants that pass the analytic test [%(default)s]")
    parser.add_argument("--k", type=int, help="nearest neighbours in the spatial graph (Spl-IsoFind default)")
    parser.add_argument("--mincells", type=int, help="minimal number of covering cells (Spl-IsoFind default)")
    parser.add_argument("--mincellspergroup", type=int,
                        help="minimal number of cells with alt fraction > 0.5 and <= 0.5 (Spl-IsoFind default)")
    parser.add_argument("--imb", type=float, help="minimal fraction of the minority group (Spl-IsoFind default)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    defaults = inspect.signature(SplIsoFind.sv.moransI_sparse).parameters
    params = {name: getattr(args, name) if getattr(args, name) is not None else defaults[name].default
              for name in SPLISOFIND_PARAMS}

    barcode_to_cell = load_two_columns(args.barcode_to_cell)
    coords = cell_coordinates(args, barcode_to_cell)
    cells = sorted(coords)
    cell_index = {c: i for i, c in enumerate(cells)}
    print("Cells with coordinates: %d" % len(cells), file=sys.stderr)
    LABELS = pandas.DataFrame({"x": [coords[c][0] for c in cells], "y": [coords[c][1] for c in cells]}, index=cells)
    COORDS = LABELS[["x", "y"]].to_numpy()
    barcodes = [b for b, c in barcode_to_cell.items() if c in cell_index]
    barcode_cell = numpy.array([cell_index[barcode_to_cell[b]] for b in barcodes], dtype=numpy.int64)
    barcode_index = pandas.Index(barcodes)
    del barcode_to_cell, coords, barcodes

    # seeds follow the variant index in the sites file, so restricted runs give the same permutation p-values
    all_ids = load_variant_ids(args.sites)
    variant_ids = load_variant_ids(args.sites, args.restrict_to)
    variant_index = pandas.Index(all_ids)
    tested = numpy.zeros(len(all_ids), dtype=bool)
    tested[variant_index.get_indexer(variant_ids)] = True
    cell_values = CellValues(tested, params)
    read_cell_values(args.barcode_counts, barcode_index, barcode_cell, variant_index, cell_values)
    print("UMIs in barcodes without a cell: %d" % cell_values.unmapped_umis, file=sys.stderr)
    eligible = sorted(cell_values.eligible)
    print("Variants: %d, passing Spl-IsoFind's filters: %d" % (len(variant_ids), len(eligible)), file=sys.stderr)

    n = len(all_ids)
    moran_i, z_scores, pvalues, perm_pvalues = (numpy.full(n, numpy.nan) for _ in range(4))
    with Pool(args.threads) as pool:
        tasks = ((v,) + cell_values.eligible[v] + (params["k"],) for v in eligible)
        for v, mi, z, p in pool.imap_unordered(analytic_moran, tasks, chunksize=16):
            moran_i[v], z_scores[v], pvalues[v] = mi, z, p
        qvalues = numpy.full(n, numpy.nan)
        qvalues[eligible] = benjamini_hochberg(pvalues[eligible])
        checked = [v for v in eligible if qvalues[v] < args.fdr]
        print("Passing the analytic test: %d; checking them with %d permutations" % (len(checked), args.check_perm),
              file=sys.stderr)
        perm_params = dict(params, nperm=args.check_perm)
        tasks = ((v, all_ids[v]) + cell_values.eligible[v] + (args.seed + v, perm_params) for v in checked)
        for v, p in pool.imap_unordered(permutation_moran, tasks):
            perm_pvalues[v] = p
    final_p = numpy.where(numpy.isnan(perm_pvalues), pvalues, perm_pvalues)
    final_q = numpy.full(n, numpy.nan)
    final_q[eligible] = benjamini_hochberg(final_p[eligible])

    passed = 0
    with open(args.output + ".spatial_sites.tsv", "w") as outf:
        outf.write("#variant_id\tcovered_cells\tcells_alt_fraction_gt_0.5\tmoransI\tz_score\tpvalue\tqvalue"
                   "\tperm_pvalue\tfinal_qvalue\tpass\n")
        for vid in variant_ids:
            v = variant_index.get_loc(vid)
            is_pass = int(not numpy.isnan(perm_pvalues[v]) and final_q[v] < args.fdr)
            passed += is_pass
            outf.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % (
                vid, cell_values.covered[v], cell_values.alt_cells[v], fmt(moran_i[v], "%.4f"),
                fmt(z_scores[v], "%.3f"), fmt(pvalues[v], "%.3g"), fmt(qvalues[v], "%.3g"),
                fmt(perm_pvalues[v], "%.3g"), fmt(final_q[v], "%.3g"), is_pass))
    print("Variants: %d, tested: %d, passing FDR %g analytically: %d, after the permutation check: %d" %
          (len(variant_ids), len(eligible), args.fdr, len(checked), passed))


if __name__ == "__main__":
    main()
