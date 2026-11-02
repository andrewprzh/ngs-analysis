#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Phase 2 region filter: is a variant enriched in a pre-defined tissue region?

For every variant from call_variants.py, alt and depth UMIs are summed over the barcodes of each region, and each
region is compared with all other regions by a one-sided Fisher's exact test (alt vs ref UMIs, region vs rest), as in
celltype_filter.py. The regions are tested twice:
  - over all barcodes with a region (cells "All");
  - within each cell type of --cell_types (or of --types): the type's barcodes in the region against the type's
    barcodes in all other regions, so that a region does not stand out only because of its cell-type composition.
    Only the variants passing --fdr over all cells are tested within types (hierarchical testing).
P-values are corrected by Benjamini-Hochberg separately over the all-cells tests (variant x region) and over the
within-type tests (passing variant x type x region). Barcodes without a region are ignored. With --restrict_to only
the variants also called in another call set are tested.

Outputs:
  PREFIX.region_tests.tsv       one line per variant, cells (All or a type) and region with depth in the region;
  PREFIX.region_sites.tsv       one line per variant: best region over all cells, its statistics, regions passing --fdr;
  PREFIX.region_type_sites.tsv  the same for the within-type tests, with regions written as type:region.
"""

import argparse

import numpy
from scipy.stats import hypergeom

from celltype_filter import add_restrict_argument, benjamini_hochberg, load_cell_types, load_variant_ids, open_text

ALL_CELLS = "All"


def load_regions(fname, column):
    regions = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if len(v) > column:
            regions[v[0]] = v[column]
    return regions


def load_variant_region_counts(fname, barcode_offsets, variant_index, block):
    """Flat alt and depth counters, block entries per variant: (all cells, region) first, then (type, region)."""
    alt = [0] * (len(variant_index) * block)
    depth = [0] * (len(variant_index) * block)
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        vid, barcode, a, d = line.rstrip("\n").split("\t")
        offsets = barcode_offsets.get(barcode)
        v = variant_index.get(vid)
        if offsets is None or v is None:
            continue
        a, d = int(a), int(d)
        all_offset, type_offset = offsets
        base = v * block
        alt[base + all_offset] += a
        depth[base + all_offset] += d
        if type_offset >= 0:
            alt[base + type_offset] += a
            depth[base + type_offset] += d
    return numpy.array(alt, dtype=numpy.int64), numpy.array(depth, dtype=numpy.int64)


def fisher_tests(alt, depth, alt_rest, depth_rest, testable):
    """One-sided Fisher's exact test of every testable (variant, cells, region) counter against the rest of the cells.
    Returns arrays with one entry per test: variant, cells, region, alt, depth, alt rest, depth rest, odds ratio,
    p-value and BH-adjusted p-value."""
    v_idx, g_idx, r_idx = numpy.nonzero(testable)
    a, d = alt[v_idx, g_idx, r_idx], depth[v_idx, g_idx, r_idx]
    ar, dr = alt_rest[v_idx, g_idx, r_idx], depth_rest[v_idx, g_idx, r_idx]
    # P(X >= a) for X ~ Hypergeom(population = d + dr, successes = a + ar, draws = d)
    pvalues = hypergeom.sf(a - 1, d + dr, a + ar, d)
    # odds ratio with 0.5 pseudocounts
    odds = ((a + 0.5) * (dr - ar + 0.5)) / ((d - a + 0.5) * (ar + 0.5))
    return v_idx, g_idx, r_idx, a, d, ar, dr, odds, pvalues, benjamini_hochberg(pvalues)


def write_sites(fname, label_column, variant_ids, labels, v_idx, stats, pvalues, qvalues, odds, fdr):
    """One line per variant: its test with the lowest p-value and all labels passing fdr. Returns the number passing."""
    a, d, ar, dr = stats
    order = numpy.lexsort((-odds, pvalues, v_idx))
    first = order[numpy.r_[True, v_idx[order][1:] != v_idx[order][:-1]]] if len(order) else order
    best = dict(zip(v_idx[first].tolist(), first.tolist()))
    significant = {}
    sig = numpy.flatnonzero(qvalues < fdr)
    for i in sig[numpy.lexsort((pvalues[sig], v_idx[sig]))].tolist():
        significant.setdefault(int(v_idx[i]), []).append(labels[i])

    with open(fname, "w") as outf:
        outf.write("#variant_id\t%s\talt_umis\tdepth_umis\tvaf\tvaf_rest\todds_ratio\tpvalue\tqvalue"
                   "\tsignificant_%ss\tpass\n" % (label_column, label_column.replace("best_", "")))
        for v, vid in enumerate(variant_ids):
            i = best.get(v)
            if i is None:
                outf.write("%s\t.\t0\t0\t.\t.\t.\t.\t.\t.\t0\n" % vid)
                continue
            passing = significant.get(v)
            outf.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.3g\t%.3g\t%s\t%d\n" %
                       (vid, labels[i], a[i], d[i], a[i] / d[i], ar[i] / dr[i], odds[i], pvalues[i], qvalues[i],
                        ",".join(passing) if passing else ".", int(bool(passing))))
    return len(significant)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="PREFIX.sites.tsv from call_variants.py")
    parser.add_argument("--barcode_counts", required=True, help="PREFIX.barcode_counts.tsv.gz from call_variants.py")
    parser.add_argument("--barcode_regions", required=True, help="TSV with barcode and region")
    parser.add_argument("--region_column", type=int, default=1,
                        help="0-based column of --barcode_regions with the region [%(default)s]")
    parser.add_argument("--cell_types", help="TSV with barcode and cell type, for the within-type tests")
    parser.add_argument("--types", help="comma-separated cell types for the within-type tests [all in --cell_types]")
    add_restrict_argument(parser)
    parser.add_argument("--fdr", type=float, default=0.05, help="BH-adjusted p-value cut-off [%(default)s]")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    barcode_regions = load_regions(args.barcode_regions, args.region_column)
    regions = sorted(set(barcode_regions.values()))
    region_index = {r: i for i, r in enumerate(regions)}
    cell_types = load_cell_types(args.cell_types) if args.cell_types else {}
    types = args.types.split(",") if args.types else sorted(set(cell_types.values()))
    type_index = {t: i for i, t in enumerate(types)}
    groups = [ALL_CELLS] + types
    n_regions = len(regions)
    # barcode -> offsets of its (all cells, region) and (type, region) counters; -1 if not of a tested type
    barcode_offsets = {}
    for barcode, region in barcode_regions.items():
        r = region_index[region]
        t = type_index.get(cell_types.get(barcode))
        barcode_offsets[barcode] = (r, -1 if t is None else (t + 1) * n_regions + r)
    del barcode_regions, cell_types

    # a site covered by two loci used to be reported twice; its barcode counts are summed
    variant_ids = load_variant_ids(args.sites, args.restrict_to)
    variant_index = {v: i for i, v in enumerate(variant_ids)}
    shape = (len(variant_ids), len(groups), n_regions)
    alt, depth = load_variant_region_counts(args.barcode_counts, barcode_offsets, variant_index,
                                            len(groups) * n_regions)
    alt, depth = alt.reshape(shape), depth.reshape(shape)
    alt_rest = alt.sum(axis=2, keepdims=True) - alt
    depth_rest = depth.sum(axis=2, keepdims=True) - depth

    # tests need depth in the region and in the rest of the cells; within types only the variants passing over all
    # cells are tested, so that the within-type BH correction is not swamped by millions of tests
    testable = (depth > 0) & (depth_rest > 0)
    over_all = (numpy.arange(len(groups)) == 0)[None, :, None]
    all_tests = fisher_tests(alt, depth, alt_rest, depth_rest, testable & over_all)
    all_passed = numpy.zeros(len(variant_ids), dtype=bool)
    all_passed[all_tests[0][all_tests[-1] < args.fdr]] = True
    type_tests = fisher_tests(alt, depth, alt_rest, depth_rest, testable & ~over_all & all_passed[:, None, None])

    with open(args.output + ".region_tests.tsv", "w") as outf:
        outf.write("#variant_id\tcells\tregion\talt_umis\tdepth_umis\talt_umis_rest\tdepth_umis_rest\tvaf\tvaf_rest"
                   "\todds_ratio\tpvalue\tqvalue\n")
        for tests in (all_tests, type_tests):
            for start in range(0, len(tests[0]), 1000000):
                for v, g, r, a_t, d_t, a_r, d_r, o, p, q in zip(*(x[start:start + 1000000].tolist() for x in tests)):
                    outf.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.3g\t%.3g\n" % (
                        variant_ids[v], groups[g], regions[r], a_t, d_t, a_r, d_r, a_t / d_t, a_r / d_r, o, p, q))

    summary = ["Variants: %d, regions: %s" % (len(variant_ids), ",".join(regions))]
    for name, tests, label_column, suffix in [
            ("all cells", all_tests, "best_region", "region_sites"),
            ("within types (variants passing over all cells)", type_tests, "best_type_region", "region_type_sites")]:
        v_idx, g_idx, r_idx, a, d, ar, dr, odds, pvalues, qvalues = tests
        if label_column == "best_region":
            labels = [regions[r] for r in r_idx.tolist()]
        else:
            labels = ["%s:%s" % (groups[g], regions[r]) for g, r in zip(g_idx.tolist(), r_idx.tolist())]
        passed = write_sites("%s.%s.tsv" % (args.output, suffix), label_column, variant_ids, labels, v_idx,
                             (a, d, ar, dr), pvalues, qvalues, odds, args.fdr)
        summary.append("%s: tests %d, variants passing FDR %g: %d" % (name, len(v_idx), args.fdr, passed))
    print("; ".join(summary))


if __name__ == "__main__":
    main()
