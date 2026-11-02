#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Phase 2 cell-type filter: is a variant enriched in a cell type?

For every variant from call_variants.py and every cell type, alt and depth UMIs are summed over the barcodes of that
type and compared with all other typed barcodes by a one-sided Fisher's exact test (alt vs ref UMIs, type vs rest;
computed as the hypergeometric upper tail). P-values are corrected by Benjamini-Hochberg over all variant x type
tests. Barcodes without a cell type are ignored. With --restrict_to only the variants also called in another call
set are tested (e.g. sites called on both sequencing platforms), and BH runs over their tests only.

Outputs:
  PREFIX.celltype_tests.tsv  one line per variant and cell type with depth in the type;
  PREFIX.celltype_sites.tsv  one line per variant: best type, its statistics, all types passing --fdr.
"""

import argparse
import gzip
from collections import defaultdict

import numpy
from scipy.stats import hypergeom


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


def load_variant_ids(sites, restrict_to=None):
    """Variant ids of a sites.tsv in file order, without repeats; with restrict_to only those also in that sites.tsv."""
    def read_ids(fname):
        return [line.rstrip("\n").split("\t")[2] for line in open(fname) if not line.startswith("#")]
    variant_ids = list(dict.fromkeys(read_ids(sites)))
    if restrict_to:
        keep = set(read_ids(restrict_to))
        variant_ids = [v for v in variant_ids if v in keep]
    return variant_ids


def add_restrict_argument(parser):
    parser.add_argument("--restrict_to", help="sites.tsv of another call set: test only the variants called there too")


def load_cell_types(fname):
    cell_types = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if len(v) >= 2:
            cell_types[v[0]] = v[1]
    return cell_types


def load_variant_type_counts(fname, cell_types, type_index, variants):
    """variant id -> (alt, depth) arrays over cell types, for the variants in the set variants."""
    n_types = len(type_index)
    counts = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        vid, barcode, alt, depth = line.rstrip("\n").split("\t")
        t = cell_types.get(barcode)
        if t is None or vid not in variants:
            continue
        if vid not in counts:
            counts[vid] = (numpy.zeros(n_types, dtype=numpy.int64), numpy.zeros(n_types, dtype=numpy.int64))
        counts[vid][0][type_index[t]] += int(alt)
        counts[vid][1][type_index[t]] += int(depth)
    return counts


def benjamini_hochberg(pvalues):
    p = numpy.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = numpy.argsort(p)
    ranked = p[order] * n / numpy.arange(1, n + 1)
    q = numpy.minimum.accumulate(ranked[::-1])[::-1]
    result = numpy.empty(n)
    result[order] = numpy.minimum(q, 1.0)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="PREFIX.sites.tsv from call_variants.py")
    parser.add_argument("--barcode_counts", required=True, help="PREFIX.barcode_counts.tsv.gz from call_variants.py")
    parser.add_argument("--cell_types", required=True, help="TSV with barcode and cell type")
    add_restrict_argument(parser)
    parser.add_argument("--fdr", type=float, default=0.05, help="BH-adjusted p-value cut-off [%(default)s]")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    cell_types = load_cell_types(args.cell_types)
    types = sorted(set(cell_types.values()))
    type_index = {t: i for i, t in enumerate(types)}
    # a site covered by two loci used to be reported twice; its barcode counts are summed
    variant_ids = load_variant_ids(args.sites, args.restrict_to)
    counts = load_variant_type_counts(args.barcode_counts, cell_types, type_index, set(variant_ids))

    # one test per variant and type with depth in the type and in the rest
    tests = []  # (variant id, type index, alt_t, depth_t, alt_rest, depth_rest)
    for vid in variant_ids:
        if vid not in counts:
            continue
        alt, depth = counts[vid]
        alt_total, depth_total = alt.sum(), depth.sum()
        for t in numpy.flatnonzero(depth):
            if depth_total - depth[t] > 0:
                tests.append((vid, t, alt[t], depth[t], alt_total - alt[t], depth_total - depth[t]))
    if tests:
        a = numpy.array([x[2] for x in tests])
        d = numpy.array([x[3] for x in tests])
        ar = numpy.array([x[4] for x in tests])
        dr = numpy.array([x[5] for x in tests])
        # P(X >= a) for X ~ Hypergeom(population = d + dr, successes = a + ar, draws = d)
        pvalues = hypergeom.sf(a - 1, d + dr, a + ar, d)
        qvalues = benjamini_hochberg(pvalues)
        # odds ratio with 0.5 pseudocounts
        odds = ((a + 0.5) * (dr - ar + 0.5)) / ((d - a + 0.5) * (ar + 0.5))
    else:
        pvalues = qvalues = odds = numpy.zeros(0)

    by_variant = defaultdict(list)
    with open(args.output + ".celltype_tests.tsv", "w") as outf:
        outf.write("#variant_id\tcell_type\talt_umis\tdepth_umis\talt_umis_rest\tdepth_umis_rest\tvaf\tvaf_rest"
                   "\todds_ratio\tpvalue\tqvalue\n")
        for i, (vid, t, a_t, d_t, a_r, d_r) in enumerate(tests):
            outf.write("%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.3g\t%.3g\n" %
                       (vid, types[t], a_t, d_t, a_r, d_r, a_t / d_t, a_r / d_r, odds[i], pvalues[i], qvalues[i]))
            by_variant[vid].append(i)

    passed = 0
    with open(args.output + ".celltype_sites.tsv", "w") as outf:
        outf.write("#variant_id\tbest_cell_type\talt_umis\tdepth_umis\tvaf\tvaf_rest\todds_ratio\tpvalue\tqvalue"
                   "\tsignificant_cell_types\tpass\n")
        for vid in variant_ids:
            idx = by_variant.get(vid)
            if not idx:
                outf.write("%s\t.\t0\t0\t.\t.\t.\t.\t.\t.\t0\n" % vid)
                continue
            best = min(idx, key=lambda i: (pvalues[i], -odds[i]))
            _, t, a_t, d_t, a_r, d_r = tests[best]
            significant = [types[tests[i][1]] for i in sorted(idx, key=lambda i: pvalues[i]) if qvalues[i] < args.fdr]
            is_pass = int(bool(significant))
            passed += is_pass
            outf.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.3g\t%.3g\t%s\t%d\n" %
                       (vid, types[t], a_t, d_t, a_t / d_t, a_r / d_r if d_r else 0, odds[best], pvalues[best],
                        qvalues[best], ",".join(significant) if significant else ".", is_pass))
    print("Variants: %d, with typed barcodes: %d, tests: %d, passing FDR %g: %d" %
          (len(variant_ids), len(by_variant), len(tests), args.fdr, passed))


if __name__ == "__main__":
    main()
