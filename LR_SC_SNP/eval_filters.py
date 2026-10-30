#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Evaluate phase 2 filters on simulated data.

Takes the labelled calls from eval_calls.py (PREFIX.sites.tsv: TP / FP) and the output of celltype_filter.py
and / or spatial_filter.py, and reports how many true and false calls pass. For the cell-type filter on data with
cell-type truth (CT_MODE / CELL_TYPES in the truth VCF) it also reports, per carrier mode, whether the best type and
the significant types are among the true carrier types.
"""

import argparse
from collections import defaultdict


def parse_info(info):
    return dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)


def load_truth_types(fname):
    types = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        info = parse_info(v[7])
        if "CT_MODE" in info:
            types[v[2]] = (info["CT_MODE"], set(info["CELL_TYPES"].split(",")))
    return types


def load_table(fname):
    rows = {}
    header = None
    for line in open(fname):
        v = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            header = [h.lstrip("#") for h in v]
            continue
        rows[v[0]] = dict(zip(header, v))
    return rows


def alt_bin(n):
    for lo, hi in [(0, 20), (20, 50), (50, 100), (100, 10 ** 9)]:
        if lo <= n < hi:
            return "%d-%d" % (lo, hi - 1) if hi < 10 ** 9 else ">=%d" % lo


def report_pass(title, calls, filt):
    print("\n%s: pass rate by label" % title)
    print("label\tcalls\ttested\tpass\tpass_rate")
    for label in ["TP", "FP"]:
        sel = [vid for vid, c in calls.items() if c["label"] == label]
        tested = [vid for vid in sel if vid in filt and filt[vid].get("pvalue", "NA") not in ("NA", ".")]
        passed = [vid for vid in sel if vid in filt and filt[vid]["pass"] == "1"]
        print("%s\t%d\t%d\t%d\t%.3f" % (label, len(sel), len(tested), len(passed), len(passed) / len(sel) if sel else 0))
    print("TP pass rate by alt UMIs")
    by_bin = defaultdict(lambda: [0, 0])
    for vid, c in calls.items():
        if c["label"] == "TP":
            b = by_bin[alt_bin(int(c["alt_umis"]))]
            b[0] += 1
            b[1] += vid in filt and filt[vid]["pass"] == "1"
    for b in ["0-19", "20-49", "50-99", ">=100"]:
        if b in by_bin:
            n, p = by_bin[b]
            print("%s\t%d\t%d\t%.3f" % (b, n, p, p / n))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval_sites", required=True, help="PREFIX.sites.tsv from eval_calls.py")
    parser.add_argument("--truth", help="truth VCF (for CT_MODE / CELL_TYPES)")
    parser.add_argument("--celltype_sites", help="PREFIX.celltype_sites.tsv from celltype_filter.py")
    parser.add_argument("--spatial_sites", help="PREFIX.spatial_sites.tsv from spatial_filter.py")
    args = parser.parse_args()

    calls = load_table(args.eval_sites)
    if args.celltype_sites:
        ct = load_table(args.celltype_sites)
        report_pass("Cell-type filter", calls, ct)
        truth_types = load_truth_types(args.truth) if args.truth else {}
        if truth_types:
            print("\nCell-type filter on true calls by carrier mode")
            print("mode\tcalls\tpass\tpass_rate\tbest_type_correct\tsignificant_types_correct\tsignificant_types_all")
            stats = defaultdict(lambda: [0, 0, 0, 0, 0])
            for vid, c in calls.items():
                if c["label"] != "TP" or c["truth_id"] not in truth_types or vid not in ct:
                    continue
                mode, true_types = truth_types[c["truth_id"]]
                s = stats[mode]
                s[0] += 1
                if ct[vid]["pass"] != "1":
                    continue
                s[1] += 1
                s[2] += ct[vid]["best_cell_type"] in true_types
                significant = ct[vid]["significant_cell_types"].split(",")
                s[3] += sum(t in true_types for t in significant)
                s[4] += len(significant)
            for mode in sorted(stats):
                n, p, best, sig_ok, sig_all = stats[mode]
                print("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d" % (mode, n, p, p / n if n else 0, best / p if p else 0,
                                                             sig_ok / sig_all if sig_all else 0, sig_all))
    if args.spatial_sites:
        report_pass("Spatial filter", calls, load_table(args.spatial_sites))


if __name__ == "__main__":
    main()
