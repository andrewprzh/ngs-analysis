#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Compare phase 1 calls (call_variants.py) with the simulation truth (simulate_barcoded_mutations.py).

Truth indels are left-aligned against the reference, as the caller does. Sites are matched by chrom, pos, ref, alt.
Recall is stratified by template-level truth: alt templates (ALT_MOL) and barcodes with alt templates;
about a third of templates are never sequenced, so these are upper bounds of what reads show.
Barcode-level accuracy for true calls: called alt barcodes vs barcodes with alt templates.

Outputs PREFIX.sites.tsv (calls with labels), PREFIX.truth.tsv (truth with call status) and a summary on stdout.
"""

import argparse
import gzip
from collections import defaultdict

import pysam

ALT_BINS = [(0, 1), (1, 10), (10, 20), (20, 50), (50, 100), (100, 10 ** 9)]
BARCODE_BINS = [(0, 1), (1, 10), (10, 20), (20, 50), (50, 100), (100, 10 ** 9)]
VAF_BINS = [(0, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 1.01)]


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


def parse_info(info):
    return dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)


def left_align(fasta, chrom, pos, ref, alt):
    """Left-align a VCF indel with a shared anchor base (1-based pos)."""
    if len(ref) == len(alt) or ref[0] != alt[0]:
        return pos, ref, alt
    while pos > 1:
        prev = fasta.fetch(chrom, pos - 2, pos - 1).upper()
        if ref[-1] != alt[-1]:
            break
        ref, alt = prev + ref[:-1], prev + alt[:-1]
        pos -= 1
    return pos, ref, alt


def load_truth(fname, fasta):
    truth = {}
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        chrom, pos, vid, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
        info = parse_info(info)
        pos, ref, alt = left_align(fasta, chrom, int(pos), ref, alt)
        truth[(chrom, pos, ref, alt)] = {"id": vid, "type": info["TYPE"], "cf": float(info["CF"]),
                                         "alt_mol": int(info["ALT_MOL"]), "cov_mol": int(info["COV_MOL"]),
                                         "ct_mode": info.get("CT_MODE", "."), "cell_types": info.get("CELL_TYPES", ".")}
    return truth


def load_truth_barcodes(fname, truth_ids):
    """variant id -> set of barcodes with alt templates, and the number of such barcodes."""
    alt_barcodes = defaultdict(set)
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        vid, barcode, is_carrier, alt, cov = line.rstrip("\n").split("\t")
        if int(alt) > 0 and vid in truth_ids:
            alt_barcodes[vid].add(barcode)
    return alt_barcodes


def load_calls(fname):
    calls = {}
    header = None
    for line in open(fname):
        v = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            header = [h.lstrip("#") for h in v]
            continue
        rec = dict(zip(header, v))
        calls[(rec["chrom"], int(rec["pos"]), rec["ref"], rec["alt"])] = rec
    return calls


def load_called_barcodes(fname, variant_ids):
    called = defaultdict(set)
    for line in open_text(fname):
        if line.startswith("#"):
            continue
        vid, barcode, alt, depth = line.rstrip("\n").split("\t")
        if vid in variant_ids and int(alt) > 0:
            called[vid].add(barcode)
    return called


def range_label(lo, hi, pct=False):
    if pct:
        return "%g-%g%%" % (lo * 100, min(hi, 1) * 100)
    if hi >= 10 ** 9:
        return ">=%d" % lo
    return "%d" % lo if hi - lo == 1 else "%d-%d" % (lo, hi - 1)


def find_bin(value, bins, pct=False):
    for lo, hi in bins:
        if lo <= value < hi:
            return range_label(lo, hi, pct)
    return "NA"


def print_table(title, rows, key, bins, pct=False):
    print("\n" + title)
    print("bin\ttruth\tcalled\trecall")
    for lo, hi in bins:
        label = range_label(lo, hi, pct)
        sel = [r for r in rows if r[key] == label]
        if sel:
            called = sum(r["called"] for r in sel)
            print("%s\t%d\t%d\t%.3f" % (label, len(sel), called, called / len(sel)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth", required=True, help="truth VCF")
    parser.add_argument("--truth_cell_counts", help="truth variant_cell_counts.tsv(.gz), for barcode-level accuracy")
    parser.add_argument("--sites", required=True, help="PREFIX.sites.tsv from call_variants.py")
    parser.add_argument("--barcode_counts", help="PREFIX.barcode_counts.tsv.gz from call_variants.py")
    parser.add_argument("--reference", required=True)
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.reference)
    truth = load_truth(args.truth, fasta)
    calls = load_calls(args.sites)
    truth_alt_barcodes = {}
    if args.truth_cell_counts:
        truth_alt_barcodes = load_truth_barcodes(args.truth_cell_counts, {t["id"] for t in truth.values()})

    labels = defaultdict(lambda: defaultdict(int))
    tp_ids = set()
    with open(args.output + ".sites.tsv", "w") as outf:
        outf.write("#variant_id\ttype\tlabel\ttruth_id\talt_umis\tdepth_umis\tvaf\talt_barcodes\thomopolymer_len"
                   "\tediting_like\terror_pvalue\tunreliable\n")
        for key, rec in calls.items():
            t = truth.get(key)
            label = "TP" if t else "FP"
            if t:
                tp_ids.add(rec["variant_id"])
            labels[rec["type"]][label] += 1
            outf.write("\t".join([rec["variant_id"], rec["type"], label, t["id"] if t else ".", rec["alt_umis"],
                                  rec["depth_umis"], rec["vaf"], rec["alt_barcodes"], rec["homopolymer_len"],
                                  rec["editing_like"], rec["error_pvalue"], rec.get("unreliable", ".")]) + "\n")

    rows = []
    with open(args.output + ".truth.tsv", "w") as outf:
        outf.write("#truth_id\tchrom\tpos\tref\talt\ttype\tcarrier_fraction\talt_templates\tcovering_templates"
                   "\talt_template_barcodes\tct_mode\tcell_types\tcalled\tcalled_alt_umis\tcalled_alt_barcodes\n")
        for key, t in truth.items():
            rec = calls.get(key)
            n_bc = len(truth_alt_barcodes.get(t["id"], ()))
            row = {"type": t["type"], "called": rec is not None, "ct_mode": t["ct_mode"],
                   "alt_bin": find_bin(t["alt_mol"], ALT_BINS), "bc_bin": find_bin(n_bc, BARCODE_BINS),
                   "vaf_bin": find_bin(t["alt_mol"] / t["cov_mol"] if t["cov_mol"] else 0, VAF_BINS, True)}
            rows.append(row)
            outf.write("\t".join(map(str, [t["id"], key[0], key[1], key[2], key[3], t["type"], t["cf"], t["alt_mol"],
                                           t["cov_mol"], n_bc, t["ct_mode"], t["cell_types"], int(rec is not None),
                                           rec["alt_umis"] if rec else 0, rec["alt_barcodes"] if rec else 0])) + "\n")

    print("type\tcalls\tTP\tFP\tprecision\ttruth\trecall")
    for vtype in ["SNV", "INS", "DEL"]:
        tp, fp = labels[vtype]["TP"], labels[vtype]["FP"]
        n_truth = sum(1 for r in rows if r["type"] == vtype)
        recall = sum(1 for r in rows if r["type"] == vtype and r["called"])
        print("%s\t%d\t%d\t%d\t%.3f\t%d\t%.3f" % (vtype, tp + fp, tp, fp, tp / (tp + fp) if tp + fp else 0,
                                                  n_truth, recall / n_truth if n_truth else 0))
    print_table("Recall by alt templates (template-level truth)", rows, "alt_bin", ALT_BINS)
    if truth_alt_barcodes:
        print_table("Recall by barcodes with alt templates", rows, "bc_bin", BARCODE_BINS)
    print_table("Recall by template VAF (alt / covering templates)", rows, "vaf_bin", VAF_BINS, True)
    modes = sorted(set(r["ct_mode"] for r in rows) - {"."})
    if modes:
        print("\nRecall by cell type mode")
        print("mode\ttruth\tcalled\trecall")
        for mode in modes:
            sel = [r for r in rows if r["ct_mode"] == mode]
            called = sum(r["called"] for r in sel)
            print("%s\t%d\t%d\t%.3f" % (mode, len(sel), called, called / len(sel)))

    if args.barcode_counts and truth_alt_barcodes:
        called_barcodes = load_called_barcodes(args.barcode_counts, tp_ids)
        id_to_truth = {calls[k]["variant_id"]: t["id"] for k, t in truth.items() if k in calls}
        tp_bc, n_called, n_truth = 0, 0, 0
        for vid, truth_id in id_to_truth.items():
            called = called_barcodes.get(vid, set())
            true_bc = truth_alt_barcodes.get(truth_id, set())
            tp_bc += len(called & true_bc)
            n_called += len(called)
            n_truth += len(true_bc)
        print("\nBarcodes of true calls: called alt barcodes %d, among them with alt templates %d (precision %.3f); "
              "barcodes with alt templates %d (recall %.3f, at most the sequenced share)" %
              (n_called, tp_bc, tp_bc / n_called if n_called else 0, n_truth, tp_bc / n_truth if n_truth else 0))


if __name__ == "__main__":
    main()
