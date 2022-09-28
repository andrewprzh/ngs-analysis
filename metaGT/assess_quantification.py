
import os
import subprocess
import argparse
import pathlib
import scipy

from collections import Counter

import numpy as np


def load_counts(inf, count_col=2, id_col=1):
    print("Loading TPM values from " + inf)
    count_dict = {}
    for l in open(inf):
        if l.startswith("#") or l.startswith("__") or l.startswith("target_id") or l.startswith("Gene") or l.startswith("Name")  or l.startswith("feature_id"):
            continue
        v = l.strip().split()
        count = float(v[count_col - 1])
        if count > 0:
            count_dict[v[id_col-1]] = count
    print("Counts loaded, checksum %.2f" % sum(count_dict.values()))
    return count_dict


def normalize_counts(count_dict):
    scale_factor = 1000000 / sum(count_dict.values())
    return {k:scale_factor * count_dict[k] for k in count_dict}


def count_deviation(joint_dict):
    print("Counting deviation histogram")
    deviation_values = []
    false_detected = 0

    for t_id in joint_dict.keys():
        exp_pair = joint_dict[t_id]
        if exp_pair[0] == 0:
            if exp_pair[1] > 0:
                false_detected += 1
            continue
        deviation_values.append(100 * exp_pair[1] / exp_pair[0])

    print("Total %d, false %d, missed %d" % (len(deviation_values), false_detected, deviation_values.count(0.0)))
    bins = [10 * i for i in range(21)]
    bins.append(10000)
    dev_vals, bins = np.histogram(deviation_values, bins)
    mid_bins = bins[:-1]
    return zip(mid_bins, dev_vals)


def count_stats(joint_dict, output, prefix="", header=""):
    ref_tpms = []
    real_tpms = []
    n_isoforms = 0
    counts_reported = 0
    full_matches = 0
    close_matches_10 = 0
    close_matches_20 = 0
    false_detected = 0
    not_detected = 0
    for t_id in joint_dict.keys():
        ref_expr = joint_dict[t_id][0]
        real_expr = joint_dict[t_id][1]
        if ref_expr == 0 and real_expr == 0:
            continue

        ref_tpms.append(ref_expr)
        real_tpms.append(real_expr)
        if ref_expr > 0:
            n_isoforms += 1
        if real_expr > 0:
            counts_reported += 1

        if real_expr == ref_expr:
            full_matches += 1
        if real_expr <= 1.1 * ref_expr and real_expr >= 0.9 * ref_expr and real_expr > 0:
            close_matches_10 += 1
        if real_expr <= 1.2 * ref_expr and real_expr >= 0.8 * ref_expr and real_expr > 0:
            close_matches_20 += 1
        if real_expr > 0 and ref_expr == 0:
            false_detected += 1
        if real_expr == 0 and ref_expr > 0:
            not_detected += 1

    outf = open(os.path.join(output, prefix + "stats.tsv"), 'w')
    outf.write(header + "\n")
    outf.write("Total reference isoforms %d\n" % n_isoforms)
    outf.write("Total non-zero counts %d\n" % counts_reported)
    print(scipy.stats.spearmanr(real_tpms, ref_tpms))
    outf.write('Correlation\t%.3f\n' % round(scipy.stats.spearmanr(real_tpms, ref_tpms).correlation, 3))
    outf.write('Full matches\t%d\t%.3f\n' % (full_matches, round(full_matches / n_isoforms, 3)))
    outf.write('Close matches (10)\t%d\t%.3f\n' % (close_matches_10, round(close_matches_10 / n_isoforms, 3)))
    outf.write('Close matches (20)\t%d\t%.3f\n' % (close_matches_20, round(close_matches_20 / n_isoforms, 3)))
    outf.write('Not detected\t%d\t%.3f\n' % (not_detected, round(not_detected / n_isoforms, 3)))
    outf.write('False detections\t%d\t%.4f\n' % (false_detected, round(false_detected / n_isoforms, 4)))
    outf.close()


def compare_transcript_counts(ref_tpm_dict, tpm_dict, output, prefix="", header=""):
    print("Filling true values")
    joint_dict = {}
    for tid in tpm_dict.keys():
        if tid in ref_tpm_dict:
            joint_dict[tid] = (ref_tpm_dict[tid], tpm_dict[tid])
        else:
            joint_dict[tid] = (0, tpm_dict[tid])
    for tid in ref_tpm_dict:
        if tid not in joint_dict:
            joint_dict[tid] = (ref_tpm_dict[tid], 0)

    print("Saving TPM values")
    with open(os.path.join(output, prefix + "tpm.values.tsv"), 'w') as out_tpms:
        for t_id in joint_dict.keys():
            if joint_dict[t_id] == (0, 0):
                continue
            out_tpms.write("%s\t%.4f\t%.4f\n" % (t_id, joint_dict[t_id][0], joint_dict[t_id][1]))

    with open(os.path.join(output, prefix + "deviation.tsv"), 'w') as out_dev:
        for hist_pairs in count_deviation(joint_dict):
            out_dev.write("%d\t%d\n" % (hist_pairs[0], hist_pairs[1]))

    count_stats(joint_dict, output, header=header, prefix=prefix)
    print("Done")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_expr', '-r', type=str, help='reference expression table, TPM', required=True)
    parser.add_argument('--ref_col', type=int, default=4, help='counts column in reference expression table')
    parser.add_argument('--counts', '-c', type=str, help='output expression table to assess, TPM')
    parser.add_argument('--counts_col', type=int, default=5, help='counts column in output expression table')

    parser.add_argument('--output', '-o', type=str, help='output folder', default="quantification_assessment")
    return parser.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    ref_count_dict = load_counts(args.ref_expr, args.ref_col)
    res_count_dict = load_counts(args.counts, args.counts_col)
    ref_tpm_dict = normalize_counts(ref_count_dict)
    print("Normalized reference, checksum %.2f" % sum(ref_tpm_dict.values()))
    res_tpm_dict = normalize_counts(res_count_dict)
    print("Normalized results, checksum %.2f" % sum(res_tpm_dict.values()))

    compare_transcript_counts(ref_count_dict, res_count_dict, args.output, prefix="RAW.")
    compare_transcript_counts(ref_tpm_dict, res_tpm_dict, args.output, prefix="NORM.")


if __name__ == '__main__':
    main()
