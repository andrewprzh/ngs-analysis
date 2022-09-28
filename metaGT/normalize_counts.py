import os
import sys
import argparse


def load_lengths(inf, col=3):
    eff_len = {}
    for l in open(inf):
        if l.startswith("Name"):
            continue

        v = l.strip().split()
        eff_len[v[0]] = float(v[col-1]) / 1000.0
    return eff_len


def convert_counts(inf, eff_len, col):
    rpk = {}
    for l in open(inf):
        if l.startswith("#") or l.startswith("Gene"):
            continue
        v = l.strip().split()
        rpk[v[0]] = float(v[col-1]) / eff_len[v[0]]

    scale_factor = 1000000.0 / sum(rpk.values())
    print(sum(rpk.values()), scale_factor)
    tpm_dict = {}
    for k in sorted(rpk.keys()):
        tpm_dict[k] = rpk[k] * scale_factor
    return tpm_dict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--len_file', type=str, help='file with gene lengths', required=True)
    parser.add_argument('--len_col', type=int, default=3, help='gene length column')
    parser.add_argument('--counts', '-c', type=str, help='output expression table', required=True)
    parser.add_argument('--counts_col', type=int, default=5, help='counts column in output expression table')
    parser.add_argument('--output', '-o', type=str, help='output file')
    return parser.parse_args()


def main():
    args = parse_args()
    eff_len = load_lengths(args.len_file, col=args.len_col)
    tpm_dict = convert_counts(args.counts, eff_len, args.counts_col)
    if not args.output:
        fname, ext = os.path.splitext(args.counts)
        args.output = fname + ".normalized.tsv"

    with open(args.output, "w") as outf:
        for tid in tpm_dict:
            outf.write("%s\t%.6f\n" % (tid, tpm_dict[tid]))


if __name__ == '__main__':
    main()
