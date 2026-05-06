#!/usr/bin/env python3
import argparse
import csv
import gzip
import os

def open_maybe_gzip(path, mode='rt'):
    """Open a file normally or with gzip depending on its extension."""
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode, newline='')

def load_positions(position_file):
    """Load (chrom, position) pairs from a TSV file (plain or gzipped)."""
    positions = set()
    with open_maybe_gzip(position_file, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or len(row) < 2:
                continue
            chrom, pos = row[0], row[1]
            positions.add((chrom, pos))
    print(f"Loaded {len(positions):,} positions from {position_file}")
    return positions

def filter_large_file(large_tsv, positions, output):
    """Filter one large TSV file by given positions (supports gz)."""
    kept = 0
    with open_maybe_gzip(large_tsv, 'rt') as fin, open_maybe_gzip(output, 'wt') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if not row or len(row) < 2:
                continue
            chrom, pos = row[0], row[1]
            if (chrom, pos) in positions:
                writer.writerow(row)
                kept += 1
    print(f"Filtered {large_tsv}: kept {kept:,} matching rows → {output}")

def main():
    parser = argparse.ArgumentParser(
        description="Filter one or more large TSV files by positions from another TSV (supports gzipped files)."
    )
    parser.add_argument(
        "--positions", "-p", required=True,
        help="TSV or TSV.GZ file with chromosome and position pairs (from collect_positions.py)."
    )
    parser.add_argument(
        "--large", "-l", nargs='+', required=True,
        help="Large TSV or TSV.GZ files to be filtered."
    )
    parser.add_argument(
        "--output-prefix", "-o", required=True,
        help="Prefix for output TSV (or TSV.GZ) files, one per input large file."
    )
    parser.add_argument(
        "--gzip-output", action="store_true",
        help="Write output files as gzipped (.gz)."
    )

    args = parser.parse_args()
    positions = load_positions(args.positions)

    for large_file in args.large:
        base_name = os.path.basename(large_file)
        out_ext = ".tsv.gz" if args.gzip_output else ".tsv"
        output_file = f"{args.output_prefix}_{base_name}{out_ext}"
        filter_large_file(large_file, positions, output_file)

if __name__ == "__main__":
    main()
