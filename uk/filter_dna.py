#!/usr/bin/env python3
import argparse
import csv

def read_positions(tsv_files):
    """Read chromosome and position pairs from multiple TSV files."""
    positions = set()
    for path in tsv_files:
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or len(row) < 2:
                    continue
                chrom, pos = row[0], row[1]
                positions.add((chrom, pos))
    return positions

def filter_large_file(large_tsv, positions, output):
    """Filter the large TSV file keeping only matching (chrom, position) pairs."""
    with open(large_tsv, 'r') as fin, open(output, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if not row or len(row) < 2:
                continue
            chrom, pos = row[0], row[1]
            if (chrom, pos) in positions:
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description="Filter a large TSV file by positions present in one or more smaller TSV files."
    )
    parser.add_argument(
        "--input", "-i", nargs='+', required=True,
        help="Input TSV files containing chromosome and position columns."
    )
    parser.add_argument(
        "--large", "-l", required=True,
        help="Large TSV file to be filtered."
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output TSV file for filtered results."
    )

    args = parser.parse_args()

    print("Reading positions from input TSV files...")
    positions = read_positions(args.input)
    print(f"Collected {len(positions):,} unique positions.")

    print("Filtering large TSV file...")
    filter_large_file(args.large, positions, args.output)
    print(f"Filtered output written to {args.output}")

if __name__ == "__main__":
    main()
