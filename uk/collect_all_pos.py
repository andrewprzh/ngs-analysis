#!/usr/bin/env python3
import argparse
import csv
import gzip

def open_maybe_gzip(path, mode='rt'):
    """Open a file normally or with gzip depending on its extension."""
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode, newline='')


def collect_positions(tsv_files, output):
    """Collect unique (chrom, position) pairs from multiple TSV files."""
    positions = set()
    for path in tsv_files:
        with open_maybe_gzip(path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or len(row) < 2:
                    continue
                chrom, pos = row[0], row[1]
                positions.add((chrom, pos))

    with open(output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for chrom, pos in sorted(positions):
            writer.writerow([chrom, pos])

    print(f"Collected {len(positions):,} unique positions and wrote them to {output}")

def main():
    parser = argparse.ArgumentParser(
        description="Collect unique chromosome-position pairs from TSV files."
    )
    parser.add_argument(
        "--input", "-i", nargs='+', required=True,
        help="Input TSV files containing chromosome and position columns."
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output TSV file with unique (chrom, position) pairs."
    )

    args = parser.parse_args()
    collect_positions(args.input, args.output)

if __name__ == "__main__":
    main()
