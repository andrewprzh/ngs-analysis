#!/usr/bin/env python3
import argparse
import csv

def load_mapping(mapping_file):
    """Load mapping from 2um to 8um spot IDs into a dictionary."""
    mapping = {}
    with open(mapping_file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            two_um, eight_um, sixteen_um = row
            mapping[two_um[:-2]] = eight_um[:-2]
    return mapping

def replace_spot_ids(input_file, mapping, output_file):
    """Replace 2um spot IDs in input file with 8um spot IDs and write output."""
    with open(input_file, newline='') as infile, open(output_file, "w", newline='') as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        header = next(reader)
        writer.writerow(header)

        for row in reader:
            two_um_id = row[1]
            if two_um_id in mapping:
                row[1] = mapping[two_um_id]
            else:
                row[1] = "*"  # or keep unchanged: row[1] = two_um_id
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description="Replace 2um spot IDs with 8um spot IDs in a TSV file."
    )
    parser.add_argument("input_tsv", help="Input TSV file with read_id and 2um spot IDs")
    parser.add_argument("mapping_tsv", help="Mapping TSV (2um, 8um, 16um)")
    parser.add_argument("output_tsv", help="Output TSV file with 8um spot IDs")
    args = parser.parse_args()

    mapping = load_mapping(args.mapping_tsv)
    replace_spot_ids(args.input_tsv, mapping, args.output_tsv)

if __name__ == "__main__":
    main()
