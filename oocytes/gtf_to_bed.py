#!/usr/bin/env python3
"""
Convert gene records from GTF format to BED format.

GTF format (1-based coordinates):
seqname source feature start end score strand frame attributes

BED format (0-based start, 1-based end):
chromosome start end name score strand
"""

import sys
import re


def parse_gtf_attributes(attr_string):
    """Parse GTF attributes field to extract gene information."""
    gene_id = None
    gene_name = None

    # Match gene_id and gene_name from attributes
    gene_id_match = re.search(r'gene_id "([^"]+)"', attr_string)
    gene_name_match = re.search(r'gene_name "([^"]+)"', attr_string)

    if gene_id_match:
        gene_id = gene_id_match.group(1)
    if gene_name_match:
        gene_name = gene_name_match.group(1)

    # Prefer gene_name, fall back to gene_id
    return gene_name if gene_name else gene_id


def gtf_to_bed(input_file, output_file=None):
    """Convert GTF gene records to BED format."""
    output = open(output_file, 'w') if output_file else sys.stdout

    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Skip comments
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                seqname, source, feature, start, end, score, strand, frame, attributes = fields

                # Only process gene features
                if feature != 'gene':
                    continue

                # Convert coordinates: GTF is 1-based, BED is 0-based start
                bed_start = int(start) - 1
                bed_end = int(end)

                # Parse gene name from attributes
                gene_name = parse_gtf_attributes(attributes)
                if not gene_name:
                    gene_name = '.'

                # Use score from GTF, or '.' if not available
                bed_score = score if score != '.' else '0'

                # Write BED format: chr start end name score strand
                output.write(f"{seqname}\t{bed_start}\t{bed_end}\t{gene_name}\t{bed_score}\t{strand}\n")

    finally:
        if output_file:
            output.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python gtf_to_bed.py <input.gtf> [output.bed]", file=sys.stderr)
        print("If output file is not specified, writes to stdout", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    gtf_to_bed(input_file, output_file)
