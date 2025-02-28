import argparse
import csv
import sys

def load_gene_transcript_mapping(mapping_file):
    """Load gene_id to transcript_id mapping from the provided TSV file."""
    gene_to_transcripts = {}
    with open(mapping_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue
            gene_id, transcript_id = row[0], row[1]
            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = set()
            gene_to_transcripts[gene_id].add(transcript_id)
    return gene_to_transcripts

def filter_counts_by_gene(mapping_file, count_file, target_genes):
    """Filter the count TSV file based on transcript IDs matching the target gene(s)."""
    gene_to_transcripts = load_gene_transcript_mapping(mapping_file)

    # Collect relevant transcript IDs
    target_transcripts = set()
    for gene in target_genes:
        if gene in gene_to_transcripts:
            target_transcripts.update(gene_to_transcripts[gene])

    # Read and filter the count file
    with open(count_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(sys.stdout, delimiter='\t')

        header_written = False
        for row in reader:
            if row[0].startswith('#'):  # Column header row
                writer.writerow(row)
                header_written = True
                continue

            if row[0] in target_transcripts:  # Filtered transcript rows
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Filter transcript count file by gene_id.")
    parser.add_argument(
        "--mapping", required=True, help="Path to the gene-transcript mapping TSV file."
    )
    parser.add_argument(
        "--counts", required=True, help="Path to the transcript count TSV file."
    )
    parser.add_argument(
        "--genes", nargs='+', required=True, help="List of gene_id values to filter."
    )

    args = parser.parse_args()

    filter_counts_by_gene(args.mapping, args.counts, args.genes)

if __name__ == "__main__":
    main()