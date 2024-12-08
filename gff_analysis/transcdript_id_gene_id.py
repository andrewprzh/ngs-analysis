import sys
import csv


def parse_gtf(gtf_file, output_tsv):
    """
    Extracts transcript IDs and gene IDs from a GTF file and writes them to a TSV file.

    Args:
        gtf_file (str): Path to the input GTF file.
        output_tsv (str): Path to the output TSV file.
    """
    with open(gtf_file, 'r') as gtf, open(output_tsv, 'w', newline='') as tsv:
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(["Transcript ID", "Gene ID"])  # Header row
        
        for line in gtf:
            if line.startswith("#"):
                continue  # Skip header/comment lines
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "transcript":
                continue  # Skip lines without sufficient fields or not transcripts
            
            attributes = fields[8]
            gene_id = None
            transcript_id = None
            
            # Parse attributes field
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split(' ')[1].strip('"')
                elif attr.startswith("transcript_id"):
                    transcript_id = attr.split(' ')[1].strip('"')
            
            if gene_id and transcript_id:
                writer.writerow([transcript_id, gene_id])

# Example usage:
# Replace 'input.gtf' with your GTF file path and 'output.tsv' with your desired output path
parse_gtf(sys.argv[1], sys.argv[2])
