import argparse
import gzip


def parse_gtf(gtf_file, output_tsv):
    if gtf_file.endswith(".gz"):
        infile = gzip.open(gtf_file, 'rt')
    else:
        infile = open(gtf_file, 'r')
    with open(output_tsv, 'w') as outfile:
        outfile.write("gene_id\ttranscript_id\n")  # Write header
        
        for line in infile:
            if line.startswith("#"):  # Skip comments
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "transcript":  # Ensure it's a transcript entry
                continue
            
            attributes = fields[8]
            gene_id = transcript_id = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split(' ')[1].strip('"')
                elif attr.startswith("transcript_id"):
                    transcript_id = attr.split(' ')[1].strip('"')
            
            if gene_id and transcript_id:
                outfile.write(f"{gene_id}\t{transcript_id}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene_id and transcript_id from a GTF file")
    parser.add_argument("--gtf_file", help="Input GTF file", required=True)
    parser.add_argument("--output_tsv", help="Output TSV file", required=True)
    
    args = parser.parse_args()
    parse_gtf(args.gtf_file, args.output_tsv)
