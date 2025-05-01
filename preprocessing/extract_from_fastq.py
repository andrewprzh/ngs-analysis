import argparse
from Bio import SeqIO
import gzip


def filter_fastq(input_fastq, read_ids_file, output_fastq):
    # Read the read IDs from the file into a set for fast lookup
    with open(read_ids_file, 'r') as f:
        read_ids = set(line.strip() for line in f)

    # Open the output FASTQ file for writing
    with open(output_fastq, 'w') as output_handle:
        # Read the input FASTQ file and filter records
        if input_fastq.endswith('.gz'):
            handle = gzip.open(input_fastq, "rt")
        else:
            handle = open(input_fastq, "r")
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in read_ids:
                SeqIO.write(record, output_handle, "fastq")

def main():
    parser = argparse.ArgumentParser(description='Filter a FASTQ file based on a list of read IDs.')
    parser.add_argument('--input', '-i', type=str, help='Path to the input FASTQ file', required=True)
    parser.add_argument('--read_ids', '-r', type=str, help='Path to the file containing read IDs', required=True)
    parser.add_argument('--output','-o', type=str, help='Path to the output FASTQ file', required=True)

    args = parser.parse_args()

    filter_fastq(args.input, args.read_ids, args.output)

if __name__ == '__main__':
    main()
