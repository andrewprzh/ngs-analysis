import argparse
import re


def process_fastq(input_fastq, output_tsv):
    # Regular expression to extract the required parts from the sequence name
    pattern = re.compile(r'@(?P<detected_barcode>[ACGT]{16})_(?P<rest>.*)')

    with open(input_fastq, 'r') as fastq_file, open(output_tsv, 'w') as tsv_file:
        # Write the header to the TSV file
        tsv_file.write("original_read_id\tground_truth_barcode\tdetected_barcode\n")

        for line in fastq_file:
            if line.startswith('@'):
                match = pattern.match(line)
                if match:
                    detected_barcode = match.group('detected_barcode')
                    rest = match.group('rest')
                    # Extract the part after the hash tag
                    original_read_id = rest.split('#', 1)[1] if '#' in rest else rest
                    original_read_id = original_read_id.split()[0][:-2]
                    # Extract the ground truth barcode from the remaining part
                    parts = original_read_id.split('_')
                    if len(parts) > 3:
                        ground_truth_barcode = parts[3]
                    else:
                        ground_truth_barcode = '*'
                    tsv_file.write(f"{original_read_id}\t{ground_truth_barcode}\t{detected_barcode}\n")


def main():
    parser = argparse.ArgumentParser(description='Process a FASTQ file to extract barcodes and IDs.')
    parser.add_argument('--input_fastq', type=str, required=True, help='Path to the input FASTQ file')
    parser.add_argument('--output_tsv', '-o', required=True, type=str, help='Path to the output TSV file')

    args = parser.parse_args()

    process_fastq(args.input_fastq, args.output_tsv)

if __name__ == '__main__':
    main()
