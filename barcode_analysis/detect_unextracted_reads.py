import argparse
import csv

def compare_tsv_files(file1, file2):
    # Read the first TSV file into a dictionary
    barcode_dict = {}
    with open(file1, 'r') as f1:
        reader = csv.DictReader(f1, delimiter='\t')
        for row in reader:
            read_id = row['original_read_id']
            ground_truth = row['ground_truth_barcode']
            detected = row['detected_barcode']
            if ground_truth == detected:
                barcode_dict[read_id] = detected

    # Read the second TSV file and compare
    with open(file2, 'r') as f2:
        reader = csv.reader(f2, delimiter='\t')
        for row in reader:
            read_id = row[0]
            detected = row[1]
            if detected == '*' and read_id in barcode_dict:
                print(read_id)

def main():
    parser = argparse.ArgumentParser(description='Compare two TSV files to find matching read IDs.')
    parser.add_argument('--blaze_tsv', type=str, help='Path to the first TSV file', required=True)
    parser.add_argument('--extracted', type=str, help='Path to the second TSV file', required=True)

    args = parser.parse_args()

    compare_tsv_files(args.blaze_tsv, args.extracted)

if __name__ == '__main__':
    main()
