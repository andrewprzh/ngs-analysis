import argparse
import csv
from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a TSV file to calculate barcode assignment statistics.")
    parser.add_argument("tsv_file", type=str, help="Path to the input TSV file.")
    return parser.parse_args()


def extract_ground_truth_barcodes(read_id):
    parts = read_id.split('_')
    ground_truth_barcodes = set([parts[i] for i in range(3, len(parts), 3) if i < len(parts) and len(parts[i]) == 25])
    return ground_truth_barcodes


def calculate_statistics(tsv_file):
    read_barcodes = defaultdict(list)

    # First pass: collect all barcodes for each read
    with open(tsv_file, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            read_id = "_".join(row[0].split('_')[:-3])
            if read_id.startswith("#"): continue
            detected_barcode = row[1]
            read_barcodes[read_id].append(detected_barcode)

    total_reads = len(read_barcodes)
    no_barcodes_assigned = 0
    all_correct = 0
    some_correct = 0
    total_assignments = 0
    correct_assignments = 0
    incorrect_assignments = 0
    total_true_barcodes = 0
    excessive_assignments = 0

    # Second pass: calculate statistics
    for read_id, detected_barcodes in read_barcodes.items():
        ground_truth_barcodes = extract_ground_truth_barcodes(read_id)
        total_true_barcodes += len(ground_truth_barcodes)

        if all(bc == '*' for bc in detected_barcodes):
            no_barcodes_assigned += 1
            continue
        
        correct_in_read = sum(db in detected_barcodes for db in ground_truth_barcodes)
        incorrect_in_read = sum(db not in ground_truth_barcodes for db in detected_barcodes if db != '*')
        #if incorrect_in_read > 0:
        #    print(read_id, ground_truth_barcodes, detected_barcodes, incorrect_in_read)
        excessive_assignments += max(0, len(detected_barcodes) - len(ground_truth_barcodes))
        total_assignments += min(len(detected_barcodes), len(ground_truth_barcodes))
        correct_assignments += correct_in_read
        incorrect_assignments += incorrect_in_read
        
        if correct_in_read == len(ground_truth_barcodes):
            all_correct += 1
        elif correct_in_read > 0:
            some_correct += 1

    precision = correct_assignments / total_assignments if total_assignments > 0 else 0
    recall = correct_assignments / total_true_barcodes if total_true_barcodes > 0 else 0

    print(f"Total reads: {total_reads}")
    print(f"Total true barcodes: {total_true_barcodes}")
    print(f"Reads with no barcodes assigned: {no_barcodes_assigned} ({no_barcodes_assigned / total_reads * 100:.2f}%)")
    print(f"Reads with all barcodes assigned correctly: {all_correct} ({all_correct / total_reads * 100:.2f}%)")
    print(f"Reads with some barcodes assigned correctly: {some_correct} ({some_correct / total_reads * 100:.2f}%)")
    print(f"Correct assignments: {correct_assignments} ({correct_assignments / total_assignments * 100:.2f}%)")
    print(f"Incorrect assignments: {incorrect_assignments} ({incorrect_assignments / total_assignments * 100:.2f}%)")
    print(f"Excessive assignments: {excessive_assignments} ({excessive_assignments / total_assignments * 100:.2f}%)")
    print(f"Overall assignment precision: {precision:.2f}")
    print(f"Overall assignment recall: {recall:.2f}")

def main():
    args = parse_arguments()
    calculate_statistics(args.tsv_file)

if __name__ == "__main__":
    main()

