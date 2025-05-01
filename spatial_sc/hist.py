import argparse
import csv
from collections import defaultdict

def read_tsv(file_path):
    """Read a TSV file into a dictionary."""
    data = {}
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            transcript_id = row['transcript_id']
            data[transcript_id] = row
    return data

def filter_non_zero_counts(lengths, counts):
    """Filter transcripts with non-zero counts."""
    non_zero_counts = []
    for transcript_id, length_info in lengths.items():
        count_info = counts.get(transcript_id, {'count': 0})
        count = int(float(count_info['count']))
        if count > 0:
            non_zero_counts.append(int(length_info['length']))
    return non_zero_counts

def generate_histogram(lengths, output_file):
    """Generate the length distribution histogram data and save to a TSV file."""
    bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
            1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500,
            6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]

    histogram_data = defaultdict(int)
    for length in lengths:
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                histogram_data[(bins[i], bins[i + 1])] += 1
                break

    # Save the histogram data to a TSV file
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['bin_start', 'bin_end', 'count'])
        for (bin_start, bin_end), count in sorted(histogram_data.items()):
            writer.writerow([bin_start, bin_end, count])

def main(lengths_file, counts_file, output_file):
    lengths = read_tsv(lengths_file)
    counts = read_tsv(counts_file)

    non_zero_counts = filter_non_zero_counts(lengths, counts)

    if not non_zero_counts:
        print("No transcripts with non-zero counts found.")
        return

    generate_histogram(non_zero_counts, output_file)
    print(f"Histogram data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate length distribution histogram data for transcripts with non-zero counts.")
    parser.add_argument("--lengths", required=True, help="TSV file with transcript lengths (transcript_id, length)")
    parser.add_argument("--counts", required=True, help="TSV file with transcript counts (transcript_id, count)")
    parser.add_argument("--output", required=True, help="Output TSV file to save the histogram data")
    args = parser.parse_args()
    main(args.lengths, args.counts, args.output)
