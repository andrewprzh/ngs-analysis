import pysam
import argparse
from collections import defaultdict
import random

def load_reads(bam_file):
    # Initialize a dictionary to store reads with their (gene, barcode) pairs
    reads_dict = {}

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over the reads in the BAM file
        for read in bam:
            # Check if the read has the required tags
            if read.has_tag('XS') and read.has_tag('XT') and read.has_tag('CR'):
                xs_tag = read.get_tag('XS')
                if xs_tag == 'Assigned':
                    gene = read.get_tag('XT')
                    barcode = read.get_tag('CR')
                    reads_dict[read.query_name] = (gene, barcode)

    return reads_dict

def count_unique_pairs(reads_dict):
    # Convert the dictionary items to a list and shuffle it
    reads_list = list(reads_dict.items())
    random.shuffle(reads_list)

    # Initialize a set to store unique (gene, barcode) pairs
    unique_pairs = set()

    # Calculate the intervals for 10%, 20%, 30%, etc.
    total_reads = len(reads_list)
    intervals = [int(total_reads * i / 10) for i in range(1, 11)]
    print(intervals)

    # Initialize a list to store the count of unique pairs at each interval
    unique_counts = []

    # Iterate over the shuffled reads list
    for i, (read_id, (gene, barcode)) in enumerate(reads_list):
        unique_pairs.add((gene, barcode))

        # Check if the current index is at an interval point
        if i+1 in intervals:
            unique_counts.append(len(unique_pairs))

    return unique_counts

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Count unique (gene, barcode) pairs in a BAM file.')
    parser.add_argument('bam_file', type=str, help='Path to the input BAM file')

    # Parse the arguments
    args = parser.parse_args()

    # Load reads into a dictionary
    reads_dict = load_reads(args.bam_file)

    # Count unique pairs
    unique_counts = count_unique_pairs(reads_dict)

    # Print the results
    for i, count in enumerate(unique_counts, start=1):
        print(f'Unique pairs at {i*10}% of reads: {count}')

# Example usage
main()

# Note: The main function is commented out to prevent execution in this environment.
# You can uncomment it and run the script in your local environment.

