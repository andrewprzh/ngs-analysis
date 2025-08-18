import pysam
import argparse
from collections import defaultdict

def count_unique_pairs(bam_file):
    # Initialize a dictionary to store unique (gene, barcode) pairs

    unique_pairs = set()
    total_reads = pysam.AlignmentFile(bam_file, "rb").count()

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:

        # Calculate the intervals for 10%, 20%, 30%, etc.
        intervals = [int(total_reads * i / 10) for i in range(1, 11)]
        current_interval = 0

        # Initialize a list to store the count of unique pairs at each interval
        unique_counts = [0] * len(intervals)

        # Iterate over the reads in the BAM file
        for i, read in enumerate(bam):
            # Check if the read has the required tags
            if read.has_tag('XS') and read.has_tag('CR'):
                xs_tag = read.get_tag('XS')
                if xs_tag == 'Assigned' and read.has_tag('XT'):
                    gene = read.get_tag('XT')
                    barcode = read.get_tag('CR')
                    unique_pairs.add((gene, barcode))

            if i + 1 == intervals[current_interval]:
                unique_counts[current_interval] = len(unique_pairs)
                current_interval += 1

    return unique_counts

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Count unique (gene, barcode) pairs in a BAM file.')
    parser.add_argument('bam_file', type=str, help='Path to the input BAM file')

    # Parse the arguments
    args = parser.parse_args()

    # Count unique pairs
    unique_counts = count_unique_pairs(args.bam_file)

    # Print the results
    for i, count in enumerate(unique_counts, start=1):
        print(f'Unique pairs at {i*10}% of reads: {count}')

# Example usage
main()

# Note: The main function is commented out to prevent execution in this environment.
# You can uncomment it and run the script in your local environment.
