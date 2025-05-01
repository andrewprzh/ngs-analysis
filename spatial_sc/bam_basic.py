import pysam
import argparse
import statistics
import numpy as np

def is_spliced(cigar):
    """Check if the CIGAR string contains 'N' (reference skipping)."""
    return any(op == 3 for op, _ in cigar)

def main(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    query_lengths = []
    spliced_count = 0
    total_count = 0

    for read in bam.fetch():
        if read.is_secondary or read.is_supplementary:
            continue

        query_lengths.append(read.query_length)
        if is_spliced(read.cigar):
            spliced_count += 1
        total_count += 1

    bam.close()

    if total_count == 0:
        print("No primary alignments found.")
        return

    average_query_length = statistics.mean(query_lengths)
    median_query_length = statistics.median(query_lengths)
    max_query_length = max(query_lengths)
    min_query_length = min(query_lengths)
    q1_query_length = np.percentile(query_lengths, 25)
    q3_query_length = np.percentile(query_lengths, 75)
    percentage_spliced = (spliced_count / total_count) * 100

    print(f"Avg\t{average_query_length}")
    print(f"Med\t{median_query_length}")
    print(f"Max\t{max_query_length}")
    print(f"Min\t{min_query_length}")
    print(f"1Q \t{q1_query_length}")
    print(f"3Q \t{q3_query_length}")
    print(f"Spl\t{percentage_spliced:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze BAM file for query lengths and spliced reads.")
    parser.add_argument("bam_file", help="Input BAM file")
    args = parser.parse_args()
    main(args.bam_file)
