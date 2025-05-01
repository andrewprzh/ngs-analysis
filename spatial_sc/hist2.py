import argparse
import pandas as pd

def read_tsv(file_path):
    """Read a TSV file into a pandas DataFrame."""
    return pd.read_csv(file_path, sep='\t')

def filter_non_zero_counts(lengths_df, counts_df):
    """Filter transcripts with non-zero counts."""
    # Drop rows with missing transcript_id in the counts DataFrame
    counts_df = counts_df.dropna(subset=['transcript_id'])

    # Merge the DataFrames using a left join
    merged_df = pd.merge(lengths_df, counts_df, on='transcript_id', how='left')

    # Fill any missing count values with zero
    merged_df['count'] = merged_df['count'].fillna(0)

    # Filter for transcripts with non-zero counts
    non_zero_counts_df = merged_df[merged_df['count'] > 0]
    return non_zero_counts_df

def generate_histogram(lengths, output_file):
    """Generate the length distribution histogram data and save to a TSV file."""
    bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
            1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500,
            6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]

    histogram_data = pd.cut(lengths, bins=bins, right=False).value_counts().sort_index()

    # Create a DataFrame for the histogram data
    histogram_df = pd.DataFrame({
        'bin_start': histogram_data.index.left,
        'bin_end': histogram_data.index.right,
        'count': histogram_data.values
    })

    # Save the histogram data to a TSV file
    histogram_df.to_csv(output_file, sep='\t', index=False)

def main(lengths_file, counts_file, output_file):
    lengths_df = read_tsv(lengths_file)
    counts_df = read_tsv(counts_file)
    counts_df.rename(index={0:"transcript_id", 1: "count"}, inplace=True)

    non_zero_counts_df = filter_non_zero_counts(lengths_df, counts_df)

    if non_zero_counts_df.empty:
        print("No transcripts with non-zero counts found.")
        return

    lengths = non_zero_counts_df['length']
    generate_histogram(lengths, output_file)
    print(f"Histogram data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate length distribution histogram data for transcripts with non-zero counts.")
    parser.add_argument("--lengths", required=True, help="TSV file with transcript lengths (transcript_id, length)")
    parser.add_argument("--counts", required=True, help="TSV file with transcript counts (transcript_id, count)")
    parser.add_argument("--output", required=True, help="Output TSV file to save the histogram data")
    args = parser.parse_args()
    main(args.lengths, args.counts, args.output)
