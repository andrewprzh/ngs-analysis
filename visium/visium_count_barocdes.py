import pysam
import argparse
import collections
import matplotlib.pyplot as plt


def count_cr_tags(bam_file, output_tsv, output_png):
    cr_counts = collections.Counter()

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            tags = dict(read.tags)
            CR = tags.get("CR", None)
            if CR:
                cr_counts[CR] += 1

    # Write counts to TSV, sorted by count descending
    with open(output_tsv, "w") as out:
        out.write("CR\tCount\n")
        for cr, count in cr_counts.most_common():
            out.write(f"{cr}\t{count}\n")

    # Prepare histogram data
    freq_counts = collections.Counter(cr_counts.values())
    x = sorted(freq_counts.keys())
    y = [freq_counts[val] for val in x]

    # Plot histogram
    plt.figure(figsize=(8, 6))
    plt.bar(x, y)
    plt.xlabel("Count (abundance of CR)")
    plt.ylabel("Number of distinct CR sequences")
    plt.title("CR Tag Abundance Histogram")
    plt.yscale("log")  # optional, in case of skewed distribution
    plt.savefig(output_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Count CR:Z tags from BAM and plot abundance histogram.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_tsv", help="Output TSV file with CR counts")
    parser.add_argument("output_png", help="Output PNG file for histogram")

    args = parser.parse_args()
    count_cr_tags(args.bam_file, args.output_tsv, args.output_png)


if __name__ == "__main__":
    main()
