import argparse
import collections
import matplotlib.pyplot as plt


def filter_sequences(input_tsv, output_tsv, output_png):
    counts = {}

    # Read TSV
    with open(input_tsv, "r") as infile:
        header = infile.readline()  # skip header
        for line in infile:
            seq, count = line.strip().split("\t")
            count = int(count)
            counts[seq] = count

    # Filter sequences
    filtered_counts = {}
    for seq, count in counts.items():
        t_fraction = seq.count("T") / len(seq)
        trailing_t = 0
        for base in reversed(seq):
            if base == "T":
                trailing_t += 1
            else:
                break

        if t_fraction > 0.5:
            continue
        if trailing_t > 10:
            continue

        filtered_counts[seq] = count

    # Write filtered TSV
    with open(output_tsv, "w") as out:
        out.write("CR\tCount\n")
        for seq, count in sorted(filtered_counts.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{seq}\t{count}\n")

    # Prepare histogram data
    freq_counts = collections.Counter(filtered_counts.values())
    x = sorted(freq_counts.keys())
    y = [freq_counts[val] for val in x]

    # Plot histogram (limit X-axis to 1000)
    plt.figure(figsize=(8, 6))
    plt.bar(x, y)
    plt.xlabel("Count (abundance of CR)")
    plt.ylabel("Number of distinct CR sequences")
    plt.title("Filtered CR Tag Abundance Histogram")
    plt.xlim(0, 1000)
    plt.yscale("log")  # optional
    plt.savefig(output_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Filter CR sequences and plot abundance histogram.")
    parser.add_argument("input_tsv", help="Input TSV file with CR counts")
    parser.add_argument("output_tsv", help="Output filtered TSV file")
    parser.add_argument("output_png", help="Output PNG histogram")

    args = parser.parse_args()
    filter_sequences(args.input_tsv, args.output_tsv, args.output_png)


if __name__ == "__main__":
    main()
