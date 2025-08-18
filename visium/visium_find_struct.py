import argparse
import collections


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(base, "N") for base in reversed(seq))


def process_tsv(input_tsv, output_tsv):
    preceding_counts = collections.Counter()
    following_counts = collections.Counter()

    with open(input_tsv, "r") as infile, open(output_tsv, "w") as outfile:
        header = [
            "read_id",
            "read_seq",
            "UR_pos",
            "CR_pos",
            "preceding_seq",
            "middle_seq",
            "following_seq",
        ]
        outfile.write("\t".join(header) + "\n")

        for line in infile:
            if line.strip() == "" or line.startswith("read_id"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            read_seq = fields[1]
            CR = fields[3].replace("CR:", "") if fields[3] != "CR:None" else None
            UR = fields[6].replace("UR:", "") if fields[6] != "UR:None" else None

            if not CR or not UR:
                continue

            # Search CR in sequence (direct)
            CR_pos = read_seq.find(CR)
            reversed_flag = False
            if CR_pos == -1:
                # Try reverse complement
                read_seq = reverse_complement(read_seq)
                CR_pos = read_seq.find(CR)
                reversed_flag = True

            if CR_pos == -1:
                continue  # CR not found at all

            UR_pos = read_seq.find(UR)
            if UR_pos == -1:
                continue  # UR not found

            # Extract sequences
            preceding_seq = read_seq[max(0, UR_pos - 20):UR_pos]
            middle_seq = read_seq[UR_pos + len(UR):CR_pos] if UR_pos < CR_pos else ""
            if not middle_seq: middle_seq = "*"
            following_seq = read_seq[CR_pos + len(CR):CR_pos + len(CR) + 20]

            # Count context sequences
            preceding_counts[preceding_seq] += 1
            following_counts[following_seq] += 1

            # Write result line
            outfile.write(
                f"{read_id}\t{read_seq}\t{UR_pos}\t{CR_pos}\t{preceding_seq}\t{middle_seq}\t{following_seq}\n"
            )

    # Report most common sequences
    print("Most common preceding sequences:")
    for seq, count in preceding_counts.most_common(5):
        print(f"{seq}\t{count}")

    print("Most common following sequences:")
    for seq, count in following_counts.most_common(5):
        print(f"{seq}\t{count}")


def main():
    parser = argparse.ArgumentParser(description="Process TSV to extract UR/CR positions and context sequences.")
    parser.add_argument("input_tsv", help="Input TSV file from BAM comparison script")
    parser.add_argument("output_tsv", help="Output TSV file with positions and context sequences")

    args = parser.parse_args()
    process_tsv(args.input_tsv, args.output_tsv)


if __name__ == "__main__":
    main()
