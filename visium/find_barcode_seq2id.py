import pysam
import sys
import numpy
import argparse
import collections


def load_whitelist(file):
    print("Loading %s" % file)
    with open(file) as f:
        return [line.strip() for line in f if line.strip()]


def analyze_barcodes(bam_file, part1_file, part2_file, out_prefix):
    part1_list = load_whitelist(part1_file)
    part2_list = load_whitelist(part2_file)

    print("Create mapping dict: concatenated sequence -> (part1, part2)")
    concat_map = {}
    for p1 in part1_list:
        for p2 in part2_list:
            concat_map[p1 + p2] = (p1, p2)

    detected_pairs = collections.Counter()

    # Dicts to collect counts
    part1_dict = {seq: collections.Counter() for seq in part1_list}
    part2_dict = {seq: collections.Counter() for seq in part2_list}

    print("Processing %s" % bam_file)
    c = 0
    reads_with_bc = 0
    reads_with_valid_bc = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            c += 1
            if c % 1000000 == 0:
                sys.stdout.write("Processed %d reads\r" % c)
            tags = dict(read.tags)
            CR = tags.get("CR")  # barcode sequence
            CB = tags.get("CB")  # barcode ID
            if not CR or not CB:
                continue

            reads_with_bc += 1
            if CR not in concat_map:
                continue  # skip if not valid

            detected_pairs[CR] += 1

            part1_seq, part2_seq = concat_map[CR]

            try:
                _, _, y, x = CB[:-2].split("_")
            except ValueError:
                continue

            reads_with_valid_bc += 1
            part1_dict[part1_seq][x] += 1
            part2_dict[part2_seq][y] += 1

    print("Done. Total reads processed %d, of them with barcodes: %d, of them valid: %d" % (c, reads_with_bc, reads_with_valid_bc))

    # Write outputs for part1
    print("Outputting files to %s" % out_prefix)
    with open(f"{out_prefix}_part1_full.tsv", "w") as out1, open(f"{out_prefix}_part1_short.tsv", "w") as out2:
        out1.write("Sequence\tTopCoords\n")
        out2.write("Sequence\tTopCoord\n")
        for seq, counter in part1_dict.items():
            if not counter:
                continue
            total = sum(counter.values())
            top_items = counter.most_common(3)
            top_str = ";".join([f"{coord}:{count}:{count/total:.2%}" for coord, count in top_items])
            out1.write(f"{seq}\t{top_str}\n")
            best_coord, _ = counter.most_common(1)[0]
            out2.write(f"{seq}\t{best_coord}\n")

    # Write outputs for part2
    with open(f"{out_prefix}_part2_full.tsv", "w") as out1, open(f"{out_prefix}_part2_short.tsv", "w") as out2:
        out1.write("Sequence\tTopCoords\n")
        out2.write("Sequence\tTopCoord\n")
        for seq, counter in part2_dict.items():
            if not counter:
                continue
            total = sum(counter.values())
            top_items = counter.most_common(3)
            top_str = ";".join([f"{coord}:{count}:{count/total:.2%}" for coord, count in top_items])
            out1.write(f"{seq}\t{top_str}\n")
            best_coord, _ = counter.most_common(1)[0]
            out2.write(f"{seq}\t{best_coord}\n")

    print("Barcode pairs observed %d / %d" % (len(detected_pairs), len(concat_map)))
    bc_counts = list(detected_pairs.values())
    print("Barcode pair counts. Min: %d, Max: %d, Mean: %d, Median: %d" % (min(bc_counts), max(bc_counts), int(numpy.mean(bc_counts)), int(numpy.median(bc_counts))))


def main():
    parser = argparse.ArgumentParser(description="Analyze BAM barcodes and map CR to CD coordinates.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("part1_whitelist", help="Whitelist file for part1 sequences")
    parser.add_argument("part2_whitelist", help="Whitelist file for part2 sequences")
    parser.add_argument("out_prefix", help="Prefix for output files")

    args = parser.parse_args()
    analyze_barcodes(args.bam_file, args.part1_whitelist, args.part2_whitelist, args.out_prefix)


if __name__ == "__main__":
    main()
