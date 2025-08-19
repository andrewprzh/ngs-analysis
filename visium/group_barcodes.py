import pysam
import argparse
import collections
import statistics


def analyze_bam(bam_file, output_tsv):
    cb_to_cr = collections.defaultdict(collections.Counter)
    ub_counter = collections.Counter()
    ur_ub_diff = 0
    total_ur_ub = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            tags = dict(read.tags)
            CR = tags.get("CR")
            CB = tags.get("CB")
            UR = tags.get("UR")
            UB = tags.get("UB")

            # Part 1: collect CR per CB
            if CB and CR:
                cb_to_cr[CB[:-2]][CR] += 1

            # Part 2: collect UB stats
            if UB:
                ub_counter[UB] += 1
                if UR:
                    total_ur_ub += 1
                    if UR != UB:
                        ur_ub_diff += 1

    # Part 1 output: sort by top CR counts descending
    results = []
    for cb, cr_counter in cb_to_cr.items():
        if not cr_counter:
            continue
        top_cr, top_count = cr_counter.most_common(1)[0]
        total = sum(cr_counter.values())
        percentage = (top_count / total) * 100 if total > 0 else 0
        best_cr = '*'
        best_cr_count = 0
        for cr in cr_counter.keys():
            if len(cr) == 30 and cr_counter[cr] > best_cr_count:
                best_cr_count = cr_counter[cr]
                best_cr = cr
        best_perc = (best_cr_count / total) * 100 if total > 0 else 0
        results.append((cb, top_cr, top_count, percentage, best_cr, best_cr_count, best_perc))

    results.sort(key=lambda x: x[2], reverse=True)

    with open(output_tsv, "w") as out:
        out.write("CB\tTop_CR\tTop_CR_Count\tPercentage\tTop_CR_30\tTop_CR_30_Count\tPercentage\n")
        for cb, top_cr, top_count, percentage, best_cr, best_cr_count, best_perc in results:
            out.write(f"{cb}\t{top_cr}\t{top_count}\t{percentage:.2f}\t{best_cr}\t{best_cr_count}\t{best_perc:.2f}\n")

    # Part 2 output: print stats
    ub_counts = list(ub_counter.values())
    if ub_counts:
        max_count = max(ub_counts)
        median_count = statistics.median(ub_counts)
        mean_count = statistics.mean(ub_counts)
        print(f"Distinct UB tags: {len(ub_counter)}")
        print(f"Max UB count: {max_count}")
        print(f"Median UB count: {median_count}")
        print(f"Mean UB count: {mean_count:.2f}")
        print(f"UR != UB in {ur_ub_diff} of {total_ur_ub} cases")
    else:
        print("No UB tags found.")


def main():
    parser = argparse.ArgumentParser(description="Analyze BAM tags: CB->CR mapping and UB stats.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_tsv", help="Output TSV file for CB->CR stats")

    args = parser.parse_args()
    analyze_bam(args.bam_file, args.output_tsv)


if __name__ == "__main__":
    main()
