import pysam
import argparse


def load_reads_to_dict(bam_file, n_reads):
    """
    Load first N reads from BAM into a dictionary.
    """
    read_dict = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i, read in enumerate(bam):
            if i >= n_reads:
                break
            tags = dict(read.tags)
            read_dict[read.query_name.split('_')[0]] = (
                read.query_sequence,
                tags.get("CR", None),
                tags.get("CB", None),
                tags.get("1R", None),
                tags.get("UR", None),
                tags.get("UB", None),
            )
    return read_dict


def compare_bams(first_bam, second_bam, n_reads):
    # Load dictionary from first BAM
    read_dict = load_reads_to_dict(first_bam, n_reads)

    # Iterate over second BAM and check for matches
    with pysam.AlignmentFile(second_bam, "rb", check_sq=False) as bam:
        reads_found = 0
        for read in bam:
            if reads_found > 1000:
                break
            if read.query_name in read_dict:
                reads_found += 1
                seq1, CR, CB, tag1R, UR, UB = read_dict[read.query_name]
                print(
                    f"{read.query_name}\t{read.query_sequence}\t{seq1}\tCR:{CR}\tCB:{CB}\t1R:{tag1R}\tUR:{UR}\tUB:{UB}"
                )


def main():
    parser = argparse.ArgumentParser(description="Compare reads between two BAM files.")
    parser.add_argument("first_bam", help="First BAM file (used to build dictionary)")
    parser.add_argument("second_bam", help="Second BAM file (to compare against dictionary)")
    parser.add_argument("-n", "--num_reads", type=int, default=1000, help="Number of reads to load from first BAM")

    args = parser.parse_args()
    compare_bams(args.first_bam, args.second_bam, args.num_reads)


if __name__ == "__main__":
    main()
