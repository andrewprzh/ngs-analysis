import pysam
import argparse


def extract_tags(bam_file, output_tsv):
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_tsv, "w") as out:
        out.write("#read_id\tCarcode\tUMI\tscore\tUMI_trusted\n")
        for read in bam:
            tags = dict(read.tags)
            CB = tags.get("CB", "*")
            if CB.endswith("-1"):
                CB = CB[:-2]
            UB = tags.get("UB", "*")

            read_id = read.query_name
            read_id = read_id[:-2]

            # Add fixed values
            col4 = 30
            col5 = True

            out.write(f"{read_id}\t{CB}\t{UB}\t{col4}\t{col5}\n")


def main():
    parser = argparse.ArgumentParser(description="Extract ReadID, CB, UB and fixed columns from BAM.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_tsv", help="Output TSV file")

    args = parser.parse_args()
    extract_tags(args.bam_file, args.output_tsv)


if __name__ == "__main__":
    main()
