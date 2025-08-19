import pysam
import argparse
import gzip
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq


def bam_to_fastq(bam_file, output_file):
    # Decide on compression
    if output_file.endswith(".gz") or output_file.endswith(".gzip"):
        handle = gzip.open(output_file, "wt")
    else:
        handle = open(output_file, "w")

    buffer = []
    processed_ids = set()

    with pysam.AlignmentFile(bam_file, "r", check_sq=False) as bam:
        for a in bam:
            seq = a.get_forward_sequence()
            if not seq:
                continue

            read_id = a.query_name
            if read_id in processed_ids:
                continue
            processed_ids.add(read_id)

            qual = a.get_forward_qualities()
            buffer.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(seq),
                    id=read_id,
                    description="",
                    letter_annotations={"phred_quality": qual},
                )
            )

            if len(buffer) > 100000:
                SeqIO.write(buffer, handle, "fastq")
                buffer = []

    SeqIO.write(buffer, handle, "fastq")
    handle.close()


def main():
    parser = argparse.ArgumentParser(description="Convert BAM to FASTQ, with optional gzip compression.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_file", help="Output FASTQ file (.gz/.gzip for compression)")

    args = parser.parse_args()
    bam_to_fastq(args.bam_file, args.output_file)


if __name__ == "__main__":
    main()
