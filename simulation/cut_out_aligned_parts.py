import sys
import pysam
from traceback import print_exc
from collections import defaultdict
import argparse
from Bio import SeqIO, Seq, SeqRecord

BIN_COUNT = 100
HIST_STEP = 1.0 / BIN_COUNT


def process(transcript_dict):
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
    new_records = []
    transcript_counts = defaultdict(int)
    for alignment in bamfile:
        if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
            continue

        new_records.append(SeqRecord.SeqRecord(seq=Seq.Seq(transcript_dict[alignment.reference_name][alignment.reference_start, alignment.reference_end + 1]),
                                               id="%s_%d_%d_%d" % (alignment.reference_name, transcript_counts[alignment.reference_name], alignment.reference_start, alignment.reference_end),
                                               description=""))
        transcript_counts[alignment.reference_name] += 1
    return new_records


def parse_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)

    required_group = parser.add_argument_group('required parameters')
    required_group.add_argument("--bam", help="BAM file", type=str, required=True)
    required_group.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    required_group.add_argument("--transcripts", "-t", help="transcriptome FASTA", type=str, required=True)
    args = parser.parse_args()

    if args.bam is None:
        parser.print_help()
        exit(-1)
    return args


def main():
    args = parse_args()
    transcript_dict = SeqIO.to_dict(SeqIO.parse(args.transcripts, "fasta"))
    fasta_records = process(transcript_dict)
    count_dict = defaultdict(int)
    for r in fasta_records:
        count_dict[r.id] = 1
    SeqIO.write(fasta_records, args.output + ".fasta", "fasta")

    scale_factor = sum(count_dict.values()) / 1000000.0
    with open(sys.argv[2], "w") as outf:
        for tid in sorted(count_dict.keys()):
            outf.write("%s\t%.2f\t%.6f\n" % (tid, count_dict[tid], count_dict[tid] / scale_factor))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)