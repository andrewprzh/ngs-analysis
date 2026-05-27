#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import Counter


FASTA_EXTS = (".fasta", ".fa", ".fna")
FASTQ_EXTS = (".fastq", ".fq")
BAM_EXTS = (".bam", ".sam", ".cram")


def open_seq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_format(path):
    name = path[:-3] if path.endswith(".gz") else path
    low = name.lower()
    if low.endswith(FASTA_EXTS):
        return "fasta"
    if low.endswith(FASTQ_EXTS):
        return "fastq"
    if low.endswith(BAM_EXTS):
        return "bam"
    raise ValueError(f"Cannot infer format from filename: {path}")


def default_output(path):
    name = path[:-3] if path.endswith(".gz") else path
    root, _ = os.path.splitext(name)
    return root + ".tsv"


def iter_read_ids(path, fmt):
    if fmt == "bam":
        try:
            import pysam
        except ImportError as e:
            raise SystemExit("pysam is required for BAM/SAM/CRAM input: pip install pysam") from e
        mode = {"bam": "rb", "sam": "r", "cram": "rc"}[path.lower().rsplit(".", 1)[-1]]
        with pysam.AlignmentFile(path, mode, check_sq=False) as af:
            for rec in af:
                if rec.is_secondary or rec.is_supplementary:
                    continue
                yield rec.query_name
        return

    with open_seq(path) as f:
        if fmt == "fasta":
            for line in f:
                if line.startswith(">"):
                    yield line[1:].split()[0]
        else:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline()
                plus = f.readline()
                qual = f.readline()
                if not qual:
                    break
                if not header.startswith("@"):
                    raise ValueError(f"Malformed FASTQ record at: {header!r}")
                yield header[1:].split()[0]


def main():
    parser = argparse.ArgumentParser(
        description="Count reads per transcript from a simulated FASTA/FASTQ/BAM "
                    "(transcript id = read id split by '_', first element). "
                    "Supports plain and gzip-compressed FASTA/FASTQ, plus SAM/BAM/CRAM."
    )
    parser.add_argument("input", help="Input FASTA/FASTQ (optionally .gz) or SAM/BAM/CRAM file")
    parser.add_argument("-o", "--output", help="Output TSV (default: <input>.tsv)")
    parser.add_argument('--polya', action='store_true', default=False, help='some transcript may have undersdcore in their id')
    args = parser.parse_args()

    fmt = detect_format(args.input)
    output = args.output or default_output(args.input)

    counts = Counter()
    for read_id in iter_read_ids(args.input, fmt):
        v = read_id.split("_", 5)
        if args.polya and len(v) > 3 and v[3] == 'aligned':
            transcript_id = v[0] + '_' + v[1]
        else:
            transcript_id = v[0]
        counts[transcript_id] += 1

    with open(output, "w") as out:
        out.write("#id\tcount\n")
        for tid in sorted(counts):
            out.write(f"{tid}\t{counts[tid]}\n")


if __name__ == "__main__":
    main()
