#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import Counter


FASTA_EXTS = (".fasta", ".fa", ".fna")
FASTQ_EXTS = (".fastq", ".fq")


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
    raise ValueError(f"Cannot infer format from filename: {path}")


def default_output(path):
    name = path[:-3] if path.endswith(".gz") else path
    root, _ = os.path.splitext(name)
    return root + ".tsv"


def iter_read_ids(path, fmt):
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
        description="Count reads per transcript from a simulated FASTA/FASTQ "
                    "(transcript id = read id split by '_', first element). "
                    "Supports plain and gzip-compressed input."
    )
    parser.add_argument("input", help="Input FASTA/FASTQ file (optionally .gz)")
    parser.add_argument("-o", "--output", help="Output TSV (default: <input>.tsv)")
    args = parser.parse_args()

    fmt = detect_format(args.input)
    output = args.output or default_output(args.input)

    counts = Counter()
    for read_id in iter_read_ids(args.input, fmt):
        counts[read_id.split("_", 1)[0]] += 1

    with open(output, "w") as out:
        out.write("#id\tcount\n")
        for tid in sorted(counts):
            out.write(f"{tid}\t{counts[tid]}\n")


if __name__ == "__main__":
    main()
