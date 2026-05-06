#!/usr/bin/env python3
"""Standalone wrapper to run IsoQuant barcode detection from the command line."""
import argparse
import sys
import os

# Add IsoQuant2 to path
ISOQUANT_DIR = os.environ.get("ISOQUANT_DIR", os.path.expanduser("~/IsoQuant2"))
sys.path.insert(0, ISOQUANT_DIR)

from isoquant_lib.modes import IsoQuantMode
from isoquant_lib.barcode_calling.detect_barcodes import process_single_thread, process_in_parallel


def parse_args():
    parser = argparse.ArgumentParser(description="Run IsoQuant barcode detection")
    parser.add_argument("--fastq", "-i", required=True, nargs="+", help="Input FASTQ/FASTA file(s)")
    parser.add_argument("--mode", required=True, help="IsoQuant mode (e.g. tenX_v3_split)")
    parser.add_argument("--barcode_whitelist", required=True, nargs="+", dest="barcodes",
                        help="Barcode whitelist file(s)")
    parser.add_argument("-o", "--output", required=True, help="Output prefix")
    parser.add_argument("--output_sequences", action="store_true", help="Output split FASTA")
    parser.add_argument("-t", "--threads", type=int, default=1)
    return parser.parse_args()


def main():
    args = parse_args()

    # Resolve mode
    try:
        iq_mode = IsoQuantMode[args.mode]
    except KeyError:
        print("Unknown mode: %s. Valid modes: %s" % (args.mode, [m.name for m in IsoQuantMode]))
        sys.exit(1)

    output_tsv = [args.output + ".barcode_calls.tsv"]
    out_fasta = [args.output + ".split_reads.fasta"] if args.output_sequences else None

    # Build args namespace expected by detect_barcodes internals
    run_args = argparse.Namespace(
        mode=iq_mode,
        barcodes=args.barcodes,
        input=args.fastq,
        output_tsv=output_tsv,
        out_fasta=out_fasta,
        molecule=None,
        threads=args.threads,
        output=args.output,   # needed for tmp_dir path in parallel mode
        tmp_dir=None,
    )

    if args.threads > 1:
        process_in_parallel(run_args)
    else:
        process_single_thread(run_args)


if __name__ == "__main__":
    main()
