#!/usr/bin/env python3
import gffutils
import pandas as pd
import argparse

def compute_transcript_lengths(db_path, output_tsv):
    db = gffutils.FeatureDB(db_path, keep_order=True)

    rows = []
    for transcript in db.features_of_type("transcript"):
        transcript_id = transcript.attributes.get("transcript_id", [""])[0]
        transcript_name = transcript.attributes.get("transcript_name", [""])[0] if "transcript_name" in transcript.attributes else ""

        # collect exons
        exons = list(db.children(transcript, featuretype="exon", order_by="start"))
        length = sum(exon.end - exon.start + 1 for exon in exons)

        rows.append((transcript_id, transcript_name, length))

    # Save TSV
    df = pd.DataFrame(rows, columns=["transcript_id", "transcript_name", "exon_length"])
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Written: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Compute exon-based transcript lengths from a gffutils database."
    )
    parser.add_argument(
        "-d", "--db",
        required=True,
        help="Path to gffutils database (.db)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )
    args = parser.parse_args()

    compute_transcript_lengths(args.db, args.output)

if __name__ == "__main__":
    main()
