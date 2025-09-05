#!/usr/bin/env python3
import gffutils
import pandas as pd
import argparse

def extract_transcript_gene_mapping(db_path, output_tsv):
    db = gffutils.FeatureDB(db_path, keep_order=True)
    rows = []

    for transcript in db.features_of_type("transcript"):
        transcript_id = transcript.attributes.get("transcript_id", [""])[0]
        transcript_name = transcript.attributes.get("transcript_name", [""])[0] if "transcript_name" in transcript.attributes else ""

        # find parent gene
        gene = list(db.parents(transcript, featuretype="gene"))
        if gene:
            gene_id = gene[0].attributes.get("gene_id", [""])[0]
            gene_name = gene[0].attributes.get("gene_name", [""])[0]
        else:
            gene_id = ""
            gene_name = ""

        rows.append((transcript_id, transcript_name, gene_id, gene_name))

    df = pd.DataFrame(rows, columns=["transcript_id", "transcript_name", "gene_id", "gene_name"])
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Written: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Extract transcript-to-gene mapping from a gffutils database."
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

    extract_transcript_gene_mapping(args.db, args.output)

if __name__ == "__main__":
    main()
