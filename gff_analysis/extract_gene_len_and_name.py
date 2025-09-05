#!/usr/bin/env python3
import gffutils
import pandas as pd
import argparse

def compute_gene_lengths_longest_transcript(db_path, output_tsv):
    db = gffutils.FeatureDB(db_path, keep_order=True)
    rows = []

    for gene in db.features_of_type("gene"):
        gene_id = gene.attributes.get("gene_id", [""])[0]
        gene_name = gene.attributes.get("gene_name", [""])[0]

        max_length = 0
        for transcript in db.children(gene, featuretype="transcript", order_by="start"):
            exons = list(db.children(transcript, featuretype="exon", order_by="start"))
            transcript_length = sum(exon.end - exon.start + 1 for exon in exons)
            if transcript_length > max_length:
                max_length = transcript_length

        rows.append((gene_id, gene_name, max_length))

    df = pd.DataFrame(rows, columns=["gene_id", "gene_name", "longest_transcript_exon_length"])
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Written: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Compute gene lengths as exonic length of the longest transcript from a gffutils database."
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

    compute_gene_lengths_longest_transcript(args.db, args.output)

if __name__ == "__main__":
    main()
