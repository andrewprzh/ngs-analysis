import argparse
import csv
import json
import os
from functools import reduce

from gtfparse import read_gtf
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Generate transcripts with alternative TSS sites")
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF annotation file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix (produces <prefix>.gtf, <prefix>_no_genes.gtf, "
                             "<prefix>_novel_transcripts.tsv, <prefix>_tss.json)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--fraction", type=float, default=0.25,
                        help="Fraction of transcripts to generate alternative TSS for (default: 0.25)")
    return parser.parse_args()


def compute_tss_positions(transcript_df):
    strand_per_tid = transcript_df.groupby('transcript_id')['strand'].first()

    exon_df = transcript_df[transcript_df['feature'] == 'exon']
    exon_max_end = exon_df.groupby('transcript_id')['end'].max()
    exon_min_start = exon_df.groupby('transcript_id')['start'].min()

    all_max_end = transcript_df.groupby('transcript_id')['end'].max()
    all_min_start = transcript_df.groupby('transcript_id')['start'].min()

    max_end = exon_max_end.reindex(strand_per_tid.index).fillna(all_max_end)
    min_start = exon_min_start.reindex(strand_per_tid.index).fillna(all_min_start)

    # TSS: min(start) for +, max(end) for - (opposite of polyA)
    position = pd.Series(np.nan, index=strand_per_tid.index)
    plus = strand_per_tid == '+'
    minus = strand_per_tid == '-'
    position[plus] = min_start[plus]
    position[minus] = max_end[minus]

    result = position.reset_index()
    result.columns = ['transcript_id', 'position']
    return result


args = parse_args()

out_gtf = args.output + ".gtf"
out_gtf_no_genes = args.output + "_no_genes.gtf"
out_novel_tr = args.output + "_novel_transcripts.tsv"
out_json = args.output + "_tss.json"

os.makedirs(os.path.dirname(out_gtf) or ".", exist_ok=True)

df = read_gtf(args.gtf).to_pandas()

np.random.seed(args.seed)
transcript_ids = df['transcript_id'].unique()
transcript_ids = np.delete(transcript_ids, 0)
transcript_ids = np.random.choice(transcript_ids, size=int(len(transcript_ids) * args.fraction))
edited_transcript = df.query('transcript_id in @transcript_ids').reset_index(drop=True)

aggregated = compute_tss_positions(edited_transcript)

new_id = dict(zip(
    aggregated['transcript_id'],
    aggregated['transcript_id'] + '_' + aggregated['position'].astype(int).astype(str)
))

edited_transcript['transcript_id'] = edited_transcript['transcript_id'].map(new_id)

# Find TSS exon: first exon for + strand, last exon for - strand (opposite of polyA)
edited_transcript['indices'] = edited_transcript.index
exons = edited_transcript[edited_transcript['feature'] == 'exon'].sort_values(['transcript_id', 'start'])

terminal_exons = exons.groupby('transcript_id', group_keys=False).apply(
    lambda g: g.iloc[-1] if g.iloc[0]['strand'] == '-' else g.iloc[0],
    include_groups=False
).reset_index(drop=True)
indices = terminal_exons['indices'].values

# Vectorized length adjustment
n = len(indices)
coin_flips = np.random.rand(n)
offsets_abs = np.random.uniform(50, 300, size=n)
signs = np.where(coin_flips < 0.5, 1.0, -1.0)
offsets = (offsets_abs * signs).astype(int)

strands = edited_transcript.loc[indices, 'strand'].values
ends = edited_transcript.loc[indices, 'end'].values.copy()
starts = edited_transcript.loc[indices, 'start'].values.copy()

plus_mask = strands == '+'
minus_mask = ~plus_mask

# + strand: adjust start, guard against start >= end
p_off = offsets[plus_mask]
p_end, p_start = ends[plus_mask], starts[plus_mask]
bad = (p_start + p_off) >= p_end
edited_transcript.loc[indices[plus_mask], 'start'] = np.where(bad, p_start - p_off, p_start + p_off)

# - strand: adjust end, guard against start >= end
m_off = offsets[minus_mask]
m_end, m_start = ends[minus_mask], starts[minus_mask]
bad = m_start >= (m_end + m_off)
edited_transcript.loc[indices[minus_mask], 'end'] = np.where(bad, m_end - m_off, m_end + m_off)

df = pd.concat([df, edited_transcript.drop('indices', axis=1)]).reset_index(drop=True)

standard_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
attr_cols = [c for c in df.columns if c not in standard_cols]

for col in attr_cols:
    df[col] = '"' + df[col].astype(str) + '"'

parts = [col + ' ' + df[col] for col in attr_cols]
df['attributes'] = reduce(lambda a, b: a + '; ' + b, parts)

output_df = df.drop(columns=attr_cols, errors='ignore')
output_df.to_csv(out_gtf, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, na_rep='.')

output_df[output_df['feature'] != 'gene'].to_csv(
    out_gtf_no_genes, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, na_rep='.')

pd.Series(edited_transcript.transcript_id.unique()).to_csv(out_novel_tr, sep='\t', index=False)

aggregated2 = compute_tss_positions(df)
aggregated2 = aggregated2.dropna(subset=['position'])
new_id2 = dict(zip(aggregated2['transcript_id'], aggregated2['position'].astype(int)))
new_id2.pop('""')

with open(out_json, 'w') as f:
    json.dump(new_id2, f)
