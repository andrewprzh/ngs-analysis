#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gffutils
from Bio import SeqIO
import gzip


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output TSV file", default="polya.tsv")
    parser.add_argument("--genedb", "-d", type=str, help="gene db", required=True)
    parser.add_argument("--reference", "-r", type=str, help="reference genome", required=True)
    parser.add_argument("--min_polya_len", type=int, help="minimal stretch to consider", default=9)
    parser.add_argument("--purely_introns", type=bool, help="considers purely intronic stretches if True, may be insised exons otherwise", default=False)

    args = parser.parse_args()
    return args


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def get_avg_transcript_len(genedb, gene):
    transcript_lengths = []
    for t in genedb.children(gene, featuretype=('transcript', 'mRNA'), order_by='start'):
        for e in genedb.children(t, order_by='start'):
            if e.featuretype == 'exon':
                transcript_lengths.append(e.end - e.start + 1)
    return int(sum(transcript_lengths) / len(transcript_lengths))


def find_polya_stretches(args, chr_record, start, end, strand):
    genic_seq = str(chr_record.seq[start:end+1])
    nucl = 'A' if strand == '=' else 'T'
    sample = nucl * args.min_polya_len
    resulting_stretches = []
    pos = 0
    pos = genic_seq.find(sample, pos)
    while pos != -1 and pos < len(genic_seq):
        i = pos + 1
        while i < len(genic_seq) and genic_seq[i] == nucl:
            i += 1
        resulting_stretches.append((start + pos, start + i - 1))
        pos = i
        pos = genic_seq.find(sample, pos)

    return resulting_stretches


def find_intronic_stretches(args, genedb, gene, polya_stretches):
    exons = set()
    introns = set()
    for t in genedb.children(gene, featuretype=('transcript', 'mRNA'), order_by='start'):
        exon_list = []
        for e in genedb.children(t, order_by='start'):
            if e.featuretype == 'exon':
                exon_list.append((e.start, e.end))
        exons.update(exon_list)
        for i in range(1, len(exon_list)):
            introns.add((exon_list[i-1][1] + 1, exon_list[i][0] - 1))

    if args.purely_introns:
        resulting_list = list(filter(lambda x: not any(overlaps(x, exon) for exon in exons), polya_stretches))
    else:
        resulting_list = list(filter(lambda x: any(overlaps(x, intron) for intron in introns), polya_stretches))
    return resulting_list


def process_genome(args, genedb, reference, outstream):
    for g in genedb.features_of_type('gene', order_by=('seqid', 'start')):
        chr_record = reference[g.seqid]
        stretches = find_polya_stretches(args, chr_record, g.start, g.end, g.strand)
        intronic_stretches = find_intronic_stretches(args, genedb, g, stretches)
        for s in intronic_stretches:
            outstream.write("%s\t%s\t%s\t%d\t%d\n" % (g.seqid, g.id, g.strand, s[0], s[1]))


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print("Loading gene database from " + args.genedb)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    print("Loading reference genome from " + args.reference)
    _, outer_ext = os.path.splitext(args.reference)
    if outer_ext.lower() in ['.gz', '.gzip']:
        with gzip.open(args.reference, "rt") as handle:
            reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        reference_record_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    with open(args.output, 'w') as outstream:
        process_genome(args, gffutils_db, reference_record_dict, outstream)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
