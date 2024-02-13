############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
import gffutils
import argparse
from collections import defaultdict
import gzip
from Bio import SeqIO
from traceback import print_exc


def analyse_intron_sites(intron, chr_record, strand, seq_size=6):
    intron_left_pos = intron[0] - 1
    intron_right_pos = intron[1] - 1

    if strand not in ['+', '-']:
        return 0, 0, 0, 0

    left_upper = chr_record[intron_left_pos - seq_size:intron_left_pos]
    left_lower = chr_record[intron_left_pos + 2:intron_left_pos + seq_size + 2]
    right_upper = chr_record[intron_right_pos - seq_size - 1:intron_right_pos - 1]
    right_lower = chr_record[intron_right_pos + 1:intron_right_pos + seq_size + 1]

    # upstream and downstream here are relative to the genome
    if strand == "+":
        donor_upstream = left_upper.rfind("GT")
        donor_downstream = left_lower.find("GT")
        acc_upstream = right_upper.rfind("AG")
        acc_downstream = right_lower.find("AG")
    else:
        acc_upstream = left_upper.rfind("CT")
        acc_downstream = left_lower.find("CT")
        donor_upstream = right_upper.rfind("AC")
        donor_downstream = right_lower.find("AC")

    donor_upstream = seq_size - donor_upstream if donor_upstream != -1 else 0
    donor_downstream = 2 + donor_downstream if donor_downstream != -1 else 0
    acc_upstream = seq_size - acc_upstream if acc_upstream != -1 else 0
    acc_downstream = 2 + acc_downstream if acc_downstream != -1 else 0

    if strand == '+':
        return donor_upstream, donor_downstream, acc_upstream, acc_downstream
    else:
        return donor_downstream, donor_upstream, acc_downstream, acc_upstream


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def process_gene_db(db, chr_dict, main_gtf_fname, microexon_gtf_fname, ambigous_site_fname,
                    micro_exon_len=50, splice_site_diff=6):
    main_gtf = open(main_gtf_fname, "w")
    microexon_gtf = open(microexon_gtf_fname, "w")
    splice_site_gtf = open(ambigous_site_fname, "w")

    genes_kept = 0
    transcripts_kept = 0
    genes_micro_exon = 0
    transcripts_micro_exon = 0
    genes_splice_site = 0
    transcripts_splice_site = 0
    genes_total = 0
    transcripts_total = 0

    current_chr_id = ""
    chr_record = None
    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        if g.seqid != current_chr_id:
            current_chr_id = g.seqid
            chr_record = str(chr_dict[g.seqid].seq)
            print("Processing chromosome %s" % current_chr_id)

        genes_total += 1
        microexon_found = set()
        ambiguous_site_found = set()
        for t in db.children(g, featuretype='transcript', order_by='start'):
            transcripts_total += 1
            exon_list = []
            for e in db.children(t, featuretype='exon', order_by='start'):
                exon_len = e.end - e.start + 1
                if exon_len <= micro_exon_len:
                    microexon_found.add(t.id)
                exon_list.append((e.start, e.end))
            introns = junctions_from_blocks(sorted(exon_list))
            for intron in introns:
                if any(x != 0 for x in analyse_intron_sites(intron, chr_record, t.strand, splice_site_diff)):
                    ambiguous_site_found.add(t.id)
                    break

        gene_str = '%s\n' % g
        if len(microexon_found) > 0 or len(ambiguous_site_found) > 0:
            genes_kept += 1
            main_gtf.write(gene_str)
        if len(microexon_found) > 0:
            genes_micro_exon += 1
            microexon_gtf.write(gene_str)
        if len(ambiguous_site_found) > 0:
            genes_splice_site += 1
            splice_site_gtf.write(gene_str)

        for t in db.children(g, featuretype='transcript', order_by='start'):
            if t.id in microexon_found or t.id in ambiguous_site_found:
                transcripts_kept += 1
                dump_transcript_to_file(main_gtf, db, t)
            if t.id in microexon_found:
                transcripts_micro_exon += 1
                dump_transcript_to_file(microexon_gtf, db, t)
            if t.id in ambiguous_site_found:
                transcripts_splice_site += 1
                dump_transcript_to_file(splice_site_gtf, db, t)

    main_gtf.close()
    microexon_gtf.close()
    splice_site_gtf.close()

    print("Total genes: %d, transcripts: %d" % (genes_total, transcripts_total))
    print("Kept genes: %d, transcripts: %d" % (genes_kept, transcripts_kept))
    print("Microexon genes: %d, transcripts: %d" % (genes_micro_exon, transcripts_micro_exon))
    print("Ambiguous splice site genes: %d, transcripts: %d" % (genes_splice_site, transcripts_splice_site))

def dump_transcript_to_file(f, db, t):
    f.write('%s\n' % t)
    for e in db.children(t, featuretype='exon', order_by='start'):
        f.write('%s\n' % e)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str, required=True)
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format", type=str, required=True)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--micro_exon_len", help="microexon threshold",  type=int, default=50)
    parser.add_argument("--splice_site_diff", help="maximal length to nearest canonical dinucleotide",
                        type=int, default=6)
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    print("Loading gene db from %s" % args.genedb)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    print("Loading reference genome from %s" % args.reference)
    _, outer_ext = os.path.splitext(args.reference)
    if outer_ext.lower() in ['.gz', '.gzip']:
        with gzip.open(args.reference, "rt") as handle:
            reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        reference_record_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    print("Excluding transcripts")
    process_gene_db(gffutils_db, reference_record_dict,
                    main_gtf_fname=args.output_prefix + ".selected.gtf",
                    microexon_gtf_fname=args.output_prefix + ".microexon.gtf",
                    ambigous_site_fname=args.output_prefix + ".ambiguous_ss.gtf",
                    micro_exon_len=args.micro_exon_len,
                    splice_site_diff=args.splice_site_diff)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
