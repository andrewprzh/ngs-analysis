############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
import gffutils
import random
import argparse
from collections import defaultdict
from Bio import SeqIO
from traceback import print_exc


def process_gene_db(db, main_gtf_fname, expressed_gtf_fname, expressed_kept_fname, excluded_gtf_fname,
                    expressed_transcripts,
                    probability=0.2,
                    min_expressed_transcripts=1,
                    min_expression=1,
                    min_excluded_expression=2,
                    remove_unspliced=True):
    main_gtf = open(main_gtf_fname, "w")
    excl_gtf = open(excluded_gtf_fname, "w")
    expr_gtf = open(expressed_gtf_fname, "w")
    expr_kept_gtf = open(expressed_kept_fname, "w")

    genes_kept = 0
    transcripts_kept = 0
    genes_reduced = 0
    transcripts_dropped = 0
    genes_expressed = 0
    transcripts_expressed = 0
    transcripts_kept_expressed = 0

    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        exon_count = {}
        expressed = set()
        all_gene_transcripts = 0
        for t in db.children(g, featuretype='transcript', order_by='start'):
            all_gene_transcripts += 1
            t_id = t.id.split('.')[0]
            if len(expressed_transcripts) == 0 or expressed_transcripts[t_id] >= min_expression:
                expressed.add(t_id)
            exon_count[t_id] = 0
            for e in db.children(t, featuretype='exon', order_by='start'):
                exon_count[t_id] += 1

        if not remove_unspliced:
            to_remove = set()
            for t_id in exon_count.keys():
                if exon_count[t_id] == 1:
                    to_remove.add(t_id)
            for t_id in to_remove:
                del exon_count[t_id]

        ignored_transcripts = set()
        kept_expressed = set()
        if len(expressed) >= min_expressed_transcripts:
            for t_id in sorted(expressed):
                if expressed_transcripts[t_id] < min_excluded_expression or \
                        len(ignored_transcripts) == all_gene_transcripts - 1:
                    kept_expressed.add(t_id)
                    continue
                r = random.random()
                if r < probability:
                    ignored_transcripts.add(t_id)
                else:
                    kept_expressed.add(t_id)

        gene_str = '%s\n' % g
        assert len(ignored_transcripts) < all_gene_transcripts
        main_gtf.write(gene_str)
        genes_kept += 1
        if expressed:
            expr_gtf.write(gene_str)
            genes_expressed += 1
        if ignored_transcripts:
            excl_gtf.write(gene_str)
            genes_reduced += 1
        if kept_expressed:
            expr_kept_gtf.write(gene_str)

        for t in db.children(g, featuretype='transcript', order_by='start'):
            t_id = t.id.split('.')[0]
            if t_id in ignored_transcripts:
                dump_transcript_to_file(excl_gtf, db, t)
                transcripts_dropped += 1
            else:
                dump_transcript_to_file(main_gtf, db, t)
                transcripts_kept += 1

            if t_id in expressed:
                transcripts_expressed += 1
                dump_transcript_to_file(expr_gtf, db, t)
            if t_id in kept_expressed:
                transcripts_kept_expressed += 1
                dump_transcript_to_file(expr_kept_gtf, db, t)

    main_gtf.close()
    excl_gtf.close()
    expr_gtf.close()
    expr_kept_gtf.close()

    print("Expressed genes: %d, transcripts: %d" % (genes_expressed, transcripts_expressed))
    print("Kept genes: %d, transcripts: %d" % (genes_kept, transcripts_kept))
    print("Reduced genes: %d, dropped transcripts: %d" % (genes_reduced, transcripts_dropped))
    print("Kept expressed trnanscripts: %d" % transcripts_kept_expressed)

def dump_transcript_to_file(f, db, t):
    f.write('%s\n' % t)
    for e in db.children(t, featuretype='exon', order_by='start'):
        f.write('%s\n' % e)


def read_expressed_isoforms(fname):
    expressed_transcripts = defaultdict(int)
    base, ext = os.path.splitext(fname)
    if ext in {'.fa', '.fasta', '.fna'}:
        for record in SeqIO.parse(fname, "fasta"):
            tokens = record.name.split('_')
            if len(tokens) < 2:
                print("Invalid fasta entry, isoform id was not found: " + record.name)
                continue
            t_id = tokens[1].split('.')[0]
            assert  t_id.startswith("E")
            expressed_transcripts[t_id] += 1
    else:
        for l in open(fname):
            t = l.strip().split()
            t_id = t[0].split('.')[0]
            assert t_id.startswith("E")
            expressed_transcripts[t_id] = int(t[1])

    return expressed_transcripts


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str, required=True)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--seed", help="randomization seed [11]", type=int, default=11)
    parser.add_argument("--min_expression", help="minimal number reads in transcript to be considered as expressed",
                        type=int, default=1)
    parser.add_argument("--min_excluded_expression", help="minimal number reads in transcript to be possibly excluded",
                        type=int, default=4)
    parser.add_argument("--min_transcripts_expressed", help="minimal number of transcripts expressed in a gene to "
                                                            "consider them for dropping", type=int, default=1)

    parser.add_argument("--fraction", help="fraction of transcripts to be removed",
                        type=float, default=0.2)
    parser.add_argument("--expressed", help="simulated reads or count table, "
                                            "will consider the entire annotation if not set", type=str)
    args = parser.parse_args()

    if args.fraction < 0 or args.fraction > 1:
        print("Invalid fraction value")
        sys.exit(-1)
    return args


def main():
    args = parse_args()
    random.seed(args.seed)

    expressed_transcripts = {}
    if args.expressed:
        print("Reading expression from " + args.expressed)
        expressed_transcripts = read_expressed_isoforms(args.expressed)
        expressed_fname=args.output_prefix + ".expression.tsv"
        print("Saving expression table to " + expressed_fname)
        f = open(expressed_fname, 'w')
        for k in expressed_transcripts:
            f.write("%s\t%d\n" % (k, expressed_transcripts[k]))
        f.close()

    print("Loading gene db")
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    print("Excluding transcripts")
    process_gene_db(gffutils_db,
                    main_gtf_fname=args.output_prefix + ".reduced.gtf",
                    expressed_gtf_fname=args.output_prefix + ".expressed.gtf",
                    excluded_gtf_fname=args.output_prefix + ".excluded.gtf",
                    expressed_kept_fname=args.output_prefix + ".expressed_kept.gtf",
                    expressed_transcripts=expressed_transcripts,
                    probability=args.fraction,
                    min_expressed_transcripts=args.min_transcripts_expressed,
                    min_expression=args.min_expression,
                    min_excluded_expression=args.min_excluded_expression)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
