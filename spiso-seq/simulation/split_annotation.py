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

def process_gene_db(db, main_gtf_fname, expressed_gtf_fname, excluded_gtf_fname, expressed_transcripts,
                    probability=0.05, min_expressed_transcripts=2,
                    min_expression=1,
                    min_excluded_expression=2):
    main_gtf = open(main_gtf_fname, "w")
    excl_gtf = open(excluded_gtf_fname, "w")
    expr_gtf = open(expressed_gtf_fname, "w")

    genes_kept = 0
    transcripts_kept = 0
    genes_reduced = 0
    transcripts_dropped = 0
    genes_expressed = 0
    transcripts_expressed = 0

    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        exon_count = {}
        expressed = set()
        for t in db.children(g, featuretype='transcript', order_by='start'):
            t_id = t.id.split('.')[0]
            if len(expressed_transcripts) > 0 and expressed_transcripts[t_id] < min_expression:
                continue
            expressed.add(t_id)
            exon_count[t_id] = 0
            for e in db.children(t, featuretype='exon', order_by='start'):
                exon_count[t_id] += 1

        to_remove = set()
        for t_id in exon_count.keys():
            if exon_count[t_id] == 1:
                to_remove.add(t_id)
        for t_id in to_remove:
            del exon_count[t_id]

        ignored_transcripts = set()
        if len(exon_count) > min_expressed_transcripts:
            for t_id in exon_count.keys():
                if expressed_transcripts[t_id] < min_excluded_expression:
                    continue
                r = random.random()
                if r < probability:
                    ignored_transcripts.add(t_id)

        gene_str = '%s\n' % g
        main_gtf.write(gene_str)
        genes_kept += 1
        if len(expressed) > len(ignored_transcripts):
            expr_gtf.write(gene_str)
            genes_expressed += 1
        if len(ignored_transcripts):
            excl_gtf.write(gene_str)
            genes_reduced += 1

        for t in db.children(g, featuretype='transcript', order_by='start'):
            transcript_str = '%s\n' % t
            t_id = t.id.split('.')[0]
            if t_id in ignored_transcripts:
                current_file = excl_gtf
                transcripts_dropped += 1
            else:
                current_file = main_gtf
                transcripts_kept += 1
            current_file.write(transcript_str)
            for e in db.children(t, featuretype='exon', order_by='start'):
                current_file.write('%s\n' % e)

            if t_id in expressed:
                transcripts_expressed += 1
                expr_gtf.write(transcript_str)
                for e in db.children(t, featuretype='exon', order_by='start'):
                    expr_gtf.write('%s\n' % e)


    main_gtf.close()
    excl_gtf.close()
    expr_gtf.close()   
    print("Expressed genes: %d, transcripts: %d" % (genes_expressed, transcripts_expressed))
    print("Kept genes: %d, transcripts: %d" % (genes_kept, transcripts_kept))
    print("Reduced genes: %d, dropped transcripts: %d" % (genes_reduced, transcripts_dropped))


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
                        type=int, default=2)
    parser.add_argument("--min_excluded_expression", help="minimal number reads in transcript to be possibly excluded",
                        type=int, default=4)
    parser.add_argument("--min_transcripts_expressed", help="minimal number of transcripts expressed in a gene to "
                                                            "consider them for dropping", type=int, default=2)

    parser.add_argument("--fraction", help="fraction of non-monoexonic transcripts to be removed",
                        type=float, default=0.05)
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
        expressed_transcripts = read_expressed_isoforms(args.expressed)
        expressed_fname=args.output_prefix + ".expression.tsv"
        print("Saving expression table to " + expressed_fname)
        f = open(expressed_fname, 'w')
        for k in expressed_transcripts:
            f.write("%s\t%d\n" % (k, expressed_transcripts[k]))
        f.close()

    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    process_gene_db(gffutils_db, main_gtf_fname=args.output_prefix + ".reduced.gtf",
                    expressed_gtf_fname=args.output_prefix + ".expressed.gtf",
                    excluded_gtf_fname=args.output_prefix + ".excluded.gtf",
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
