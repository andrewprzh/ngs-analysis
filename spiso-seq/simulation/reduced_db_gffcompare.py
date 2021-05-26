#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import subprocess
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
import gffutils
from enum import Enum, unique


@unique
class TranscriptType(Enum):
    known = 1
    novel = 2
    undefined = 0


def separate_isoquant(l):
    if l.find("known"):
        return TranscriptType.known
    elif l.find("nic"):
        return TranscriptType.novel
    return TranscriptType.undefined


def separate_stringtie(l):
    if l.find("reference_id"):
        return TranscriptType.known
    else:
        return TranscriptType.novel


def separate_tid(l):
    if l.find('transcript_id "ENSMUST'):
        return TranscriptType.known
    else:
        return TranscriptType.novel


SEPARATE_FUNCTORS = {'isoquant':separate_isoquant,
                     'stringtie':separate_stringtie,
                     'flair':separate_tid,
                     'talon':separate_tid,
                     'bambu':separate_tid}


def split_gtf(ingtf_path, seaprate_functor, out_known_path, out_novel_path):
    out_known = open(out_known_path, "w")
    out_novel = open(out_novel_path, "w")
    for l in open(ingtf_path):
        if l.startswith("#"):
            continue
        ttype = seaprate_functor(l)
        if ttype == TranscriptType.novel:
            out_novel.write(l)
        elif ttype == TranscriptType.known:
            out_known.write(l)

    out_novel.close()
    out_known.close()


def run_gff_compare(reference_gtf, compared_gtf, output):
    result = subprocess.run(["gffcompare", "-r", reference_gtf, "-o", output, compared_gtf])

    if result.returncode != 0:
        print("gffcompare faile for " + compared_gtf)
        return


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", default="gtf_stats")
    parser.add_argument("--genedb", "-d", type=str, help="prefix to reduced gene db")
    parser.add_argument("--gtf", "-g", type=str, help="output gtf")
    parser.add_argument("--tool", type=str, choices=['isoquant', 'talon', 'sqanti', 'flair', 'bambu', 'stringtie'],
                        help="tool used for generating GTF")

    args = parser.parse_args()
    if not args.genedb or not args.gtf or not args.tool:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    out_known_path = os.path.join(args.output, args.tool + ".known.gtf")
    out_novel_path= os.path.join(args.output, args.tool + ".known.gtf")
    print("Seprating known and novel transcripts")
    split_gtf(args.gtf, SEPARATE_FUNCTORS[args.tool], out_known_path, out_novel_path)
    print("Running gffcompare for entire GTF")
    expressed_gtf = args.genedb + ".expressed.gtf"
    run_gff_compare(expressed_gtf, args.gtf, os.path.join(args.output, args.tool + ".full.stats"))
    print("Running gffcompare for known transcripts")
    expressed_gtf = args.genedb + ".expressed_kept.gtf"
    run_gff_compare(expressed_gtf, out_known_path, os.path.join(args.output, args.tool + ".known.stats"))
    print("Running gffcompare for novel transcripts")
    expressed_gtf = args.genedb + ".excluded.gtf"
    run_gff_compare(expressed_gtf, out_novel_path, os.path.join(args.output, args.tool + ".novel.stats"))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
