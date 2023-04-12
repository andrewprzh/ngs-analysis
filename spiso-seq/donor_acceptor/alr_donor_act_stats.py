#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import gffutils
from collections import defaultdict
import sys
import argparse
from traceback import print_exc
from glob import glob

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix")
    parser.add_argument("--genedb", "-g", type=str, help="GTF converted to .db (can be found in IsoQuant output)", required=True)
    parser.add_argument("--isoquant_output", "-i", type=str, nargs="+", help="IsoQuant output folders (must contain intron_coutns.tsv)")

    args = parser.parse_args()
    return args


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def detect_alt_introns(genedb, chrid):
    start_pos_map = defaultdict(set)
    end_pos_map = defaultdict(set)
    intron_to_gene_dict = {}

    for g in genedb.region(seqid=chrid, featuretype="gene"):
        for t in genedb.children(g, featuretype=('transcript', 'mRNA'), order_by='start'):
            exons = []
            for e in genedb.children(t, order_by='start'):
                if e.featuretype == 'exon':
                    exons.append((e.start, e.end))

            introns = junctions_from_blocks(exons)
            for i in range(len(introns)):
                intron = introns[i]
                if i > 0:
                    prev_intron = introns[i-1]
                    end_pos_map[(intron[1], prev_intron[1])].add((intron[0], intron[1], t.strand))
                if i < len(introns) - 1:
                    next_intron = introns[i + 1]
                    start_pos_map[(intron[0], next_intron[0])].add((intron[0], intron[1], t.strand))
                intron_to_gene_dict[(intron[0], intron[1], t.strand)] = g.id

    alt_donor_dict = defaultdict(set)
    alt_acceptor_dict = defaultdict(set)

    for pos in end_pos_map.keys():
        if len(end_pos_map[pos]) <= 1:
            continue
        plus_strand_introns = set()
        minus_strand_introns = set()
        for i in end_pos_map[pos]:
            if i[2] == "+":
                plus_strand_introns.add(i)
            else:
                minus_strand_introns.add(i)
        if len(plus_strand_introns) > 1:
            for i in plus_strand_introns:
                alt_donor_dict[i] = plus_strand_introns
        if len(minus_strand_introns) > 1:
            for i in plus_strand_introns:
                alt_acceptor_dict[i] = plus_strand_introns

    for pos in start_pos_map.keys():
        if len(start_pos_map[pos]) <= 1:
            continue
        plus_strand_introns = set()
        minus_strand_introns = set()
        for i in start_pos_map[pos]:
            if i[2] == "+":
                plus_strand_introns.add(i)
            else:
                minus_strand_introns.add(i)
        if len(plus_strand_introns) > 1:
            for i in plus_strand_introns:
                alt_acceptor_dict[i] = plus_strand_introns
        if len(minus_strand_introns) > 1:
            for i in start_pos_map[pos]:
                alt_donor_dict[i] = start_pos_map[pos]

    return alt_donor_dict, alt_acceptor_dict, intron_to_gene_dict


class IntronCountIterator:
    def __init__(self, inf_path):
        self.inf_path = inf_path
        self.inf = open(inf_path, "r")
        self.current_chr_id = ""

    def __del__(self):
        self.inf.close()

    def next(self):
        l = self.inf.readline()
        if not l:
            return None, None, None
        if l.startswith("#"):
            return self.next()

        v = l.split("\t")
        chr_id = v[0]
        if self.current_chr_id != chr_id:
            self.current_chr_id = chr_id
            return chr_id, None, None
        return chr_id, (int(v[1]), int(v[2]), v[3]), int(v[7])


def process_chromosome(genedb, chrid, iterator_list, donor_outf, acc_outf):
    alt_donor_dict, alt_acceptor_dict, intron_to_gene_dict = detect_alt_introns(genedb, chrid)
    #print(alt_donor_dict)
    #print(alt_acceptor_dict)
    alt_donor_counts = {}
    alt_acceptor_counts = {}
    n_samples = len(iterator_list)
    next_chr = None
    for i, it in enumerate(iterator_list):
        intron = [0, 0, 0]
        while intron[1] is not None:
            intron = it.next()
            if chrid != intron[0]:
                if not next_chr:
                    next_chr = intron[0]
                else:
                    assert next_chr == intron[0]
                break
            if intron[1] in alt_donor_dict:
                if intron[1] not in alt_donor_counts:
                    alt_donor_counts[intron[1]] = [0] * n_samples
                alt_donor_counts[intron[1]][i] += intron[2]
            if intron[1] in alt_acceptor_dict:
                if intron[1] not in alt_acceptor_counts:
                    alt_acceptor_counts[intron[1]] = [0] * n_samples
                alt_acceptor_counts[intron[1]][i] = intron[2]

    dump_introns_for_chrid(chrid, alt_donor_dict, alt_donor_counts, intron_to_gene_dict, donor_outf)
    dump_introns_for_chrid(chrid, alt_acceptor_dict, alt_acceptor_counts, intron_to_gene_dict, acc_outf)
    return next_chr


def dump_introns_for_chrid(chrid, alt_dict, count_dict, intron_to_gene_dict, outf):
    processed_introns = set()
    group_id = 0
    for intron in alt_dict:
        if intron in processed_introns:
            continue
        for i in alt_dict[intron]:
            processed_introns.add(i)
            if i not in count_dict:
                continue
            count_data = "\t".join(list(map(str, count_dict[i])))
            gid = intron_to_gene_dict[i]
            outf.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\n" % (chrid, chrid + "_" + str(group_id), gid, i[0], i[1], i[2], count_data))
        group_id += 1

def main():
    args = parse_args()
    print("Loading gene DB from %s" % args.genedb)
    genedb = gffutils.FeatureDB(args.genedb)
    iterator_list = []
    sample_list = []
    current_chrid = ""
    for iq_out in args.isoquant_output:
        intron_count_files = glob(os.path.join(iq_out, "*.intron_counts.tsv"))
        sample_list.append(os.path.basename(os.path.normpath(iq_out)))
        assert len(intron_count_files) == 1
        print("Found %s" % intron_count_files[0])
        iterator_list.append(IntronCountIterator(intron_count_files[0]))
        if not current_chrid:
            current_chrid = iterator_list[-1].next()[0]
        else:
            assert current_chrid == iterator_list[-1].next()[0]

    # print(sample_list)
    donor_outf = open(args.output + "_alt_donors.tsv", "w")
    acc_outf = open(args.output + "_alt_acceptors.tsv", "w")
    header = "#chrid\tgroup\tgene\tdonor\tacceptor\tstrand\t%s\n" % "\t".join(sample_list)
    donor_outf.write(header)
    acc_outf.write(header)
    while current_chrid is not None:
        print("Processing chromosome %s" % str(current_chrid))
        current_chrid = process_chromosome(genedb, current_chrid, iterator_list, donor_outf, acc_outf)
    donor_outf.close()
    acc_outf.close()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
