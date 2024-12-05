#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import os
import sys
import argparse
from traceback import print_exc
from collections import defaultdict
from pyfaidx import Fasta
from Bio import SeqIO, Seq, SeqRecord


DELTA = 10

base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverse_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


def interval_len(interval):
    return interval[1] - interval[0] + 1


def intervals_total_length(sorted_range_list):
    total_len = 0
    for r in sorted_range_list:
        total_len += interval_len(r)
    return total_len


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def interval_list_overlaps(range_list1, range_list2):
    pos1 = 0
    pos2 = 0
    while pos1 < len(range_list1) and pos2 < len(range_list2):
        if overlaps(range_list1[pos1], range_list2[pos2]):
            return True
        elif range_list1[pos1][1] < range_list2[pos2][0]:
            pos1 += 1
        else:
            pos2 += 1

    return False


def contains(bigger_range, smaller_range):
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def get_interval_until(interval_list, pos):
    res = []
    for i in interval_list:
        if i[0] >= pos:
            break
        elif i[1] <= pos:
            res.append(i)
        else:
            res.append((i[0], pos))
    return res


def get_interval_from(interval_list, pos):
    res = []
    for i in interval_list:
        if i[0] >= pos:
            res.append(i)
        elif i[1] <= pos:
            continue
        else:
            res.append((pos, i[1]))
    return res


def get_trailing_unique(range_list1, range_list2):
    range1 = (range_list1[0][0], range_list1[-1][1])
    range2 = (range_list2[0][0], range_list2[-1][1])
    if not overlaps(range1, range2):
        print("Non-overlapping regions")
        return None, None

    if contains(range1, range2):
        utr1_extra = get_interval_until(range_list1, range2[0]) + get_interval_from(range_list1, range2[1])
        utr2_extra = None
    elif contains(range2, range1):
        utr2_extra = get_interval_until(range_list2, range1[0]) + get_interval_from(range_list2, range1[1])
        utr1_extra = None
    elif range1[0] < range2[0]:
        utr1_extra = get_interval_until(range_list1, range2[0])
        utr2_extra = get_interval_from(range_list2, range1[1])
    else:
        utr2_extra = get_interval_until(range_list2, range1[0])
        utr1_extra = get_interval_from(range_list1, range2[1])

    return utr1_extra, utr2_extra


def to_coors(reg_str):
    utrs = reg_str.split(";%;")
    if not utrs[0]: utrs = utrs[1:]
    chr_id = ""
    utr_coords = []
    strand = ""
    for u in utrs:
        v = u.split('_')
        chr_id = v[0]
        utr_coords.append((int(v[1]), int(v[2])))
        strand = v[3]

    return chr_id, utr_coords, strand


def load_genes(inf):
    gene_dict = {}
    for l in open(inf):
        v = l.split()
        gene_id = v[0]
        pa1 = v[6]
        pa2 = v[7]
        if pa2 == "NA": continue

        gene_dict[gene_id] = {pa1, pa2}
    return gene_dict


# 00000cbd-aa5d-4bf5-b193-3cead774ea12    ENSG00000240972.2       OldL5-6 GTGGGTCAGACTTT  ATGGAAGCA       ;%;chr22_23894583_23894771_+;%;chr22_23894945_23895039_+        chr22_23894402_23894402_+       chr22_23895224_23895224_+       ;%;chr22_23894402_23894582_+;%;chr22_23894772_23894944_+;%;chr22_23895040_23895223_+       known   2       ENST00000215754.8       protein_coding  chr22_23895103_23895103_+       121     ;%;chr22_23895104_23895224_+
def process_allinfo_with_utr(allinfo_utr, gene_polya_dict):
    gene_utr_counts = defaultdict(int)

    for l in open(allinfo_utr):
        v = l.split()
        gene_id = v[1]
        if gene_id not in gene_polya_dict: continue
        polya = v[7]
        if polya not in gene_polya_dict[gene_id]: continue
        group = "Old" if v[2].startswith("Old") else "Young"
        utr = v[-1]
        if utr == "*": continue
        gene_utr_counts[(gene_id, utr, group)] += 1

    return gene_utr_counts


def extrac_seq(chr_dict, chr_id, intervals, strand, start_delta=-1, end_delta=0):
    if chr_id not in chr_dict:
        return None
    seq = ""
    for interval in intervals:
        s = str(chr_dict[chr_id][interval[0]-start_delta:interval[1]+end_delta]).upper()
        if strand == '-':
            s = reverse_complement(s)
        seq += s
    return seq


def get_utr_classification(gene_utr_counts, chr_dict):
    gene_dict = defaultdict(list)
    for gene_id, utr, group in gene_utr_counts.keys():
        gene_dict[gene_id].append((utr, group))

    utr_sequences = defaultdict(list)
    for gene_id in gene_dict.keys():
        print("Processing gene %s" % gene_id)
        utrs = set([x[0] for x in gene_dict[gene_id]])
        if len(utrs) != 2:
            print("This gene has %d UTRs and will be ignored" % len(utrs))
            continue
        utrs = list(utrs)
        utr_str1 = utrs[0]
        chr_id, utr1, strand = to_coors(utr_str1)
        utr_str2 = utrs[1]
        utr2 = to_coors(utr_str2)[1]

        if gene_utr_counts[(gene_id, utr_str1, "Young")] > gene_utr_counts[(gene_id, utr_str1, "Old")]:
            utr1_dominant = "Young"
        else:
            utr1_dominant = "Old"
        if gene_utr_counts[(gene_id, utr_str2, "Young")] > gene_utr_counts[(gene_id, utr_str2, "Old")]:
            utr2_dominant = "Young"
        else:
            utr2_dominant = "Old"

        if not interval_list_overlaps(utr1, utr2):
            utr1_seq = extrac_seq(chr_dict, chr_id, utr1, strand)
            utr2_seq = extrac_seq(chr_dict, chr_id, utr2, strand)
            if utr1_seq:
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr1_seq), id="%s_%s_%s" % (gene_id, utr_str1, utr1_dominant), description="")
                utr_sequences[("no_overlap", utr1_dominant)].append(seq_rec)
            if utr2_seq:
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr2_seq), id="%s_%s_%s" % (gene_id, utr_str2, utr2_dominant), description="")
                utr_sequences[("no_overlap", utr2_dominant)].append(seq_rec)
        else:
            utr1_seq = extrac_seq(chr_dict, chr_id, utr1, strand)
            utr2_seq = extrac_seq(chr_dict, chr_id, utr2, strand)
            if utr1_seq:
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr1_seq), id="%s_%s_%s" % (gene_id, utr_str1, utr1_dominant), description="")
                utr_sequences[("overlap_full_UTR", utr1_dominant)].append(seq_rec)
            if utr2_seq:
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr2_seq), id="%s_%s_%s" % (gene_id, utr_str2, utr2_dominant), description="")
                utr_sequences[("overlap_full_UTR", utr2_dominant)].append(seq_rec)
            utr1_extra, utr2_extra = get_trailing_unique(utr1, utr2)
            if utr1_extra:
                utr1_extra_seq = extrac_seq(chr_dict, chr_id, utr1_extra, strand)
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr1_extra_seq), id="%s_%s_%s" % (gene_id, utr1_extra, utr1_dominant),
                                              description="")
                utr_sequences[("overlap_extra_UTR", utr1_dominant)].append(seq_rec)
            if utr2_extra:
                utr2_extra_seq = extrac_seq(chr_dict, chr_id, utr2_extra, strand)
                seq_rec = SeqRecord.SeqRecord(seq=Seq.Seq(utr2_extra_seq), id="%s_%s_%s" % (gene_id, utr2_extra, utr2_dominant),
                                              description="")
                utr_sequences[("overlap_extra_UTR", utr2_dominant)].append(seq_rec)

    return utr_sequences


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", required=True)
    parser.add_argument("--allinfo", "-a", type=str, help="input ALLINFO with UTRs", required=True)
    parser.add_argument("--gene_list", "-l", type=str, help="gene list with polyAs and FDRs", required=True)
    parser.add_argument("--genome", "-g", type=str, help="genome in FASTA", required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    gene_polya_dict = load_genes(args.gene_list)
    gene_utr_counts = process_allinfo_with_utr(args.allinfo, gene_polya_dict)
    reference_record_dict = Fasta(args.genome)
    utr_sequences = get_utr_classification(gene_utr_counts, reference_record_dict)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    for seq_type, dominant in utr_sequences.keys():
        out_fasta = os.path.join(args.output, "%s_%s.fasta" % (dominant, seq_type))
        SeqIO.write(utr_sequences[(seq_type, dominant)], out_fasta, 'fasta')



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
