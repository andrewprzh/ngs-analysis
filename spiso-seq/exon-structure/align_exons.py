import sys
import os
import argparse
from Bio import SeqIO
import gzip
from Bio import pairwise2
from collections import defaultdict
import numpy
from traceback import print_exc


class ExonPair:
    def __init__(self, chr_id, strand, exon1, exon2, exon_type="", gene_id=""):
        self.chr_id = chr_id
        self.stand = strand
        self.exon1 = exon1
        self.exon2 = exon2
        self.exon_type = exon_type
        self.gene_id = gene_id
        self.global_similarity = 0.0
        self.local_similarity = 0.0


def load_exon_pairs(inf):
    print("Loading exon pairs")
    exon_storage = defaultdict(list)
    for l in open(inf):
        t = l.strip().split()
        exon1_t = t[0].split("_")
        exon2_t = t[0].split("_")
        assert exon2_t[0] == exon1_t[0]
        assert exon2_t[3] == exon1_t[3]
        chr_id = exon1_t[0]
        exon1 = (int(exon1_t[1]), int(exon1_t[2]))
        exon2 = (int(exon2_t[1]), int(exon2_t[2]))
        exon_pair = ExonPair(chr_id, exon1_t[3], exon1, exon2, t[3], t[2])
        exon_storage[chr_id].append(exon_pair)
    return exon_storage


class ExonComparator:
    def __init__(self, args):
        self.args = args
        print("Loading reference genome from " + self.args.reference)
        _, outer_ext = os.path.splitext(self.args.reference)
        if outer_ext.lower() in ['.gz', '.gzip']:
            with gzip.open(self.args.reference, "rt") as handle:
                self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        else:
            self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(self.args.reference, "fasta"))

    def map_exon_pairs(self, exon_storage):
        print("Aligning exon pairs")
        for chr_id, exon_list in exon_storage.items():
            chr_record = self.reference_record_dict[chr_id]
            for exon_pair in exon_list:
                exon1_seq = chr_record[exon_pair.exon1[0]:exon_pair.exon1[1] + 1]
                exon2_seq = chr_record[exon_pair.exon2[0]:exon_pair.exon2[1] + 1]
                exon_min_len = min(len(exon1_seq), len(exon2_seq))
                exon_pair.global_similarity = pairwise2.align.globalxx(exon1_seq, exon2_seq, score_only = True) / exon_min_len
                exon_pair.local_similarity = pairwise2.align.localxx(exon1_seq, exon2_seq, score_only=True) / exon_min_len

    def compute_hist(self, exon_storage, property_func=lambda x:x.global_similarity):
        print("Computing histograms")
        exon_similarity_map = defaultdict(list)
        for chr_id, exon_list in exon_storage.items():
            for exon_pair in exon_list:
                exon_similarity_map[exon_pair.exon_type].append(property_func(exon_pair))

        hist_map = {}
        bins = [0.1 * i for i in range(11)]
        for exon_type, similarity_values in exon_similarity_map.items():
            hist_map[exon_type] = numpy.histogram(similarity_values, bins)[0]
        return hist_map


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", "-r", help="reference genome", type=str)
    parser.add_argument("--exon_pairs", "-e", help="exon pairs in 4-column TSV", type=str)
    parser.add_argument("--output", "-o", help="output file", type=str, default="exon_similarity.tsv")
    args = parser.parse_args()
    if not args.reference or not args.exon_pairs:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    exon_storage = load_exon_pairs(args.exon_pairs)
    exon_comparator = ExonComparator(args)
    exon_comparator.map_exon_pairs(exon_storage)
    outf = open(args.output, "w")
    hist_map = exon_comparator.compute_hist(exon_storage, property_func=lambda x:x.global_similarity)
    outf.write("Global similarity histograms")
    for exon_type, hist in hist_map:
        outf.write(exon_type + "\t" + "\t".join([str(x) for x in hist]))
    hist_map = exon_comparator.compute_hist(exon_storage, property_func=lambda x: x.local_similarity)
    outf.write("Local similarity histograms")
    for exon_type, hist in hist_map:
        outf.write(exon_type + "\t" + "\t".join([str(x) for x in hist]))
    print("Done")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)








