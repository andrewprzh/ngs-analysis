import logging
import sys
import os
from Bio import Seq
from traceback import print_exc
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from Bio import pairwise2
from kmer_indexer import KmerIndexer
from collections import defaultdict

logger = logging.getLogger('IsoQuant')
BARCODE_LENGTH = 16
UMI_LENGTH = 12
PARTIAL_R1 = "CTACACGACGCTCTTCCGATCT"
PARTIAL_R1_LENGTH = 22
NOT_FOUND = "not_found"
range_size = 10


class Finder:
    def __init__(self, window_size=20, min_polya_fraction=0.8):
        self.window_size = window_size
        self.min_polya_fraction = min_polya_fraction
        self.polyA_count = int(self.window_size * self.min_polya_fraction)
        self.partial_R1 = Seq(revComp(PARTIAL_R1))
        self.partialR1_indexer = KmerIndexer([PARTIAL_R1], kmer_size=4)

    # poly A tail detection
    def find_polya(self, seq):
        if len(seq) < self.window_size:
            return -1
        i = 0
        a_count = seq[0:self.window_size].count('A')
        while i < len(seq) - self.window_size:
            if a_count >= self.polyA_count:
                break
            first_base_a = seq[i] == 'A'
            new_base_a = i + self.window_size < len(seq) and seq[i + self.window_size] == 'A'
            if first_base_a and not new_base_a:
                a_count -= 1
            elif not first_base_a and new_base_a:
                a_count += 1
            i += 1

        if i >= len(seq) - self.window_size:
            return -1

        return i + max(0, seq[i:].find('AA'))

    def find_partial_R1(self, sequence, start):
        seq1 = Seq(sequence[start:])
        alignments = pairwise2.align.localms(seq1, self.partial_R1, 1, -1, -1, -1)  # penalty for gaps
        alignment = max(alignments, key=lambda i: (i.score, -i.start))
        return start + alignment.start

    def find_partial_R1_kmer(self, sequence, start, min_kmers=2):
        offset = 23  # barcode + umi length approximately
        seq = sequence[start + offset:]
        str_list = [seq[i:i + PARTIAL_R1_LENGTH] for i in range(20)]

        best_match = -1
        pos = -1
        for i, s in enumerate(str_list):
            r1_dict = self.partialR1_indexer.get_occurrences(s)
            if len(r1_dict) > 0 and r1_dict[0][1] >= min_kmers and r1_dict[0][1] > best_match:
                pos = i
                best_match = r1_dict[0][1]

        if best_match == -1:
            return -1, -1
        return best_match, start + offset + pos


def revComp(my_seq):  ## obtain reverse complement of a sequence
    base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}
    lms = list(my_seq[:-1])  ## parse string into list of components
    lms.reverse()
    try:
        lms = [base_comp[base] for base in lms]
    except TypeError:
        pass
    lms = ''.join(lms)  ## get string from reversed and complemented list
    return lms


# Given start and end position of scanning sequence, look for special strings using kmers
def kmer_scan(kmer_indexer, start_scan, end_scan, seq, str_length=BARCODE_LENGTH + 2, min_matches=3, gap=1):
    str_list = [seq[i:i + str_length] for i in range(start_scan, end_scan)]
    candidates = set()
    candidates_with_score = list()
    max_score = 0
    for s in set(str_list):
        kmer_dict = kmer_indexer.get_occurrences(s)
        if len(kmer_dict) == 0:
            continue
        score = kmer_dict[0][1]

        if score >= min_matches and score > max_score:
            max_score = score
        for cws in kmer_dict:
            if cws[1] >= min_matches and cws[1] >= max_score - gap:
                candidates_with_score.append(cws)

    candidates_with_score.sort(key=lambda elem: elem[1], reverse=True)
    for cws in candidates_with_score:
        if cws[1] >= max_score - gap:
            candidates.add(cws[0])

    return candidates


def find_candidate_with_max_score(candidates, s, min_score=13):
    best = [0, 0]
    found = NOT_FOUND
    seq1 = Seq(s)
    for current in candidates:
        seq2 = Seq(current)
        alignments = pairwise2.align.localms(seq1, seq2, 1, -1, -1.5, -0.1)
        alignment = max(alignments, key=lambda i: (i.score, -i.start))
        if alignment.score > best[0]:
            found = current
            best[0] = alignment.score
            best[1] = alignment.start
        elif alignment.score == best[0] and alignment.start < best[1]:
            found = current
            best[1] = alignment.start

    if best[0] < min_score:
        return NOT_FOUND, -1
    return found, best[1]


def get_file_name(input_file_name):
    base = os.path.basename(input_file_name)
    return os.path.splitext(base)[0]


def process(input_file_name, output_dir, data):
    umi_indexers, barcode_and_umi_to_gene, read_to_gene = data
    finder = Finder()
    if input_file_name[-1] == 'q':
        fasta_sequences = SeqIO.parse(open(input_file_name), 'fastq')
    else:
        fasta_sequences = SeqIO.parse(open(input_file_name), 'fasta')
    file_name = get_file_name(input_file_name)
    output_file = output_dir + "/" + file_name + ".tsv"
    barcode_indexer = KmerIndexer(barcodes, kmer_size=5)  # allows two mistakes

    using_gene = False
    with open(output_file, 'w') as out_file:
        for fasta in tqdm(fasta_sequences, total=200000, desc="processing fasta sequences"):
            name, sequence = fasta.id, str(fasta.seq)
            start = finder.find_polya(sequence)

            if start == -1:
                sequence = revComp(sequence)
                start = finder.find_polya(sequence)
                if start == -1:
                    out_file.write(name + "\t" + NOT_FOUND + "\t" + NOT_FOUND + "\tbad_polyA" + "\n")
                    continue

            start += 30
            best_match, partial_start = finder.find_partial_R1_kmer(sequence, start)  # allow 3 possible mistakes
            if partial_start == -1:
                out_file.write(name + "\t" + NOT_FOUND + "\t" + NOT_FOUND + "\tbad_partial" + "\n")
                continue

            hypothetical_barcode_start = partial_start - BARCODE_LENGTH
            offset = 4
            bc_candidates = kmer_scan(barcode_indexer, hypothetical_barcode_start - offset,
                                      hypothetical_barcode_start + offset, sequence)

            found_bc, found_bc_pos = find_candidate_with_max_score(bc_candidates,
                                                                   sequence[hypothetical_barcode_start - offset:
                                                                            hypothetical_barcode_start + BARCODE_LENGTH + offset])
            found_bc_pos += hypothetical_barcode_start - offset

            ## if we found barcode we will try to find UMI
            if found_bc == NOT_FOUND:
                out_file.write(name + "\t" + NOT_FOUND + "\t" + NOT_FOUND + "\n")
            else:
                offset = 5
                umi_candidates = kmer_scan(umi_indexers[found_bc], start - offset,
                                           found_bc_pos - UMI_LENGTH - 1, sequence, str_length=UMI_LENGTH + 1, gap=3)

                if using_gene:
                    found_umi = NOT_FOUND
                    if name not in read_to_gene:
                        break
                    gene_id = read_to_gene[name]

                    for umi_cnd in umi_candidates:
                        key = get_key(found_bc, umi_cnd)
                        if key in barcode_and_umi_to_gene:
                            if barcode_and_umi_to_gene[key] == gene_id:
                                found_umi = umi_cnd
                                break
                else:
                    found_umi, _ = find_candidate_with_max_score(umi_candidates,
                                                                 sequence[start - offset:
                                                                          start + UMI_LENGTH + offset])

                out_file.write(name + "\t" + found_bc + "\t" + found_umi + "\n")


def prelim(args):
    global barcodes

    barcode_and_umi_to_gene = dict()
    barcode_to_umi_indexer = dict()
    read_to_gene = dict()
    barcode_to_umi_set = defaultdict(set)

    gene_file = args.geneFile.replace(u'\xa0', u'')
    with open(gene_file) as f:
        for line in tqdm(f, total=10000, desc="reading reads to genes map"):
            read_id, gene_id = line.strip('\n').split('\t')[:2]
            read_to_gene[read_id] = gene_id

    sd_file = args.sdFile.replace(u'\xa0', u'')
    barcodes = set()

    with open(sd_file) as f:
        for line in tqdm(f, total=68205160, desc="reading spatial data"):
            gene_id, barcode, umi = line.strip('\n').split('\t')
            barcodes.add(barcode)
            barcode_to_umi_set[barcode].add(umi)
            barcode_and_umi_to_gene[get_key(barcode, umi)] = gene_id

    for key, umis in tqdm(barcode_to_umi_set.items(), desc="creating indexers"):
        barcode_to_umi_indexer[key] = KmerIndexer(umis, kmer_size=7)



    return barcode_to_umi_indexer, barcode_and_umi_to_gene, read_to_gene


def get_key(str1, str2):
    return str1 + '-' + str2


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fasta', type=str, help='fasta file')
    parser.add_argument('sdFile', type=str, help='genes, barcodes and umis')
    parser.add_argument('geneFile', type=str, help='read ids to gene ids map')
    parser.add_argument('--outDir', default="./", type=str, help='directory to put output in')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    data = prelim(args)
    process(args.fasta, args.outDir, data)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
