############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from Bio import pairwise2
from Bio import Seq


window_size = 16
polyA_count = int(window_size * 0.75)


def find_polyt_start(seq):
    if len(seq) < window_size:
        return -1
    i = 0
    a_count = seq[0:window_size].count('A')
    while i < len(seq) - window_size:
        if a_count >= polyA_count:
            break
        first_base_a = seq[i] == 'T'
        new_base_a = i + window_size < len(seq) and seq[i + window_size] == 'T'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        i += 1

    if i >= len(seq) - window_size:
        return -1

    return i + max(0, seq[i:].find('TT'))


base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}


def reverese_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)


def align_pattern(sequence, start, end, pattern):
    seq1 = Seq.Seq(sequence[start:end])
    alignments = pairwise2.align.localms(seq1, pattern, 2, -2, -3, -2)  # penalty for gaps
    alignment = max(alignments, key=lambda i: (i.score, -i.start))
    return start + alignment.start, start + alignment.end


def find_candidate_with_max_score(barcode_matches, read_sequence, min_score=10):
    best_match = [0, 0, 0]
    best_barcode = None
    read_seq = Seq.Seq(read_sequence)
    for barcode in barcode_matches.keys():
        barcode_seq = Seq.Seq(barcode)
        alignments = pairwise2.align.localms(read_seq, barcode_seq, 1, -1, -1.5, -1)
        alignment = max(alignments, key=lambda i: (i.score, -i.start))
        if alignment.score < min_score:
            continue

        if alignment.score > best_match[0]:
            best_barcode = barcode
            best_match[0] = alignment.score
            best_match[1] = alignment.start
            best_match[2] = alignment.end
        elif alignment.score == best_match[0] and alignment.start < best_match[1]:
            best_barcode = barcode
            best_match[1] = alignment.start
            best_match[2] = alignment.end

    return best_barcode, best_match[1], best_match[2]