############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from Bio import pairwise2
from Bio import Seq
from ssw import AlignmentMgr


window_size = 16
polyA_count = int(window_size * 0.75)


def find_polyt_start(seq):
    if len(seq) < window_size:
        return -1
    i = 0
    a_count = seq[0:window_size].count('T')
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


def align_pattern(sequence, start, end, pattern, min_score=0):
    seq1 = Seq.Seq(sequence[start:end])
    alignments = pairwise2.align.localms(pattern, seq1, 1, -1, -1.5, -1)
    alignment = max(alignments, key=lambda i: (i[2], -i[3]))
    if alignment[2] < min_score:
        return None, None
    return start + alignment[3], start + alignment[4] - 1


def align_pattern_ssw(sequence, start, end, pattern, min_score=0):
    seq = sequence[start:end]
    align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
    align_mgr.set_read(pattern)
    align_mgr.set_reference(seq)
    alignment = align_mgr.align(gap_open=1, gap_extension=1)
    if alignment.optimal_score < min_score:
        return None, None, None, None, None
    return start + alignment.reference_start, start + alignment.reference_end, \
        alignment.read_start, alignment.read_end, alignment.optimal_score


def find_candidate_with_max_score(barcode_matches, read_sequence, min_score=10):
    best_match = [0, 0, 0]
    best_barcode = None
    read_seq = Seq.Seq(read_sequence)
    for barcode in barcode_matches.keys():
        barcode_seq = Seq.Seq(barcode)
        alignments = pairwise2.align.localms(read_seq, barcode_seq, 1, -1, -1.5, -1)
        alignment = max(alignments, key=lambda i: (i[2], -i[3]))
        if alignment[2] < min_score:
            continue

        if alignment[2] > best_match[0]:
            best_barcode = barcode
            best_match[0] = alignment[2]
            best_match[1] = alignment[3]
            best_match[2] = alignment[4]
        elif alignment[2] == best_match[0] and alignment[3] < best_match[1]:
            best_barcode = barcode
            best_match[1] = alignment[3]
            best_match[2] = alignment[4]

    return best_barcode, best_match[1], best_match[2] - 1


def find_candidate_with_max_score_ssw(barcode_matches, read_sequence, min_score=10):
    best_match = [0, 0, 0]
    best_barcode = None
    align_mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
    align_mgr.set_reference(read_sequence)
    for barcode in barcode_matches.keys():
        align_mgr.set_read(barcode)
        alignment = align_mgr.align(gap_open=1, gap_extension=1)
        if alignment.optimal_score < min_score:
            continue

        if alignment.optimal_score > best_match[0]:
            best_barcode = barcode
            best_match[0] = alignment.optimal_score
            best_match[1] = alignment.reference_start
            best_match[2] = alignment.reference_end
        elif alignment.optimal_score  == best_match[0] and alignment.reference_start < best_match[1]:
            best_barcode = barcode
            best_match[1] = alignment.reference_start
            best_match[2] = alignment.reference_end

    return best_barcode, best_match[0], best_match[1], best_match[2]