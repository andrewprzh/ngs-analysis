############################################################################
# Authors: Hagen Tilgner, Andrey Prjibelski
############################################################################

import os
import sys
import math
from traceback import print_exc
from functools import partial

EULER_C = 0.58
# can be done in linear time
HARMONIC_NUMBERS = [sum([1 / float(j) for j in range(1, i + 1)]) for i in range(1, 100)]


# introns -- ordered list of coordinate pairs
def detect_matching_intron_chains(predicted_introns, reference_introns):
    pred_pos = 0
    ref_pos = 0
    # list of pairs, indicating 0-based indices of matching introns
    # e.g. [(1,3), (5,5)] means that introns 1,2,3 and 5 from predicted_introns match reference_introns
    matching_chains = []
    current_matching_chain = None

    while pred_pos < len(predicted_introns) and ref_pos < len(reference_introns):
        if predicted_introns[pred_pos] == reference_introns[ref_pos]:
            # equal introns
            if current_matching_chain is None:
                current_matching_chain = (pred_pos, pred_pos)
            else:
                current_matching_chain = (current_matching_chain[0], pred_pos)
            pred_pos += 1
            ref_pos += 1
        elif predicted_introns[pred_pos][1] < reference_introns[ref_pos][1]:
            # predicted_intron is on the left from reference_intron
            if current_matching_chain is not None:
                matching_chains.append(current_matching_chain)
                current_matching_chain = None
            pred_pos += 1
        else:
            # reference_intron is on the left from predicted_intron
            if current_matching_chain is not None:
                matching_chains.append(current_matching_chain)
                current_matching_chain = None
            ref_pos += 1

    if current_matching_chain is not None:
        matching_chains.append(current_matching_chain)

    return matching_chains


# count matching intron chain using Hagen's scheme
# that adds 1/n to the score for every matching sub-chain on length n
# thus, for a chain of length 4 it gives score:
# (1 + 1 + 1 + 1) + (1/2 + 1/2 + 1/2) + (1/3 + 1/3) + 1/4 = 6.42
#
# basically, this scheme scores each chain of length N
# as (N + 1) * (H_N - 1), where H_N is N-th Harmonic number, i.e.
# H_N = sum (1 + 1/2 + ... + 1/N)
# furthermore,  H_N be approximated by ln N + 1/2N + C function, # where C is Euler's constant
# see simplified scoring schemes below
def score_hagen(matching_chains, min_chain_len = 1):
    total_score = 0.0
    # consider all matching chains
    for chain in matching_chains:
        matching_chain_len = chain[1] - chain[0] + 1
        chain_score = 0.0
        for i in range(min_chain_len, matching_chain_len + 1):
            # for each sub-chain of length i add a score of 1/i
            chain_score += float(matching_chain_len - i + 1) / float(i)

        total_score += chain_score

    return total_score


# count using pre-computed Harmonic numbers
def score_hagen_fast(matching_chains):
    total_score = 0.0
    # consider all matching chains
    for chain in matching_chains:
        matching_chain_len = chain[1] - chain[0] + 1
        harmonic_number = HARMONIC_NUMBERS[matching_chain_len]
        total_score += float(matching_chain_len + 1) * (harmonic_number - 1)

    return total_score


# count n log n approximation
def score_nlogn(matching_chains):
    total_score = 0.0
    # consider all matching chains
    for chain in matching_chains:
        matching_chain_len = chain[1] - chain[0] + 1
        total_score += float(matching_chain_len + 1) * (math.log(matching_chain_len) - 1 + EULER_C + 1 / float(matching_chain_len * 2))

    return total_score


# count quadratic approximation
def score_quadratic(matching_chains):
    total_score = 0.0
    # consider all matching chains
    for chain in matching_chains:
        matching_chain_len = float(chain[1] - chain[0] + 1)
        total_score += matching_chain_len * matching_chain_len

    return total_score


def score_prediction(predicted_introns, reference_introns):
    # best possible predicted chain is all introns matched
    best_possible_chain = [(0, len(reference_introns) - 1)]
    matching_chains = detect_matching_intron_chains(predicted_introns, reference_introns)
    print("Matching chains are " + str(matching_chains))

    hagen_score = score_hagen(matching_chains)
    max_hagen_score = score_hagen(best_possible_chain)
    print("Hagen absolute score: {0:.2f}, normalized Hagen score: {1:.2f}".format(hagen_score, hagen_score / max_hagen_score))

    hagen_score = score_hagen_fast(matching_chains)
    max_hagen_score = score_hagen_fast(best_possible_chain)
    print("Fast Hagen absolute score: {0:.2f}, normalized fast Hagen score: {1:.2f}".format(hagen_score, hagen_score / max_hagen_score))

    nlogn_score = score_nlogn(matching_chains)
    max_nlogn_score  = score_nlogn(best_possible_chain)
    print("nlogn absolute score: {0:.2f}, normalized nlong score: {1:.2f}".format(nlogn_score, nlogn_score / max_nlogn_score))

    q_score = score_quadratic(matching_chains)
    max_q_score = score_quadratic(best_possible_chain)
    print("Quadratic absolute score: {0:.2f}, normalized quadratic score: {1:.2f}".format(q_score, q_score / max_q_score))


# parse a line from file
def parse_introns(line):
    # introns -- ordered list of coordinate pairs
    introns = []
    coordinate_pairs = line.strip().split(',')
    for pair in coordinate_pairs:
        coordinates = list(map(int, pair.split('-')))
        assert len(coordinates) == 2
        assert coordinates[0] <= coordinates[1]
        if len(introns) > 0:
            assert introns[-1][1] < coordinates[0]

        introns.append((coordinates[0], coordinates[1]))

    return introns


def process_file(fname):
    inf = open(fname)

    predicted_introns = None
    l = inf.readline()
    while l:
        if l.startswith("#") or len(l.strip()) == 0:
            l = inf.readline()
            continue

        if predicted_introns == None:
            # if nothing is parsed yet
            predicted_introns = parse_introns(l)
        else:
            # if predicted introns are already scanned
            reference_introns = parse_introns(l)
            score_prediction(predicted_introns, reference_introns)
            predicted_introns = None

        l = inf.readline()

    inf.close()


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Provide input file with introns\n")
        sys.exit(-1)

    process_file(sys.argv[1])


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)