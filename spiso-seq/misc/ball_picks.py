import sys
import scipy.special
import random
from collections import defaultdict
import itertools
from traceback import print_exc

ITERATIONS = 300000


def simulate(picks, color_counts):
    balls = []
    color = 0
    for count in color_counts:
        for i in range(count):
            balls.append(color)
        color += 1
    #print(balls)

    pick_counts = defaultdict(int)
    for i in range(ITERATIONS):
        different_colors = len(set(random.sample(balls, picks)))
        pick_counts[different_colors] += 1

    #print(pick_counts)
    return pick_counts


def expected_num_from_simulation(pick_counts):
    total_sum = 0
    expected = 0
    for k in pick_counts.keys():
        v = pick_counts[k]
        total_sum += v
        expected += k * v
    return float(expected) / float(total_sum)


def subset_probability(picks, color_counts, color_ids):
    total_balls = sum(color_counts)
    interesting_balls = sum([color_counts[i] for i in color_ids])
    if interesting_balls < picks:
        return 0.0
    return float(scipy.special.perm(interesting_balls, picks)) / float(scipy.special.perm(total_balls, picks))


def simplified_formula(picks, color_counts):
    different_colors = len(color_counts)
    expected_colors = float(different_colors)
    index_list = [i for i in range(different_colors)]
    for combination in itertools.combinations(index_list, different_colors - 1):
        expected_colors -= subset_probability(picks, color_counts, combination)
    return expected_colors


def only_subset_probabilty(picks, color_counts, color_ids):
    sign = 1.0
    probability = 0.0
    for i in range(len(color_ids), 0, -1):
        for combination in itertools.combinations(color_ids, i):
            probability += sign * subset_probability(picks, color_counts, combination)
        sign *= -1.0
    return probability


def honest_analytical_formula(picks, color_counts):
    expected = 0.0
    index_list = [i for i in range(len(color_counts))]
    for i in range(1, len(color_counts) + 1):
        for combination in itertools.combinations(index_list, i):
            expected += i * only_subset_probabilty(picks, color_counts, combination)
    return expected


def read_input(inf):
    f = open(inf)
    picks = int(f.readline())
    color_counts = list(map(int, f.readline().split()))
    return picks, color_counts


def main():
    if len(sys.argv) < 2:
        print("Specifiy input file!")
        sys.exit(-1)

    picks, color_counts = read_input(sys.argv[1])
    print("Balls of different colors: " + ", ".join(map(str, color_counts)))
    print("Number of picks: %d" % picks)
    print("Honest analytical formula: %2.6f" % honest_analytical_formula(picks, color_counts))
    print("Simplified analytical formula: %2.6f" % simplified_formula(picks, color_counts))
    print("Simulation: %2.6f" % expected_num_from_simulation(simulate(picks, color_counts)))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


