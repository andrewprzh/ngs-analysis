import sys
import pysam
from traceback import print_exc
import numpy as np


def main():
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
    left_truncations = []
    right_truncations = []
    for alignment in bamfile:
        if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
              continue

        transcript_len = bamfile.get_reference_length(bamfile.get_reference_name(alignment.reference_id))
        left_truncations.append(float(alignment.reference_start) / float(transcript_len))
        right_truncations.append(float(transcript_len - alignment.reference_end) / float(transcript_len))

    coverage_hist = [0] * 100
    for i in range(len(left_truncations)):
        for pos in range(int(round(left_truncations[i] * 100)), int(round((1 - right_truncations[i]) * 100))):
            coverage_hist[pos] += 1

    print("Coverage hist")
    print("\t".join(map(str, coverage_hist)))
    print("")

    bins = [0.0] + [0.005 + 0.01 * i for i in range(0,100)] + [1.0]
    print(len(bins), bins)
    res_left = np.histogram(left_truncations, bins=bins, density=True)
    probs_l = []
    for i, w in enumerate(res_left[0]):
        probs_l.append(w * (bins[i+1] - bins[i]))
    print(sum(probs_l), probs_l)

    res_right = np.histogram(right_truncations, bins=bins, density=True)
    probs_r = []
    for i, w in enumerate(res_right[0]):
        probs_r.append(w * (bins[i+1] - bins[i]))

    print(sum(probs_r), probs_r)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
