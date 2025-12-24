#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

def load_simple_mapping(file):
    """Load mapping from key -> value (2-column TSV)."""
    mapping = {}
    with open(file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            mapping[row[0]] = row[1]
    return mapping

def load_spot_mapping(file):
    """Load mapping from 2um spot ID -> (8um, 16um)."""
    mapping = {}
    with open(file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            spot2um, spot8um, spot16um = row
            mapping[spot2um[:-2]] = {"8um": spot8um[:-2], "16um": spot16um[:-2]}
    return mapping

def parse_true_barcode(read_id):
    """Extract true barcode (part1|part2) from read_id (between 3rd and 4th underscore)."""
    parts = read_id.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected read_id format: {read_id}")
    if read_id.startswith('PacBio'):
        return parts[7]
    return parts[3]

def compute_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0.0
    return precision, recall

def get_spot_ids(barcode, part1_to_coord2, part2_to_coord1, spot_map, is_barcode_seq=True):
    """Convert barcode (p1|p2) into spot IDs at 2um, 8um, 16um."""
    if is_barcode_seq:
        try:
            p1, p2 = barcode.split("|")
        except ValueError:
            return None  # malformed

        if p1 not in part1_to_coord2 or p2 not in part2_to_coord1:
            return None

        coord1 = part2_to_coord1[p2]
        coord2 = part1_to_coord2[p1]
        spot2um = f"s_002um_{coord1}_{coord2}"
    else:
        spot2um = barcode

    if spot2um not in spot_map:
        return None

    return {
        "2um": spot2um,
        "8um": spot_map[spot2um]["8um"],
        "16um": spot_map[spot2um]["16um"],
    }

def evaluate(reads_file, part1_to_coord2, part2_to_coord1, spot_map, assignment_type='seq'):
    counts = {level: Counter() for level in ["2um", "8um", "16um"]}

    with open(reads_file, newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        c = 0

        for row in reader:
            read_id, pred_barcode = row[0], row[1]
            if pred_barcode == '*':
                for level in ["2um", "8um", "16um"]:
                    counts[level]["fn"] += 1
                continue

            true_barcode = parse_true_barcode(read_id)

            true_spots = get_spot_ids(true_barcode, part1_to_coord2, part2_to_coord1, spot_map)
            pred_spots = get_spot_ids(pred_barcode, part1_to_coord2, part2_to_coord1, spot_map,
                                      is_barcode_seq=(assignment_type == 'seq'))

            if not true_spots or not pred_spots:
                for level in ["2um", "8um", "16um"]:
                    counts[level]["fn"] += 1
                continue

            for level in ["2um", "8um", "16um"]:
                if true_spots[level] == pred_spots[level]:
                    counts[level]["tp"] += 1
                else:
                    counts[level]["fp"] += 1
                    c += 1
                    #if c < 20: print("%s %s %s %s %s" % (read_id, pred_barcode, true_spots[level], pred_spots[level], level))
                    #counts[level]["fn"] += 1

    results = {}
    for level in ["2um", "8um", "16um"]:
        tp = counts[level]["tp"]
        fp = counts[level]["fp"]
        fn = counts[level]["fn"]
        precision, recall = compute_metrics(tp, fp, fn)
        results[level] = {"precision": precision, "recall": recall, "tp": tp, "fp": fp, "fn": fn}
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate precision/recall at 2um, 8um, 16um spot ID levels."
    )
    parser.add_argument("reads_tsv", help="TSV with read_id and predicted barcodes")
    parser.add_argument("part1_to_coord2", help="TSV mapping part1 -> coord2")
    parser.add_argument("part2_to_coord1", help="TSV mapping part2 -> coord1")
    parser.add_argument("spot_map", help="TSV mapping 2um -> 8um,16um spot IDs")
    parser.add_argument("--assignment", choices=['spot', 'seq'], default='seq',
                        help="barcode assignment type [spot | seq]")
    args = parser.parse_args()

    part1_to_coord2 = load_simple_mapping(args.part1_to_coord2)
    part2_to_coord1 = load_simple_mapping(args.part2_to_coord1)
    spot_map = load_spot_mapping(args.spot_map)

    results = evaluate(args.reads_tsv, part1_to_coord2, part2_to_coord1, spot_map, args.assignment)

    print("Level\tPrecision\tRecall\tTP\tFP\tFN")
    for level in ["2um", "8um", "16um"]:
        r = results[level]
        print(f"{level}\t{r['precision']:.4f}\t{r['recall']:.4f}\t{r['tp']}\t{r['fp']}\t{r['fn']}")

if __name__ == "__main__":
    main()
