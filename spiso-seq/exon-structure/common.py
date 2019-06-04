import os
import sys
import pysam


dbg_barcodes = set()

MIN_READS = 10
MIN_READS_FOR_EXON = 2
RELATIVE_COUNT_COEFF = 5

def debug(bc, s):
    if bc in dbg_barcodes:
        print(s)


class ReadRecord:
    block_coordinates = []
    exon_present = []
    support = 1
    coverage = 1

    def __init__(self, coords, present, support = 1, coverage = 1):
        self.block_coordinates = coords
        self.exon_present = present
        self.support = support
        self.coverage = coverage

    def FillGaps(self):
        i = 0
        while i < len(self.exon_present):
            if self.exon_present[i] == 0:
                i += 1
                continue

            j = i + 1
            while j < len(self.exon_present) and self.exon_present[j] == 0:
                j += 1

            if j == len(self.exon_present):
                break

            k = i + 1
            while k < j:
                self.exon_present[k] = -1
                k += 1
            i = j

    def Merge(self, other):
        if len(self.exon_present) != len(other.exon_present):
            return False

        new_exons = [0 for i in range(len(self.exon_present))]

        for i in range(len(self.exon_present)):
            if self.exon_present[i] == other.exon_present[i]:
                new_exons[i] = self.exon_present[i]
            elif self.exon_present[i] == 0:
                new_exons[i] = other.exon_present[i]
            elif other.exon_present[i] == 0:
                new_exons[i] = self.exon_present[i]
            else:
                return False

        self.block_coordinates = [] #not supported yet
        self.exon_present = new_exons
        self.support += other.support
        return True

    def Complete(self):
        has_zero = False
        for i in range(1, len(self.exon_present) - 1):
            if self.exon_present[i] == 0:
                has_zero = True
                break
        return not has_zero

    def InformativeExons(self):
        count = 0
        for i in range(1, len(self.exon_present) - 1):
            if self.exon_present[i] == 0:
                count += 1
        return count

    def Supported(self):
        return self.coverage >= MIN_READS

class BarcodeRecord:
    barcode = ""
    read_infos = []

    def __init__(self, bc, reads = []):
        self.barcode = bc
        self.read_infos = reads

    def AddAlignment(self, read_record):
        self.read_infos.append(read_record)


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def left_of(range1, range2):
    return range1[1] < range2[0]

def equal_ranges(range1, range2):
    return abs(range1[0] - range2[0]) <= 2 and abs(range1[1] - range2[1]) <= 2 

def covers_end(bigger_range, smaller_range):
    return bigger_range[1] <= smaller_range[1] and bigger_range[0] <= smaller_range[0]

def covers_start(bigger_range, smaller_range):
    return bigger_range[0] >= smaller_range[0] and bigger_range[1] >= smaller_range[1]

def contains(bigger_range, smaller_range):
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]


def contains_approx(bigger_range, smaller_range):
    return bigger_range[1] + 2 >= smaller_range[1] and bigger_range[0] - 2 <= smaller_range[0]

def overlaps_to_left(bigger_range, smaller_range):
    return smaller_range[1] >= bigger_range[0] and smaller_range[1] <= bigger_range[1]

def overlaps_to_right(bigger_range, smaller_range):
    return smaller_range[0] >= bigger_range[0] and smaller_range[0] <= bigger_range[1]

def hamming(l1, l2):
    if len(l1) != len(l2):
        return -1
    d = 0
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            d += 1

    return d


def diff_only_present(l1, l2):
    if len(l1) != len(l2):
        return -1
    d = 0
    for i in range(len(l1)):
        if l1[i] == 0 or l2[i] == 0:
            continue
        if l1[i] != l2[i]:
            d += 1
    return d


def count_diff(cov_profile, sign_profile):
    if len(cov_profile) != len(sign_profile):
        return -1
    d = 0
    for i in range(len(sign_profile)):
        if cov_profile[i] == 0 or sign_profile[i] == 0:
            continue
        if cov_profile[i] * sign_profile[i] < 0:
            d += abs(cov_profile[i])
    return d


def sign(i):
    return 0 if i == 0 else (-1 if i < 0 else 1)

class BarcodeStats:
    contig_alignment_count = 0
    contig_count = 0
    distinct_barcodes_in_contigs = 0
    total_matrix_barcodes = 0

    correct = []
    incorrect = []
    incorrect1 = []
    matrix_only = []
    contigs_only = []
    contigs_only_complete = []
    incomplete = []
    incomplete5 = []
    incomplete3 = []
    incomplete_middle = []
    incomplete_in_dict = []
    incomplete_correct = []
    incomplete_incorrect = []
    incomplete_missing = []
    low_covered = []
    ambiguous = []
    with_secondary = []
    missed = []
    used = set()

    def __init__(self):
        self.contig_alignment_count = 0
        self.contig_count = 0
        self.distinct_barcodes_in_contigs = 0
        self.total_matrix_barcodes = 0
        self.correct = []
        self.incorrect = []
        self.incorrect1 = []
        self.matrix_only = []
        self.contigs_only = []
        self.contigs_only_complete = []
        self.incomplete = []
        self.incomplete5 = []
        self.incomplete3 = []
        self.incomplete_middle = []
        self.incomplete_in_dict = []
        self.incomplete_correct = []
        self.incomplete_incorrect = []
        self.incomplete_missing = []
        self.low_covered = []
        self.ambiguous = []
        self.with_secondary = []
        self.missed = []
        self.used = set()

