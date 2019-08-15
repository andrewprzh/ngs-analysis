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

def equal_ranges(range1, range2, delta = 3):
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta 

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


def find_matching_positions(l1, l2):
    if len(l1) != len(l2):
        return -1
    matches = [0 for i in range(len(l1))]
    for i in range(len(l1)):
        if l1[i] == l2[i]:
            matches[i] = 1
    return matches


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


def is_subprofile(short_isoform_profile, long_isoform_profile):
    if len(short_isoform_profile) != len(long_isoform_profile):
        return None

    if all(el == -1 for el in short_isoform_profile):
        return None

    short_range_start = None
    short_range_end = None
    for i in range(len(short_isoform_profile)):
        if short_isoform_profile[i] == 1:
            if short_range_start is None:
                short_range_start = i
            short_range_end = i

    for i in range(short_range_start, short_range_end + 1):
        if short_isoform_profile[i] != long_isoform_profile[i]:
            return False
    return True


def sign(i):
    return 0 if i == 0 else (-1 if i < 0 else 1)


def table_to_str(d, write_coordinates = False, delim = '\t'):
        vertical_keys = set()
        horisontal_keys = set()
        for k in sorted(d.keys()):
            vertical_keys.add(k[0])
            horisontal_keys.add(k[1])

        vertical_keys = sorted(list(vertical_keys))
        horisontal_keys = sorted(list(horisontal_keys))

        header = ""
        if write_coordinates:
            header = delim + delim.join(map(str, horisontal_keys)) + "\n"
        res = header
        for x1 in vertical_keys:
            row_els = []
            if write_coordinates:
                row_els.append(str(x1))

            for x2 in horisontal_keys:
                v = d.get((x1, x2), 0)
                row_els.append(str(v))
            res += delim.join(row_els)  + "\n"
        return res


#check whether genes overlap and should be processed together
def genes_overlap(gene_db1, gene_db2):
    if (gene_db1.seqid != gene_db2.seqid):
        print("Processing chromosome " + gene_db2.seqid)
        return False
    return overlaps((gene_db1.start, gene_db1.end), (gene_db2.start, gene_db2.end))

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

