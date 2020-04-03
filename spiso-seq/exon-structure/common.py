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

def equal_ranges(range1, range2, delta = 1):
    return abs(range1[0] - range2[0]) <= delta and abs(range1[1] - range2[1]) <= delta

def equal_right_border(range1, range2, delta = 1):
    return abs(range1[1] - range2[1]) <= delta

def equal_left_border(range1, range2, delta = 1):
    return abs(range1[0] - range2[0]) <= delta

def covers_end(bigger_range, smaller_range):
    return bigger_range[1] <= smaller_range[1] and bigger_range[0] <= smaller_range[0]

def covers_start(bigger_range, smaller_range):
    return bigger_range[0] >= smaller_range[0] and bigger_range[1] >= smaller_range[1]

def contains(bigger_range, smaller_range):
    return bigger_range[1] >= smaller_range[1] and bigger_range[0] <= smaller_range[0]

def contains_approx(bigger_range, smaller_range, delta = 1):
    return bigger_range[1] + delta >= smaller_range[1] and bigger_range[0] - delta <= smaller_range[0]

def overlaps_to_left(bigger_range, smaller_range):
    return smaller_range[1] >= bigger_range[0] and smaller_range[1] <= bigger_range[1]

def overlaps_to_right(bigger_range, smaller_range):
    return smaller_range[0] >= bigger_range[0] and smaller_range[0] <= bigger_range[1]

# list of of non-overlapping blocks
def total_nonoverlaped_blocks_length(blocks):
    total_len = 0
    for block in sorted(blocks):
        total_len += block[1] - block[0] + 1
    return total_len

# list of of non-overlapping blocks, list of block coverage of which is interesting (i.e. exons)
def total_covering_length(read_blocks, region_of_interested):
    read_pos = 0
    ref_pos = 0
    total_len = 0

    while read_pos < len(read_blocks) and ref_pos < len(region_of_interested):
        while read_pos < len(read_blocks) and left_of(read_blocks[read_pos], region_of_interested[ref_pos]):
            read_pos += 1
        if read_pos == len(read_blocks):
            break

        while ref_pos < len(region_of_interested) and left_of(region_of_interested[ref_pos], read_blocks[read_pos]):
            ref_pos += 1
        if ref_pos == len(region_of_interested):
            break

        if overlaps(region_of_interested[ref_pos], read_blocks[read_pos]):
            total_len += min(region_of_interested[ref_pos][1], read_blocks[read_pos][1]) - max(region_of_interested[ref_pos][0], read_blocks[read_pos][0]) + 1

            if (read_blocks[read_pos][1] < region_of_interested[ref_pos][1]):
                read_pos += 1
            else:
                ref_pos += 1
        elif left_of(region_of_interested[ref_pos], read_blocks[read_pos]):
            ref_pos += 1
        else:
            read_pos += 1

    return total_len


def junctions_from_blocks(blocks):
    junctions = []
    if len(blocks) >= 2:
        for i in range(0, len(blocks) - 1):
            if blocks[i][1] + 1 < blocks[i + 1][0]:
                junctions.append((blocks[i][1] + 1, blocks[i + 1][0] - 1))
    return junctions


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

def intersection(set1, set2):
    res = set()
    for el in set1:
        if el in set2:
            res.add(el)
    return  res


def is_subprofile(short_isoform_profile, long_isoform_profile):
    if len(short_isoform_profile) != len(long_isoform_profile):
        return None

    if all(el == -1 for el in short_isoform_profile):
        return None

    short_range_start = None
    short_range_end = None
    for i in range(len(short_isoform_profile)):
        if short_isoform_profile[i] != 0:
            if short_range_start is None:
                short_range_start = i
            short_range_end = i

    if short_range_start is None or short_range_end is None:
        return False
    for i in range(short_range_start, short_range_end + 1):
        if short_isoform_profile[i] != long_isoform_profile[i]:
            return False
    return True


def concat_gapless_blocks(blocks, cigar_list):
    cigar_index = 0
    block_index = 0
    resulting_blocks = []

    current_block = None
    deletions_before_block = 0

    while cigar_index < len(cigar_list) and block_index < len(blocks):
        #init new block
        if current_block is None:
            #init new block from match
            if cigar_list[cigar_index][0] == 0:
                current_block = (blocks[block_index][0] - deletions_before_block, blocks[block_index][1])
                deletions_before_block = 0
                block_index += 1
            # keep track of deletions before matched block
            elif cigar_list[cigar_index][0] == 2:
                deletions_before_block = cigar_list[cigar_index][1]
        # found intron, add current block
        elif cigar_list[cigar_index][0] == 3:
            resulting_blocks.append(current_block)
            current_block = None
        # add deletion to block
        elif cigar_list[cigar_index][0] == 2:
            current_block = (current_block[0], current_block[1] + cigar_list[cigar_index][1])
        # found match - merge blocks
        elif cigar_list[cigar_index][0] == 0:
            current_block = (current_block[0], blocks[block_index][1])

            block_index += 1
        cigar_index += 1

    if current_block is not None:
        resulting_blocks.append(current_block)

    return resulting_blocks


def normalize_alignment_blocks(alignment):
    aligned_blocks = concat_gapless_blocks(alignment.get_blocks(), alignment.cigartuples)
    return list(map(lambda x: (x[0] + 1, x[1]), aligned_blocks))


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

