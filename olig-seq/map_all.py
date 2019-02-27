import os
import sys
from pyfasta import Fasta
import ssw
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


SCORE_THRESHOLD = 300

MIN_LEN = 160
MAX_LEN = 175

class Stats:
    oligs_coverage = {}
    lengths = {}
    deletions = {}
    insertions = {}
    mismatches = {}
    low_score = 0
    too_short = 0
    too_long = 0

    def add_long(self):
        self.too_long += 1

    def add_short(self):
        self.too_short += 1

    def add(self, olig_id, alignment, read):
        if alignment.score < SCORE_THRESHOLD:
            self.low_score += 1
            return
        
        length = len(read)
        if length not in self.lengths:
            self.lengths[length] = 0
        self.lengths[length] += 1

        mm = alignment.mismatch_count
        if mm not in self.mismatches:
            self.mismatches[mm] = 0
        self.mismatches[mm] += 1

        ic = alignment.insertion_count
        if ic not in self.insertions:
            self.insertions[ic] = 0
        self.insertions[ic] += 1

        dc = alignment.deletion_count
        if dc not in self.deletions:
            self.deletions[dc] = 0
        self.deletions[dc] += 1

        if olig_id not in self.oligs_coverage:
            self.oligs_coverage[olig_id] = 0
        self.oligs_coverage[olig_id] += 1

    def print_report(self):
        print("Too short: " + str(self.too_short) + ", too long " + str(self.too_long) + ", low score " + str(self.low_score))
        print("Mismatches: " + str(self.mismatches))
        print("Deletions: " + str(self.deletions))
        print("Insertions: " + str(self.insertions))
        print("Lengths")
        for l in sorted(self.lengths):
            print(str(l) + '\t' + str(self.lengths[l]))
        print("Olig coverage")
        for o in sorted(self.oligs_coverage):
            print(o + '\t' + str(self.oligs_coverage[o]))


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <oligs in FASTA> <processed reads in FASTA> > <output.txt>")
    exit(0)

oligs = Fasta(sys.argv[1])
reads = Fasta(sys.argv[2])

aligner = ssw.Aligner()
stats = Stats()

for r in reads.keys():
    read = str(reads[r][:])
    if len(read) < MIN_LEN:
        stats.add_short()
        continue
    if len(read) > MAX_LEN:
        stats.add_long()
        continue

    max_score = 0
    best_al = None
    best_seqs = None
    for o in sorted(oligs.keys()):
        alignment = aligner.align(reference=str(oligs[o][:]), query=read)
        score = alignment.score
        if max_score < score:
            max_score = score
            best_al = alignment
            best_seqs = (r, o)

    stats.add(best_seqs[1], best_al, read)
    print(best_seqs[0] + " -> " + best_seqs[1])
    print(best_al.alignment_report())

stats.print_report()


