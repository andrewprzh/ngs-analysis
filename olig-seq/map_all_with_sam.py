import os
import sys
from Bio import SeqIO
import ssw
import pysam

SCORE_THRESHOLD = 300
MAX_INDELS = 0

MIN_LEN = 169
MAX_LEN = 169


def print_dict(d):
    for k in sorted(d.keys()):
        print(str(k) + '\t' + str(d[k]))

class Stats:
    good_reads = 0
    oligs_covered = 0
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

    def InitOligs(self, oligs):
        for o in oligs:
            self.oligs_coverage[o] = 0
        

    def add(self, olig_id, alignment, read):
        if alignment.score < SCORE_THRESHOLD or alignment.insertion_count + alignment.deletion_count > MAX_INDELS:
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

        self.good_reads += 1


    def print_report(self):
        olig_cov_stat = {}
        for o,c in self.oligs_coverage.iteritems():
            if c > 0:
                self.oligs_covered += 1
            if c not in olig_cov_stat:
                olig_cov_stat[c] = 0
            olig_cov_stat[c] += 1

        print("Too short: " + str(self.too_short) + ", too long " + str(self.too_long) + ", low score " + str(self.low_score))
        print("Mismatch length distribution")
        print_dict(self.mismatches)
        print("Deletion length distribution")
        print_dict(self.deletions)
        print("Insertion length distribution")
        print_dict(self.insertions)
        print("Olig coverage distribution")
        print_dict(olig_cov_stat)
        print("Read length distribution")
        print_dict(self.lengths)
        print("Per olig coverage")
        print_dict(self.oligs_coverage)


def map_from_sam(samfile_in):
    read_mappings = {}
    for r in samfile_in:
        read_id = r.query_name
        olig_id = None if r.reference_id < 0 else r.reference_name
        if read_id is None or olig_id is None:
            continue
        if read_id not in read_mappings:
            read_mappings[read_id] = []
        read_mappings[read_id].append(olig_id)

    return read_mappings


def map_reads(oligs, reads, read_mappings):
    sys.stderr.write("   Mapping with ssw\n")
    aligner = ssw.Aligner()
    stats = Stats()
    stats.InitOligs(oligs)

    i = 0
    prev_percent = 0
    for r in read_mappings.keys():
        i += 1
        percent = 100 * i / len(read_mappings)
        if percent - 5 >= prev_percent:
            sys.stderr.write("   " + str(percent) + "% done\r")
            prev_percent = percent

        read = str(reads[r].seq)
        if len(read) < MIN_LEN:
            stats.add_short()
            continue
        if len(read) > MAX_LEN:
            stats.add_long()
            continue

        max_score = 0
        best_al = None
        best_seqs = None
        for o in read_mappings[r]:
            alignment = aligner.align(reference=str(oligs[o].seq), query=read)
            score = alignment.score
            if max_score < score:
                max_score = score
                best_al = alignment
                best_seqs = (r, o)

        stats.add(best_seqs[1], best_al, read)
        #print(best_seqs[0] + " -> " + best_seqs[1])
        #print(best_al.alignment_report())

    return stats


if len(sys.argv) < 5:
    print("Usage: " + sys.argv[0] + " <oligs in FASTA> <processed reads in FASTA> <SAM file> <stats file> > <output.txt>")
    exit(0)

oligs = SeqIO.index(sys.argv[1], "fasta")
reads = SeqIO.index(sys.argv[2], "fasta")
samfile_in = pysam.AlignmentFile(sys.argv[3], "rb")

read_mappings = map_from_sam(samfile_in)
stats = map_reads(oligs, reads, read_mappings)
stats.print_report()

statf = open(sys.argv[4], "a+")
statf.write(str(stats.good_reads) + "\t" + str(stats.oligs_covered) + "\n")
statf.close()





