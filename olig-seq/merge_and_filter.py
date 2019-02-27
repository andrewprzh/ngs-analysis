import os
import sys
import ssw

FULL_F_ADAPTERS = ['CCCAGCCGGCCATGGCC', 'CCGCCTCCACCACTGCTAGC']
FULL_R_ADAPTERS = ['GCTAGCAGTGGTGGAGGCGG', 'GGCCATGGCCGGCTGGG']
RESTRICTED_F_ADAPTERS = ['CATGGCC', 'CTAGC']
RESTRICTED_R_ADAPTERS = ['GCTAG', 'GGCCATG']
                                    
F_ADAPTERS = ['CCCAGCCGGCCATGGCC', 'CCGCCTCCACCACTGCTAGC']
R_ADAPTERS = ['GCTAGCAGTGGTGGAGGCGG', 'GGCCATGGCCGGCTGGG']
                                       
GOOD_ADAPTER_THRESHOLD = 0.7
RESTRICTED_GOOD_ADAPTER_THRESHOLD = 1

LENGTH_THRESHOLD = 0
UNEQUAL_BASES_THRESHOLD = 0.1

OLIG_LEN = 132

def print_dict(d):
    for k in sorted(d.keys()):
        print(str(k) + '\t' + str(d[k]))


def comp(letter):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[letter.upper()]


def rev_comp(seq):
    return ''.join(map(comp, seq[::-1]))


def hamming(l1, l2):
    if len(l1) != len(l2):
        return -1
    d = 0
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            d += 1

    return d


class Stats:
    total_reads = 0
    with_adapters = 0
    lengths = {}
    non_equal_bases = {}
    only_first_read = 0
    only_second_read = 0
    both = 0
    dropped = 0
    too_much_unequal_bases = 0
    short = 0
    clipped_left = {}
    clipped_right = {}
    mm_in_l_adapters = {}
    mm_in_r_adapters = {}
    

    def add(self, length, non_equal, which_read):
        self.with_adapters += 1
        if length not in self.lengths:
            self.lengths[length] = 0
        self.lengths[length] += 1
        if non_equal not in self.non_equal_bases:
            self.non_equal_bases[non_equal] = 0
        self.non_equal_bases[non_equal] += 1
        if which_read == 'b':
            self.both += 1
        elif which_read == 'l':
            self.only_first_read += 1
        else:
            self.only_second_read += 1

    def print_stats(self):
        print("Both reads used: " + str(self.both) + ", only left " + str(self.only_first_read) + ", only right " + str(self.only_second_read))
        print("Dropped: no adapters - " + str(self.dropped) + ", too short " + str(self.short) + ", too discordant " + str(self.too_much_unequal_bases))
        print("Distribution of bases clipped in left adapters")
        print_dict(self.clipped_left)
        print("Distribution of bases clipped in right adapters")
        print_dict(self.clipped_right)
        print("Distribution of mismatches in left adapters")
        print_dict(self.mm_in_l_adapters)
        print("Distribution of mismatches in right adapters")
        print_dict(self.mm_in_r_adapters)
        print("Length distribution")
        print_dict(self.lengths)
        print("Number of unequal bases distribution")
        print_dict(self.non_equal_bases)

class Read:
    read_id = ''
    seq = ''
    qual = ''

    def __init__(self, read_id = '', seq = '', qual = ''):
        self.read_id = read_id
        self.seq = seq
        self.qual = qual

    def read_from_file(self, f):
        self.read_id = f.readline()
        if not self.read_id:
            return False
        self.seq = f.readline()
        if not self.seq:
            return False
        id2 = f.readline()
        if not id2:
            return False
        self.qual = f.readline()
        if not self.qual:
            return False
        self.read_id = self.read_id.strip()
        self.seq = self.seq.strip()
        self.qual = self.qual.strip()
        return True

    def trim(self):
        for a in R_ADAPTERS:
            rpos = self.seq.find(a)
            if rpos != -1:
                self.seq = self.seq[:rpos]
                self.qual = self.qual[:rpos]
                break
        for a in F_ADAPTERS:
            fpos = self.seq.find(a)
            if fpos != -1:
                self.seq = self.seq[fpos + len(a):]
                self.qual = self.qual[fpos + len(a):]
                break

    def length(self):
        return len(self.seq)

    def rc(self):
        self.seq = rev_comp(self.seq)
        self.qual = self.qual[::-1]
        return self

    def append_adapters(self, left, right):
        self.seq = left + self.seq + right
        self.qual = ''.join(['G' for l in range(len(left))]) + self.qual +  ''.join(['G' for l in range(len(right))]) 

    def write_to_file(self, outf, out_format):
        if self.length() > 0:
            if out_format == "FASTQ":
                outf.write('@' + self.read_id[1:] + '\n')
                outf.write(self.seq + '\n')
                outf.write('+' + self.read_id[1:] + '\n')
                outf.write(self.qual + '\n')
            else:
                outf.write('>' + self.read_id[1:] + '\n')
                outf.write(self.seq + '\n')


def MergeReads(read1, read2, rc):
    if read1.length() > read2.length():
        return read2, 0
    elif read1.length() < read2.length():
        return read1, 0

    r = Read()
    r.read_id = read1.read_id
    r1 = read1
    r2 = read2
    if rc:
        r1.rc()
    else:
        r2.rc()
    non_equal_bases = 0
    for i in range(0, r1.length()):
        if r2.seq[i] != r1.seq[i]:
            non_equal_bases += 1
        if r1.qual[i] >= r2.qual[i]:
            r.seq += r1.seq[i]
            r.qual += r1.qual[i]
        else:
            r.seq += r2.seq[i]
            r.qual += r2.qual[i]
    return r, non_equal_bases


def IsGoodAdapter(adapter, adapters):
    aligner = ssw.Aligner()
    max_score = 0
    adapter_len = 0
    for a in adapters:
        alignment = aligner.align(reference=a, query=adapter)
        if max_score < alignment.score:
            max_score = alignment.score
            adapter_len = len(a)

    return max_score / 2 >= adapter_len * GOOD_ADAPTER_THRESHOLD


def HasGoodMainAdapter(original, trimmed):
    main_f_adapter = F_ADAPTERS[0]
    pos = original.seq.find(trimmed.seq)
    read_adapter = original.seq[max(0, pos - 20):pos]
     
    return IsGoodAdapter(read_adapter, [main_f_adapter])
    

def RigthReadHasFirstAdapter(original1, original2, trimmed1, trimmed2):
    d1 = HasGoodMainAdapter(original1, trimmed1)
    d2 = HasGoodMainAdapter(original2, trimmed2)

    if not d1 and not d2:
        return None
    elif d1:
        return False
    else:
        return True

def ConfirmAdapters(original, trimmed):
    pos_left = original.seq.find(trimmed.seq)
    pos_right = pos_left + len(trimmed.seq)
    left_adapter = original.seq[:pos_left]
    right_adapter = original.seq[pos_right:pos_right + 20]
    
    return IsGoodAdapter(left_adapter, F_ADAPTERS) and IsGoodAdapter(right_adapter, R_ADAPTERS)


def CheckAdapters(original, trimmed, rc, stats):
    l_ref_adapter = F_ADAPTERS[1 if rc else 0]
    r_ref_adapter = R_ADAPTERS[1 if rc else 0]

    pos_left = original.seq.find(trimmed.seq)
    pos_right = pos_left + len(trimmed.seq)
    left_adapter = original.seq[max(0, pos_left - len(l_ref_adapter)) : pos_left]
    right_adapter = original.seq[pos_right : pos_right + len(r_ref_adapter)]

    aligner = ssw.Aligner()
    l_alignment = aligner.align(reference=l_ref_adapter, query=left_adapter)
    r_alignment = aligner.align(reference=r_ref_adapter, query=right_adapter)
    #print(l_alignment.alignment_report())
    #print(r_alignment.alignment_report())

    clipped = len(l_ref_adapter) - (l_alignment.match_count + l_alignment.mismatch_count + l_alignment.insertion_count + l_alignment.deletion_count)
    if clipped not in stats.clipped_left:
        stats.clipped_left[clipped] = 0
    stats.clipped_left[clipped] += 1

    clipped = len(r_ref_adapter) - (r_alignment.match_count + r_alignment.mismatch_count + r_alignment.insertion_count + r_alignment.deletion_count)
    if clipped not in stats.clipped_right:
        stats.clipped_right[clipped] = 0
    stats.clipped_right[clipped] += 1

    mm = l_alignment.mismatch_count
    if mm not in stats.mm_in_l_adapters:
        stats.mm_in_l_adapters[mm] = 0
    stats.mm_in_l_adapters[mm] += 1

    mm = r_alignment.mismatch_count
    if mm not in stats.mm_in_r_adapters:
        stats.mm_in_r_adapters[mm] = 0
    stats.mm_in_r_adapters[mm] += 1



def FindOlig(original1, original2, trimmed1, trimmed2, rc, stats, switched = False):
    ldiff = abs(trimmed1.length() - OLIG_LEN)
    rdiff = abs(trimmed2.length() - OLIG_LEN)
    if ldiff > rdiff:
        return FindOlig(original2, original1, trimmed2, trimmed1, not rc, stats, True)

    r = Read()
    if ConfirmAdapters(original1, trimmed1):
        CheckAdapters(original1, trimmed1, rc, stats)
        r = trimmed1.rc() if rc else trimmed1
        which_read = 'l' if not switched else 'r'
        r.read_id = r.read_id.split()[0]
        if r.length() < LENGTH_THRESHOLD:
            stats.short += 1
            r = Read()
        else:
            stats.add(r.length(), -1, which_read)
    elif ConfirmAdapters(original2, trimmed2):
        CheckAdapters(original2, trimmed2, not rc, stats)
        r = trimmed2.rc() if not rc else trimmed2
        which_read = 'l' if switched else 'r'
        r.read_id = r.read_id.split()[0]
        if r.length() < LENGTH_THRESHOLD:
            stats.short += 1
            r = Read()
        else:
            stats.add(r.length(), -1, which_read)
    else:
        stats.dropped += 1
    return r
        

def ProduceOlig(original1, original2, trimmed1, trimmed2, stats):
    trimmed1.trim()
    trimmed2.trim()
    rc = RigthReadHasFirstAdapter(original1, original2, trimmed1, trimmed2)
    r = Read()
    if rc == None:
        stats.dropped += 1
        return r

    if trimmed1.length() == trimmed2.length():
        CheckAdapters(original1, trimmed1, rc, stats)
        CheckAdapters(original2, trimmed2, not rc, stats)

        r, non_equal_bases = MergeReads(trimmed1, trimmed2, rc)
        if r.length() < LENGTH_THRESHOLD:
            stats.short += 1
            r = Read()
        elif non_equal_bases > UNEQUAL_BASES_THRESHOLD * trimmed1.length():
            stats.too_much_unequal_bases += 1
            r = Read()
        else:
            stats.add(r.length(), non_equal_bases, 'b')
            r.read_id = r.read_id.split()[0]
    else:
        r = FindOlig(original1, original2, trimmed1, trimmed2, rc, stats)
        stats.only_first_read += 1
#        print(r.read_id + " " + r.seq + " " + str(rc) )
#        raw_input('')         



    if r.length() > 0:
        r.append_adapters(FULL_F_ADAPTERS[0], FULL_R_ADAPTERS[0])
    return r


def ProcessAll(orig1, orig2, trim1, trim2, outf, out_format):
    stats = Stats()
    while True:
        original1 = Read()
        original2  = Read()
        trimmed1 = Read() 
        trimmed2  = Read()

        if not original1.read_from_file(orig1) or not original2.read_from_file(orig2)  or not trimmed1.read_from_file(trim1)  or not trimmed2.read_from_file(trim2):
            break

        stats.total_reads += 1
        read = ProduceOlig(original1, original2, trimmed1, trimmed2, stats)
        read.write_to_file(outf, out_format)

    outf.close()
    return stats
        

if len(sys.argv) < 9:
    print("Usage: " + sys.argv[0] + " <original reads left> <original reads right> <trimmed left> <trimmed right> <out file> <FASTA | FASTQ> <full | restricted> <stat file>")
    exit(0)


out_format = sys.argv[6] 
if sys.argv[7].startswith('r'):
    F_ADAPTERS = RESTRICTED_F_ADAPTERS
    R_ADAPTERS = RESTRICTED_R_ADAPTERS
    GOOD_ADAPTER_THRESHOLD = RESTRICTED_GOOD_ADAPTER_THRESHOLD

stats = ProcessAll(open(sys.argv[1]), open(sys.argv[2]), open(sys.argv[3]), open(sys.argv[4]), open(sys.argv[5], 'w'), out_format)
stats.print_stats() 

statf = open(sys.argv[8], "a+")
statf.write(str(stats.total_reads) + "\t" + str(stats.with_adapters) + "\t")
statf.close()

