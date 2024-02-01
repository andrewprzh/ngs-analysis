from ctypes import cdll
lib = cdll.LoadLibrary('./cseqlib.so')
import ctypes
from collections import defaultdict


class KmerIndexer:
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_strings, kmer_size=6):
        self.seq_list = list(known_strings)
        self.k = kmer_size
        self.index = defaultdict(list)
        self._index()

    def _get_kmers(self, seq):
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    def _index(self):
        for i, barcode in enumerate(self.seq_list):
            for kmer in self._get_kmers(barcode):
                self.index[kmer].append(i)

    def append(self, barcode):
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer in self._get_kmers(barcode):
            self.index[kmer].append(index)

    def empty(self):
        return len(self.seq_list) == 0

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)
        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index[kmer]:
                barcode_counts[i] += 1
                barcode_positions[i].append(pos)

        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
                continue
            if ignore_equal and self.seq_list[i] == sequence:
                continue
            result.append((self.seq_list[i], count, barcode_positions[i]))

        if not result:
            return {}

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        if max_hits == 0:
            return {x[0]:x for x in result}
        return {x[0]:x for x in list(result)[:max_hits]}


BARCODES = ["CCGAATGTGTGGTT", "CGAATGTTAAGGTG", "CCTAAACACTTTAT", "CTTTCTACCAAGCA", "CTATAACGGGCTCG", "CGACACGGTATAGC", "CATTAAAGACTATC", "CGGAAACTAATGCG", "CTGAACCCGGGGAG", "CACCAAATTCTGAA"]

TEST_SEQ = [
"GTGTACTTCGTTCAGTTACCAATTTGGGTGTTTAGCATGGTCATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGCCGAATGTGTGGTTTTGGTGCTGATATTGTGGGGG0GAAAATGTCCTCGGCATAAAAGCGCCATTTTAATTTAAGAAAACGGGGAACTATGC", 
"TCAGGATTGGATTTATATGACTGATCAGTTTCCTCTGCTGTTATCGAAAGCAGATATCAAATGGCTGTGGACGAATGTTAAGGTGGGAATGCAGGTGATTGGAGTTGGTCCAAAGGAAGTTGTGAGTTCTGGGAGAGGCAGAAGGAAAGCAGCTGCCATGTTCTGAAGGTTATCAGCACCTG", 
"TTGGTGATAGAAACAGGCACAGAAGGTGTTAGCAGGTTCCTGTTTGTCCTCTCGCACCCCCTCCCTCCTGATGGTGACCTTGTCCCAGGTCTTCTACCAGGCCTTCTACCAGATGCCTTGCTGAGAATTCACAGAGGCCGTAGACCTGAAGAACCAACAACCTTCCAT", 
"AGGCTTTGAGGTCCCACTCCGGCAGCAGGACGTGCTGGCTCCCAGGAATCACTGACATCAAGGCGTGTAAAATAACACAAAGAGGTTTGCAAAGTCCACAGCCCAAGAGGACCTAAACACTTATGCCGGAGCGTCCTGTTTTTATTTTTTTTTTTTTTCCATCACGCTTCGTTGTTCACGTTGCTTGTGTGTGTA", 
"AAACAGGGTACAATATTTAGAAACGGTGAAAGGAAATGACTGGCCTAAAACTCTTGTGATTCAGTGACTCAAGGATGATTGACACTGTGTAAAAACAGGCACATTAGACCAAGAGATAATTTGCCTTTTTTTTTGTTTTTAAAAATTAAACCTATGGAC", 
"ATGAATAAAATGAAAAACTTTATAAGCACAGCTTATTGAACTAAAGGGCCTTAAGACCATCACTTCCAACGGCCTCATTTTACAGGAGAGAAATCTGAGGACCAGGGAAACCAAGTGATTTTTTTTCTTTTCTGGGGTCACATGGAATGTCAACAGCAGAGCTGTGAAGGCGTTCA", 
"GGTCTGCTAACCTCCCGGCTATGCTCATTCATGGAGAGTGCTTCGAAGAGTGTTTGCAACATTTAGTCACAGTTTATCTTTGGTTCCAATTCCACTTTATTTTTATTTTTAATGTGCATTAAAGACTATTGTGAAAATGGCCCAGATTCATATGATTTGTTGCAGGTCAAACAGGTATTAG", 
"AAAAACTGCCAAGCTTGCACCGCTTATGTATAGTTATTTGTTGTGTATGTGCAAGTGTTTGTATGTGTGTGAGCACATAAGCATAATCTCTTTTTTTTTTTTACACACACACACACACACCATTCCTACACTGGACCCGGGGAGTCAAAAAGCTCTGAAAATTAAACTTTTTCATAAATTTGTGACAAATT"]

'''
kmer_index = lib.KmerIndexer_new(7)
for b in BARCODES:
    lib.KmerIndexer_add(kmer_index, b.encode(),  ctypes.c_int(len(b)))

s = TEST_SEQ[0]
ArrayType = ctypes.c_int * 10000
result = ArrayType()
lib.KmerIndexer_get_occurrences(kmer_index, s.encode(), ctypes.c_int(len(s)), ctypes.c_int(0), ctypes.c_int(1), ctypes.c_int(1), result, ctypes.c_int(100))

print(result[0])
print(BARCODES[result[1]])
pos_count = result[2]
print(pos_count)
for i in range(pos_count):
    print(result[3 + i])



py_kmer_index = KmerIndexer(BARCODES, 7)
print(py_kmer_index.get_occurrences(s))
'''

#py_kmer_index = KmerIndexer(BARCODES, 7)
kmer_index = lib.KmerIndexer_new(7)
for b in BARCODES:
    lib.KmerIndexer_add(kmer_index, b.encode(),  ctypes.c_int(len(b)))
    
for i in range(500000):
    s = TEST_SEQ[i % 8]
    ArrayType = ctypes.c_int * 10000
    result = ArrayType()
    lib.KmerIndexer_get_occurrences(kmer_index, s.encode(), ctypes.c_int(len(s)), ctypes.c_int(0), ctypes.c_int(1), ctypes.c_int(1), result, ctypes.c_int(100))
    #py_kmer_index.get_occurrences(s)
    #kmer_index = lib.KmerIndexer_new(7)
    #for b in BARCODES:
    #    lib.KmerIndexer_add(kmer_index, b.encode(),  ctypes.c_int(len(b)))

    

