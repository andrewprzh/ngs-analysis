from collections import defaultdict


class KmerIndexer:
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_strings, kmer_size=6):
        self.barcode_list = list(known_strings)
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
        for i, barcode in enumerate(self.barcode_list):
            for kmer in self._get_kmers(barcode):
                self.index[kmer].append(i)

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1):
        barcode_counts = defaultdict(int)
        for kmer in self._get_kmers(sequence):
            for i in self.index[kmer]:
                barcode_counts[i] += 1

        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
                continue
            result.append((self.barcode_list[i], count))

        result = sorted(result, reverse=True, key=lambda x: x[1])
        if max_hits == 0:
            return result
        return result[:max_hits]

