#include<string>
#include<cstring>
#include<vector>
#include<unordered_map>
#include<iostream>
#include<algorithm>


extern "C" {

void rev_comp(const char* s, char* rc, int len);

int find_polyt_start(const char* seq, int seq_end, int seq_start, int window_size6, int min_t_count);

}


char rev_comp_arr[256] = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T', ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', 'N', ' ', ' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 't', ' ', 'g', ' ', ' ', ' ', 'c', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'a', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};


void rev_comp(const char* s, char* rc, int len) {
    for (int i = len - 1; i >= 0; --i) {
        int ri =  len - 1 - i;
        rc[ri] = rev_comp_arr[s[i]];
    }
}


int find_polyt_start(const char* seq, int seq_end, int seq_start, int window_size, int min_t_count) {
    if (seq_end - seq_start < window_size) { return -1; }
    
    int curr_t_count = 0;
    for (int j = seq_start; j < seq_start + window_size; ++j) {
        if (seq[j] == 'T') { ++curr_t_count; }
    }
    
    int idx = seq_start;
    int max_idx = seq_end - window_size;
    while (curr_t_count < min_t_count && idx < max_idx) {
        if (seq[idx] == 'T') { --curr_t_count; }
        if (seq[idx + window_size] == 'T') { ++curr_t_count; }
        ++idx;
    }
    
    if (curr_t_count < min_t_count) { return -1; }
    
    int small_window_size = 3;
    int min_polyt_stat_count = 3;
    curr_t_count = 0;
    int i = idx;
    for (; i < idx + small_window_size; ++i) {
        if (seq[i] == 'T') {
            ++curr_t_count;
        }
    }
    
    while (curr_t_count < min_polyt_stat_count && i <= idx + window_size - small_window_size) {
        if (seq[i] == 'T') { --curr_t_count; }
        if (seq[i + small_window_size] == 'T') { ++curr_t_count; }
        ++i;
    }
    if (curr_t_count >= min_polyt_stat_count) { return i; }

    return -1;
}


class KmerIndexer {
    int kmer_size_;
    int seq_counter_;
    std::unordered_map<std::string, std::vector<int>> kmer_index_;

public:
    KmerIndexer(int kmer_size): kmer_size_(kmer_size), seq_counter_(0), kmer_index_() {}
    
    void add(const char* seq, int len) {
        for (int i = 0; i <= len - kmer_size_; ++i) {
            auto kmer = std::string(seq + i, kmer_size_);
            auto kmer_it = kmer_index_.find(kmer);
            if (kmer_it == kmer_index_.end()) {
                kmer_it = kmer_index_.emplace(std::make_pair(kmer, std::vector<int>())).first;
            }
            kmer_it->second.push_back(seq_counter_);
        }
        ++seq_counter_;
    }

    bool empty() const {
        return seq_counter_ == 0;
    }
    
    void get_occurrences(const char* seq, int len, int max_hits, int min_kmers, int hits_delta, int* result, int max_res_size) const {
        std::unordered_map<int, int> kmer_counts;
        std::unordered_map<int, std::vector<int>> kmer_positions;
        int max_count = 0;
        
        for (int i = 0; i <= len - kmer_size_; ++i) {
            auto kmer = std::string(seq + i, kmer_size_);
            auto kmer_it = kmer_index_.find(kmer);
            if (kmer_it == kmer_index_.end()) {
                continue;
            }
            
            for (const auto& idx : kmer_it->second) {
                kmer_counts[idx] += 1;
                kmer_positions[idx].push_back(i);
                if (max_count < kmer_counts[idx]) {
                    max_count = kmer_counts[idx];
                }
            }
        }
        
        int count_cutoff = std::max(min_kmers, max_count - hits_delta);
        std::vector<int> result_seqs;
        result_seqs.reserve(kmer_counts.size());
        for (const auto& kmer_count : kmer_counts) {
            if (kmer_count.second >= count_cutoff) {
                result_seqs.push_back(kmer_count.first);
            }
        }
        auto& kmer_counts_ = kmer_counts;
        std::sort(result_seqs.begin(), result_seqs.end(), [&kmer_counts_](int seq_idx1, int seq_idx2) { return kmer_counts_[seq_idx1] < kmer_counts_[seq_idx2]; } );
        

        int top_hits = std::min((int) result_seqs.size(), max_hits);
        if (max_hits == 0) {
            top_hits = (int) result_seqs.size();
        }
        
        result[0] = top_hits;
        int current_result_index = 1;
        for (int i = 0; i < top_hits; ++i) {
            result[current_result_index++] = result_seqs[i];
            const auto& positions = kmer_positions[result_seqs[i]];
            result[current_result_index++] = int(positions.size());

            std::memcpy(result + current_result_index, positions.data(), positions.size() * sizeof(int));
            current_result_index += positions.size();
            if (current_result_index >= max_res_size) {
                break;
            }
        }
        
    }

};


extern "C" {

KmerIndexer* KmerIndexer_new(int kmer_size){ return new KmerIndexer(kmer_size); }

void KmerIndexer_add(KmerIndexer* kmer_indexer, const char* seq, int len) { kmer_indexer->add(seq, len); }

bool KmerIndexer_empty(const KmerIndexer* kmer_indexer) { return kmer_indexer->empty(); }

void KmerIndexer_get_occurrences(const KmerIndexer* kmer_indexer, const char* seq, int len, int max_hits, int min_kmers, int hits_delta, int* result, int max_res_size) { 
    kmer_indexer->get_occurrences(seq, len, max_hits, min_kmers, hits_delta, result, max_res_size); 
}

}



int main() {
    /* std::string test[8] = {
"GTGTACTTCGTTCAGTTACCAATTTGGGTGTTTAGCATGGTCATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGTGGGGGTTTTGAAAATGTCCTCGGCATAAAAGCGCCATTTTAATTTAAGAAAACGGGGAACTATGC", 
"TCAGGATTGGATTTATATGACTGATCAGTTTCCTCTGCTGTTATCGAAAGCAGATATCAAATGGCTGTGGAGGAATGCAGGTGATTGGAGTTGGTCCAAAGGAAGTTGTGAGTTCTGGGAGAGGCAGAAGGAAAGCAGCTGCCATGTTCTGAAGGTTATCAGCACCTG", 
"TTGGTGATAGAAACAGGCACAGAAGGTGTTAGCAGGTTCCTGTTTGTCCTCTCGCACCCCCTCCCTCCTGATGGTGACCTTGTCCCAGGTCTTCTACCAGGCCTTCTACCAGATGCCTTGCTGAGAATTCACAGAGGCCGTAGACCTGAAGAACCAACAACCTTCCAT", 
"AGGCTTTGAGGTCCCACTCCGGCAGCAGGACGTGCTGGCTCCCAGGAATCACTGACATCAAGGCGTGTAAAATAACACAAAGAGGTTTGCAAAGTCCACAGCCCAAGAGGAGCCGGAGCGTCCTGTTTTTACCATCACGCTTCGTTGTTCACGTTGCTTGTGTGTGTA", 
"AAACAGGGTACAATATTTAGAAACGGTGAAAGGAAATGACTGGCCTAAAACTCTTGTGATTCAGTGACTCAAGGATGATTGACACTGTGTAAAAACAGGCACATTAGACCAAGAGATAATTTGAAACCTTATTATTGGGTATTGTTTTTAAAAATTAAACCTATGGAC", 
"ATGAATAAAATGAAAAACTTTATAAGCACAGCTTATTGAACTAAAGGGCCTTAAGACCATCACTTCCAACGGCCTCATTTTACAGGAGAGAAATCTGAGGACCAGGGAAACCAAGTGACTTTTCTGGGGTCACATGGAATGTCAACAGCAGAGCTGTGAAGGCGTTCA", 
"GGTCTGCTAACCTCCCGGCTATGCTCATTCATGGAGAGTGCTTCGAAGAGTGTTTGCAACATTTAGTCACAGTTTATCTTTGGTTCCAATTCCACATTTACTCTATTTTTAATGTGTGTGAAAATGGCCCAGATTCATATGATTTGTTGCAGGTCAAACAGGTATTAG", 
"AAAAACTGCCAAGCTTGCACCGCTTATGTATAGTTATTTGTTGTGTATGTGCAAGTGTTTGTATGTGTGTGAGCACATAAGCATAATCTCTTTACACACACACACACACACCATTCCTACATCAAAAAGCTCTGAAAATTAAACTTTTTCATAAATTTGTGACAAATT"};

    int len = test[0].length();
    for (int i = 0; i <= 1000000; ++i) {
        char* rc = rev_comp(test[i % 8].c_str(), len);
        delete rc;
    }*/
    
    std::string s = std::string("AAAAACTGCCAAGCTTGCACCGCTTATGTATAGTTATTTGTTGTGTATGTGCAAGTGTTTGTATGTGTGTGAGCACATAAGCATAATCTCTTTACACACACACACACACACCATTCCTACATCAAAAAGCTCTGAAAATTAAACTTTTTCATAAATTTGTGACAAATT");
    int pos = find_polyt_start(s.c_str(), s.length(), 80, 16, 12);
    std::cout  << pos << std::endl;
    return 0;
}



// GTGTACTTCGTTCAGTTACCAATTTGGGTGTTTAGCATGGTCATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGTGGGGGTTTTGTTTTTTTTTTTTAAAATGTCCTCGGCATAAAAGCGCCATTTTAATTTAAGAAAACGGGGAACTATGC

