#include<string>
#include<vector>
#include<iostream>



char rev_comp_arr[256] = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T', ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', 'N', ' ', ' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 't', ' ', 'g', ' ', ' ', ' ', 'c', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'a', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

extern "C" {

void rev_comp(const wchar_t* s, wchar_t* rc, int len) {
    for (int i = len - 1; i >= 0; --i) {
        int ri =  len - 1 - i;
        rc[ri] = rev_comp_arr[s[i]];
        /// std::cout << char(s[i]) << " " << char(rc[ri]) << std::endl;
    }
}

}

/*

int main() {
    std::string test[8] = {
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
    }
    return 0;
}


*/
