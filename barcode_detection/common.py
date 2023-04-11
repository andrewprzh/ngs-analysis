def find_polya(self, seq):
    if len(seq) < self.window_size:
        return -1
    i = 0
    a_count = seq[0:self.window_size].count('A')
    while i < len(seq) - self.window_size:
        if a_count >= self.polyA_count:
            break
        first_base_a = seq[i] == 'A'
        new_base_a = i + self.window_size < len(seq) and seq[i + self.window_size] == 'A'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        i += 1

    if i >= len(seq) - self.window_size:
        return -1

    return i + max(0, seq[i:].find('AA'))


base_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', " ": " "}

def reverese_complement(my_seq):  ## obtain reverse complement of a sequence
    lms = list(map(lambda x: base_comp[x], my_seq))[::-1]
    return ''.join(lms)