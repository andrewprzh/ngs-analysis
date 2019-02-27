import os
import sys


def range_in(small_range, large_range):
    return small_range[0] > large_range[0] and small_range[1] < large_range[1]

def range_at_start(small_range, large_range):
    return small_range[0] == large_range[0] and small_range[1] <= large_range[1]

def range_at_end(small_range, large_range):
    return small_range[0] >= large_range[0] and small_range[1] == large_range[1]

def print_dict(d):
    for k in sorted(d.keys()):
        print(str(k) + '\t' + str(d[k]))

def print_table(d):
    max_x = 0
    max_y = 0
    for k in sorted(d.keys()):
        if k[0] > max_x:
            max_x = k[0]
        if k[1] > max_y:
            max_y = k[1]

    s = '\t'
    for y in range(max_y + 1):
        s += str(y) + '\t'
    print(s)
    for x in range(max_x + 1):
        s = str(x)
        for y in range(max_y + 1):
            v = d.get((x,y), 0)
            s += '\t' + str(v)
        print(s)
            

class Feature:
    chromosome = ""
    coords = (0, 0)
    gene_id = ""
    transcript_id = ""
    feature_type = ""

    def __init__(self, l):
        tokens = l.split()
        self.chromosome = tokens[0]
        self.feature_type = tokens[2]
        self.coords = (int(tokens[3]), int(tokens[4]))
        self.gene_id = tokens[9][1:-2] if tokens[8] == "gene_id" else ""
        self.transcript_id = tokens[13][1:-2] if tokens[12] == "transcript_id" else ""


class Gene:
    chromosome = ""
    gene_id = ""
    start_codons = []
    stop_codons = []
    annotated_codon_pairs = []
    cds = []

    def __init__(self):
        self.gene_id = "N/A"
        self.chromosome = "0"
        self.start_codons = []
        self.stop_codons = []
        self.annotated_codon_pairs = []
        self.cds = []

    def __init__(self, feature):
        self.gene_id = feature.gene_id
        self.chromosome = feature.chromosome
        self.start_codons = []
        self.stop_codons = []
        self.annotated_codon_pairs = []
        self.cds = []

    def add_feature(self, feature):
        if feature.gene_id != self.gene_id:
            print("Wrong gene id")
            return

        if feature.feature_type == "CDS":
            self.cds.append(feature.coords)
        elif feature.feature_type == "start_codon":
            self.start_codons.append(feature)
        elif feature.feature_type == "stop_codon":
            self.stop_codons.append(feature)
            if len(self.start_codons) > 0 and feature.transcript_id == self.start_codons[-1].transcript_id:
                self.annotated_codon_pairs.append((self.start_codons[-1].coords, feature.coords))
            else:
                pass#print("Cannot find corresponding start codon")


class CodonInfo:
    chomosome = ""
    gene_id = ""
    start_codons = set()
    stop_codons = set()
    annotated_codon_pairs = set()
    start_codons_within_exons = set()

    def __init__(self, gene):
        self.chromosome = gene.chromosome
        self.gene_id = gene.gene_id
        self.start_codons = set()
        self.stop_codons = set()
        self.annotated_codon_pairs = set()
        self.start_codons_within_exons = set()

        for codon in gene.start_codons:
            self.start_codons.add(codon.coords)
        for codon in gene.stop_codons:
            self.stop_codons.add(codon.coords)
        self.annotated_codon_pairs = set(gene.annotated_codon_pairs)

        cds_set = set(gene.cds)
        for codon in self.start_codons:
            for cds in cds_set:
                if range_in(codon, cds):
                    self.start_codons_within_exons.add(codon)
                    break

    def to_string(self):
        if len(self.start_codons) * len(self.stop_codons) <= 1:
            return ""
        s = self.chromosome + '\t' + self.gene_id + '\n'
        for c in sorted(self.start_codons):
            s += str(c[0]) + '\t' + str(c[1]) + '\t' 
        s += '\n'
        for c in sorted(self.stop_codons):
            s += str(c[0]) + '\t' + str(c[1]) + '\t' 
        s += '\n'
        for c in sorted(self.annotated_codon_pairs):
            s += str(c[0][0]) + '\t' + str(c[0][1]) + '\t' + str(c[1][0]) + '\t' + str(c[1][1]) + '\t'
        s += '\n' 
        return s
        

class Stats:
    codon_count_matrix = {}
    annotated_pairs_histogram = {}
    annotated_vs_all_pairs = {}
    codons_in_exons = 0
    
    def add_stats(self, codon_info_list):
        for codon_info in codon_info_list:
            start_codons = len(codon_info.start_codons)
            stop_codons = len(codon_info.stop_codons)
            count_pair = (start_codons, stop_codons)
            self.codon_count_matrix[count_pair] = self.codon_count_matrix.get(count_pair, 0) + 1
            
            annotated_codons = len(codon_info.annotated_codon_pairs)
            self.annotated_pairs_histogram[annotated_codons] = self.annotated_pairs_histogram.get(annotated_codons, 0) + 1

            total_pairs = start_codons * stop_codons
            ann_vs_all = (annotated_codons,total_pairs)
            self.annotated_vs_all_pairs[ann_vs_all] = self.annotated_vs_all_pairs.get(ann_vs_all, 0) + 1

            self.codons_in_exons += len(codon_info.start_codons_within_exons)

    def print_report(self):
        print("Codon count maxtrx")
        print_table(self.codon_count_matrix)
        print("Histogram of distinct codon pairs count")
        print_dict(self.annotated_pairs_histogram)
        print("Annotated codon pairs vs all possible codon pairs")
        print_table(self.annotated_vs_all_pairs)
        print("Number of codons within exons")
        print(self.codons_in_exons)


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <GTF> <output codon list file> > <statistics>")
    exit(0)

gtf = open(sys.argv[1])
out_name = sys.argv[2]
out_f = open(out_name, 'w')

codon_info_list = []
current_gene = None
count = 0
stats = Stats()

for l in gtf:
    count += 1
    if count % 10000 == 0:
        sys.stderr.write("\r   " + str(count) + " lines processed")

    if l.startswith("#"):
        continue
    feature = Feature(l)

    if feature.feature_type == "gene":
        if current_gene is not None:
            codon_info_list.append(CodonInfo(current_gene))
        current_gene = Gene(feature)
    else:
        current_gene.add_feature(feature)

    if len(codon_info_list) > 1000:
        stats.add_stats(codon_info_list)
        for c in codon_info_list:
            out_f.write(c.to_string())
        del codon_info_list[:]

codon_info_list.append(CodonInfo(current_gene))
stats.add_stats(codon_info_list)
for c in codon_info_list:
    out_f.write(c.to_string())
out_f.close()
stats.print_report()








