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
    strand = ""

    def __init__(self, l):
        tokens = l.split()
        self.chromosome = tokens[0]
        self.feature_type = tokens[2]
        self.coords = (int(tokens[3]), int(tokens[4]))
        self.gene_id = ""
        self.transcript_id = ""
        for i in range(8, len(tokens)):
            if tokens[i]  == "gene_id":
                self.gene_id = tokens[i + 1][1:-2]
            elif tokens[i] == "transcript_id":
                self.transcript_id = tokens[i + 1][1:-2]

            if self.transcript_id != "" and self.gene_id != "":
                break
        self.strand = tokens[6]


class Gene:
    chromosome = ""
    gene_id = ""
    strans = "
    features = {}

    def __init__(self, feature = None):
        self.gene_id = feature.gene_id if feature is not None else "N/A"
        self.chromosome = feature.chromosome if feature is not None else "0"
        self.strand = feature.strand if feature is not None else "N/A"
        self.features = {}

    def add_feature(self, feature):
        if feature.gene_id != self.gene_id:
            print("Wrong gene id")
            return
        if feature.transcript_id == "":
            print("Wrong transcript id")
            return
        if feature.transcript_id not in self.features:
            self.features[feature.transcript_id] = []
        self.features[feature.transcript_id].append(feature)


class CodonInfo:
    chomosome = ""
    gene_id = ""
    strand = ""
    start_codons = set()
    stop_codons = set()
    annotated_codon_pairs = set()
    start_codons_within_exons = set()
    stop_codons_within_exons = set()

    def __init__(self, gene):
        self.chromosome = gene.chromosome
        self.gene_id = gene.gene_id
        self.strand = gene.strand
        self.start_codons = set()
        self.stop_codons = set()
        self.annotated_codon_pairs = set()
        self.start_codons_within_exons = set()
        self.stop_codons_within_exons = set()

        all_cds = set()
        for transcript in gene.features.keys():
            cds = []
            start_codon = None
            stop_codon = None
            for feature in gene.features[transcript]:
                if feature.feature_type == "CDS":
                    cds.append(feature)
                elif feature.feature_type == "start_codon":
                    start_codon = feature
                elif feature.feature_type == "stop_codon":
                    stop_codon = feature

            cds = sorted(cds, key=lambda x: x.coords)
            if start_codon is None and len(cds) > 0:
                if cds[0].strand == '+' and cds[-1].strand == '+':
                    start_codon = cds[0]
                    start_codon.coords = (start_codon.coords[0], start_codon.coords[0] + 2)
                elif cds[0].strand == '-' and cds[-1].strand == '-':
                    start_codon = cds[-1]
                    start_codon.coords = (start_codon.coords[1] - 2, start_codon.coords[1])
                else:
                    print("Incompatible strands")

            if stop_codon is None and len(cds) > 0:
                if cds[0].strand == '+' and cds[-1].strand == '+':
                    stop_codon = cds[-1]
                    stop_codon.coords = (stop_codon.coords[1] - 2, stop_codon.coords[1])
                elif cds[0].strand == '-' and cds[-1].strand == '-':
                    stop_codon = cds[0]
                    stop_codon.coords = (stop_codon.coords[0], stop_codon.coords[0] + 2)
                else:
                    print("Incompatible strands")

            if start_codon is not None and stop_codon is not None:
                self.annotated_codon_pairs.add((start_codon.coords, stop_codon.coords))
            if start_codon is not None:
                self.start_codons.add(start_codon.coords)
            if stop_codon is not None:
                self.stop_codons.add(stop_codon.coords)

            for c in cds:
                all_cds.add(c.coords)

        for codon in self.start_codons:
            for c in all_cds:
                if range_in(codon, c):
                    self.start_codons_within_exons.add(codon)
                    #print("Start codon " + str(codon) + " in CDS " + str(c))
                    break

        for codon in self.stop_codons:
            for c in all_cds:
                if range_in(codon, c):
                    self.stop_codons_within_exons.add(codon)
                    #print("Stop codon " + str(codon) + " in CDS " + str(c))
                    break

    def to_string(self):
        if len(self.start_codons) * len(self.stop_codons) < 1:
            return ""
        s = self.chromosome + '\t' + self.gene_id  + '\t' + self.strand + '\n'
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
    start_codons_in_exons = 0
    stop_codons_in_exons = 0

    
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

            self.start_codons_in_exons += len(codon_info.start_codons_within_exons)
            self.stop_codons_in_exons += len(codon_info.stop_codons_within_exons)


    def print_report(self):
        print("Codon count maxtrx")
        print_table(self.codon_count_matrix)
        print("Histogram of distinct codon pairs count")
        print_dict(self.annotated_pairs_histogram)
        print("Annotated codon pairs vs all possible codon pairs")
        print_table(self.annotated_vs_all_pairs)
        print("Number of start codons within exons")
        print(self.start_codons_in_exons)
        print("Number of stop codons within exons")
        print(self.stop_codons_in_exons)


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








