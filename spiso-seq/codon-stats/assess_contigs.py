import os
import sys
import pysam

DELIM = ','

def to_pairs(arr):
    res = []
    i = 0
    while i < len(arr):
        res.append((arr[i], arr[i+1]))
        i += 2
    return res
        

def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def left_of(range1, range2):
    return range1[1] < range2[0]


def right_of(range1, range2):
    return range1[0] > range2[1]


def contained(range1, range2):
    return range1[1] <= range2[1] and range1[0] >= range2[0]


def dic_to_str(d):
    s = ""
    for k in sorted(d.keys()):
        s += str(k) + '\t' + str(d[k]) + '\n'

def table_to_str(d):
    vertical_keys = set()
    horisontal_keys = set()
    for k in sorted(d.keys()):
        vertical_keys.add(k[0])
        horisontal_keys.add(k[1])

    res = ""
    s = '\t'
    for x2 in horisontal_keys:
        s += str(x2) + '\t'
    res += s + '\n'
    for x1 in vertical_keys:
        s = str(x1)
        for x2 in horisontal_keys:
            v = d.get((x1,x2), 0)
            s += '\t' + str(v)
        res += s + '\n'

    return res


def dict_sum(d):
    s = 0
    for v in d.values():
        s += v
    return s


def table_to_simple_str(d):
    vertical_keys = set()
    horisontal_keys = set()
    for k in sorted(d.keys()):
        vertical_keys.add(k[0])
        horisontal_keys.add(k[1])

    res = ""
    for x1 in vertical_keys:
        s = ""
        for x2 in horisontal_keys:
            v = d.get((x1,x2), 0)
            s += str(v) + '\t'
        res += s + '\n'

    return res


class CodonInfo:
    gene_id = ""
    strand = ""
    start_codons = []
    stop_codons = []
    annotated_codon_pairs = []

    def __init__(self, gene_id, strand, start_codon_line, stop_codon_line, codon_pairs_line):
        self.gene_id = gene_id
        self.strand = strand
        tokens = start_codon_line.strip().split('\t')
        self.start_codons = to_pairs(map(lambda x: int(x), start_codon_line.strip().split('\t')))
        self.stop_codons = to_pairs(map(lambda x: int(x), stop_codon_line.strip().split('\t')))

        tokens = codon_pairs_line.strip().split('\t')
        if len(tokens) > 0 and len(tokens) % 4 == 0:
            self.annotated_codon_pairs = to_pairs(to_pairs(map(lambda x: int(x), tokens)))
        else:
            self.annotated_codon_pairs = []

    def get_gene_coords(self):
        start = self.start_codons[0][0]
        end = self.stop_codons[-1][1]
        return (start, end) if start < end else (end, start)


class ChromosomeWideCodonStructure:
    gene_coords = []
    codon_infos = []
    skipped_genes = 0

    def __init__(self):
        self.gene_coords = []
        self.codon_infos = []
        self.skipped_genes = 0


    def add_gene(self, gene):
        coords = gene.get_gene_coords()

        if len(self.gene_coords) > 0 and not right_of(coords, self.gene_coords[-1]):
            self.skipped_genes += 1
            return False

        self.gene_coords.append(coords)
        if len(self.gene_coords) > 1 and self.gene_coords[-2] > self.gene_coords[-1]:
            sys.stderr.write("Unordered genes " + str(self.gene_coords[-2]) + " " + str(self.gene_coords[-1]))
        self.codon_infos.append(gene)
        return True


    def find_gene_index(self, to_find):
        left = 0
        right = len(self.gene_coords) - 1

        if left_of(to_find, self.gene_coords[left]) or right_of(to_find, self.gene_coords[right]):
            return (-1, -1)

        while True:
            center = (left + right) / 2
            if overlaps(to_find, self.gene_coords[center]):
                return self.detect_bounds(to_find, center)
            elif left >= right:
                #print("Gene not found " + str(to_find))
                #if center > 0:
                #    print(self.gene_coords[center - 1])
                #print(self.gene_coords[center])
                #if center < len(self.gene_coords) - 1:
                #    print(self.gene_coords[center + 1])
                #x = raw_input()
                return (-1, -1)
            elif left_of(to_find, self.gene_coords[center]):
                right = center - 1
            elif right_of(to_find, self.gene_coords[center]):
                left = center + 1
            else:
                sys.stderr.write("Error in binary search\n")

    def get_gene_id(self, index):
        return self.codon_infos[index].gene_id

    def get_gene_coodrs(self, index):
        return self.gene_coords[index]

    def get_strand(self, index):
        return self.codon_infos[index].strand

    def get_annotated_codons(self, index):
        return self.codon_infos[index].annotated_codon_pairs

    def get_all_codons(self, index):
        start_codons = self.codon_infos[index].start_codons
        stop_codons = self.codon_infos[index].stop_codons
        return start_codons, stop_codons

    def detect_bounds(self, to_find, overlap_coord):
        left = overlap_coord    
        right = overlap_coord

        while right < len(self.gene_coords) - 1 and to_find[1] > self.gene_coords[right + 1][0]:
            right += 1
        while left > 0 and to_find[0] < self.gene_coords[left - 1][1]:
            left -= 1

        return (left, right)


    def find_covered_codons(self, index, alignment_blocks):
        start_codons = self.codon_infos[index].start_codons
        stop_codons = self.codon_infos[index].stop_codons

        start_codon_index = 0
        stop_codon_index = 0
        alignment_block_index = 0

        covered_start_codons = []
        covered_stop_codons = []

        while alignment_block_index < len(alignment_blocks):
            while start_codon_index < len(start_codons) and left_of(start_codons[start_codon_index], alignment_blocks[alignment_block_index]):
                start_codon_index += 1
            while stop_codon_index < len(stop_codons) and left_of(stop_codons[stop_codon_index], alignment_blocks[alignment_block_index]):
                stop_codon_index += 1

            if start_codon_index < len(start_codons) and contained(start_codons[start_codon_index], alignment_blocks[alignment_block_index]):
                covered_start_codons.append(start_codons[start_codon_index])
            if stop_codon_index < len(stop_codons) and contained(stop_codons[stop_codon_index], alignment_blocks[alignment_block_index]):
                covered_stop_codons.append(stop_codons[stop_codon_index])

            alignment_block_index += 1

        #if len(covered_start_codons) == 0 and len(covered_stop_codons) == 0:
            #print(alignment_blocks)
            #print(start_codons)
            #print(stop_codons)
            #x = raw_input()
        return covered_start_codons, covered_stop_codons


class GenoWideCodonStructure:
    chr_dict = {}

    def __init__(self, inf_name):
        self.chr_dict = {}
        inf = open(inf_name)

        count = 0
        gene_len = 0
        while True:
            l = inf.readline()
            if not l:
                break

            tokens = l.strip().split('\t')
            if len(tokens) != 3:
                sys.stderr.write("Wrong format\n")
            if tokens[0] not in self.chr_dict:
                self.chr_dict[tokens[0]] = ChromosomeWideCodonStructure()
            
            start_l = inf.readline()
            stop_l = inf.readline()
            pairs_l = inf.readline()
            gene_info = CodonInfo(tokens[1], tokens[2], start_l, stop_l, pairs_l)
            gene_coords = gene_info.get_gene_coords()
            if self.chr_dict[tokens[0]].add_gene(gene_info):
                count += 1
                gene_len += gene_coords[1] - gene_coords[0]

        sys.stderr.write("loaded " + str(count) + " genes of total length " + str(gene_len) + "\n")
        #print(self.chr_dict['17'].gene_coords)
        #print(self.chr_dict['17'].codon_infos[0].gene_id)
        #print(self.chr_dict['17'].codon_infos[0].start_codons)
        #print(self.chr_dict['17'].codon_infos[0].stop_codons)



class ContigCodonStats:
    codon_table = {}
    annotated_codon_table = {}

    def __init__(self, start_codons = [], stop_codons = []):
        self.codon_table = {}
        for start_c in start_codons:
            for stop_c in stop_codons:
                self.codon_table[(start_c, stop_c)] = 0
        self.annotated_codon_table = {}

    def add_contig(self, start_codon, stop_codon, barcode_count = 1):
        if (start_codon, stop_codon) not in self.codon_table:
            self.codon_table[(start_codon, stop_codon)] = 0
        self.codon_table[(start_codon, stop_codon)] += barcode_count

    def add_annotated_codons(self, annotated_codon_pairs):
        for codon_pair in annotated_codon_pairs:
            self.annotated_codon_table[codon_pair] = 1


class Stats:

    cover_multiple_genes = 0
    cover_zero_genes = 0
    cover_zero_codons = 0
    cover_only_start_codon = 0
    cover_only_stop_codon = 0
    cover_multiple_start_codons = 0
    cover_multiple_stop_codons = 0
    cover_multiple_both_codons = 0
    cover_two_codons = 0
    skipped_genes = 0

    codon_stat_tables = {}

    def init_gene(self, gene_id, start_codons, stop_codons):
        if gene_id not in self.codon_stat_tables:
            self.codon_stat_tables[gene_id] = ContigCodonStats(start_codons, stop_codons)

    def add_contig(self, gene_id, strand, start_codons, stop_codons, barcode_count = 1):
        if len(start_codons) == 0 and len(stop_codons) > 0:
            self.cover_only_stop_codon += 1
        elif len(start_codons) > 0 and len(stop_codons) == 0:
            self.cover_only_start_codon += 1
        elif len(start_codons) == 0 and len(stop_codons) == 0:
            self.cover_zero_codons += 1
        else:
            if gene_id not in self.codon_stat_tables:
                self.codon_stat_tables[gene_id] = ContigCodonStats()

            if len(start_codons) * len(stop_codons) > 1:
                if len(start_codons) > 1 and len(stop_codons) > 1:
                    self.cover_multiple_both_codons += 1
                elif len(start_codons) == 1:
                    self.cover_multiple_stop_codons += 1
                else:
                    self.cover_multiple_start_codons += 1
                #remove this break to count ambiguous
                return
            else:
                self.cover_two_codons += 1
            if strand == "+":
                self.codon_stat_tables[gene_id].add_contig(sorted(start_codons)[0], sorted(stop_codons)[0], barcode_count)
            elif strand == "-":
                self.codon_stat_tables[gene_id].add_contig(sorted(start_codons)[-1], sorted(stop_codons)[-1], barcode_count)
            else:
                sys.stderr.write("Undefined strand")


    def add_annotated_codons(self, gene_id, annotated_codon_pairs):
        if gene_id not in self.codon_stat_tables:
            self.codon_stat_tables[gene_id] = ContigCodonStats()
        self.codon_stat_tables[gene_id].add_annotated_codons(annotated_codon_pairs)

    def count_skipped(self, genome_wide_codons):
        for c in genome_wide_codons.chr_dict.values():
            self.skipped_genes += c.skipped_genes


    def print_report(self):
        print("Skipped genes\t" + str(self.skipped_genes))
        print("Contigs covering multiple genes (ambiguous)\t" + str(self.cover_multiple_genes))
        print("Contigs covering zero genes\t" + str(self.cover_zero_genes))
        print("Contigs covering only start codons\t" + str(self.cover_only_start_codon))
        print("Contigs covering only stop codons\t" + str(self.cover_only_stop_codon))
        print("Contigs covering no codons\t" + str(self.cover_zero_codons))
        print("Contigs covering multiple start codons\t" + str(self.cover_multiple_start_codons))
        print("Contigs covering multiple stop codons\t" + str(self.cover_multiple_stop_codons))
        print("Contigs covering multiple start & stop codons\t" + str(self.cover_multiple_both_codons))
        print("Contigs covering one pair of codons\t" + str(self.cover_two_codons))

        print("Exon coordination tables with barcode count")
        for k in self.codon_stat_tables.keys():
            if dict_sum(self.codon_stat_tables[k].codon_table) > 0:
                print("=" + k)
                print(table_to_simple_str(self.codon_stat_tables[k].codon_table))



def get_barcode_map(sam_file_name):
    barcode_map = {}
    contigs_name, ext = os.path.splitext(sam_file_name)
    barcode_map_file = contigs_name + "_map.txt"
    for line in open(barcode_map_file):
        tokens = line.strip().split("_barcodeIDs_")
        if len(tokens) != 2:
            sys.stderr.write("Wrong format, _barcodeIDs_ was not found in " + line)
            continue
        barcode_map[tokens[0]] = tokens[1].replace("_", "-").split(DELIM)
    return barcode_map


def process_sam(samfile_in, codon_info, barcode_map):
    stats = Stats()
    counter = 0

    for r in samfile_in:
        counter += 1
        if counter % 10000 == 0:
           sys.stderr.write("\r   " + str(counter) + " lines processed")

        if r.reference_id == -1:
            continue
        chr_name = r.reference_name.strip()
        if chr_name not in codon_info.chr_dict:
            stats.cover_zero_genes += 1
            #sys.stderr.write("Wrong chromosome " + chr_name + "\n")
            continue

        blocks = sorted(r.get_blocks())
        contig_coord_range = [blocks[0][0], blocks[-1][1]]
        gene_index_range = codon_info.chr_dict[chr_name].find_gene_index(contig_coord_range)

        if gene_index_range == (-1, -1):
            stats.cover_zero_genes += 1
        elif gene_index_range[1] > gene_index_range[0]:
            stats.cover_multiple_genes += 1
        else:
            gene_index = gene_index_range[0]
            contig_name = r.query_name.strip()
            barcode_count = len(barcode_map[contig_name])
            gene_id = codon_info.chr_dict[chr_name].get_gene_id(gene_index)
            strand = codon_info.chr_dict[chr_name].get_strand(gene_index)
            start_codons, stop_codons =  codon_info.chr_dict[chr_name].find_covered_codons(gene_index, blocks)

            stats.init_gene(gene_id, codon_info.chr_dict[chr_name].codon_infos[gene_index].start_codons, codon_info.chr_dict[chr_name].codon_infos[gene_index].stop_codons)
            stats.add_contig(gene_id, strand, start_codons, stop_codons, barcode_count)
            stats.add_annotated_codons(gene_id, codon_info.chr_dict[chr_name].get_annotated_codons(gene_index))

    return stats


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <Codon coordinates file> <BAM/SAM file>\n")
        exit(0)

    sys.stderr.write("Reading codon structure... ")
    codon_info = GenoWideCodonStructure(sys.argv[1])
    sys.stderr.write("done.\nReading barcode file... ")
    barcode_map = get_barcode_map(sys.argv[2])
    sys.stderr.write("done.\nProcessing SAM file...\n")
    samfile_in = pysam.AlignmentFile(sys.argv[2], "rb")
    stats = process_sam(samfile_in, codon_info, barcode_map)
    sys.stderr.write("\ndone.\nPrinting report... ")
    stats.count_skipped(codon_info)
    stats.print_report()
    sys.stderr.write(" done.\n")

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()

