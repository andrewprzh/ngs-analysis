import os
import sys
from assess_contigs import *

class Feature:
    coords = (0, 0)
    feature_id = ""
    feature_type = ""
    strand = ""
    coverage = 0

    def __init__(self, feature_id, feature_type, strand, coords, cov = 0):
        self.coords = coords
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.strand = strand
        self.coverage = cov


class Gene:
    chromosome = ""
    gene_feature = None
    start_codons = set()
    stop_codons = set()
    features = {}
    contigs_info = []

    def __init__(self, feature, chromosome, start_codons, stop_codons):
        self.chromosome = chromosome
        self.gene_feature = feature
        self.start_codons = start_codons
        self.stop_codons = stop_codons
        self.features = {}
        self.contigs_info = []


    def add_contig_to_gene(self, contig_feature, covered_start_codons, covered_stop_codons, contig_alignment):
        contig_id = contig_feature.feature_id
        if contig_id not in self.features:
            self.features[contig_id] = []
        contig = self.features[contig_id]
        contig.append(contig_feature)

        strand = contig_feature.strand
        blocks = sorted(contig_alignment.get_blocks())
        for b in blocks:
            contig.append(Feature("block_" + str(b[0]), "block", strand, b))
        for c in covered_start_codons:
            contig.append(Feature("start_" + str(c[0]), "start_codon", strand, c))
        for c in covered_stop_codons:
            contig.append(Feature("stop_" + str(c[0]), "stop_codon", strand, c))

        self.contigs_info.append((covered_start_codons[0][0], covered_stop_codons[0][0], blocks[0][0], blocks[-1][1], contig_id))


    def to_string(self):
        s = "gene" + '\t' +  self.gene_feature.feature_id  + '\t' +  self.chromosome + '\t' + self.gene_feature.strand  + '\t'  + str(self.gene_feature.coords[0]) + '\t'  + str(self.gene_feature.coords[1]) + '\t' + str(len(self.features.keys())) + '\n'
        for c in self.start_codons:
            s += "start_codon"  + '\t' + str(c[0])  + '\t' +  str(c[1]) + '\n'
        for c in self.stop_codons:
            s += "stop_codon"  + '\t' + str(c[0])  + '\t' +  str(c[1]) + '\n'
        for start_c, stop_c, start_coord, stop_coord, contig_id in sorted(self.contigs_info):
            for feature in self.features[contig_id]:
                if feature.feature_type == "contig":
                    s += feature.feature_type  + '\t' + str(feature.coords[0])  + '\t' +  str(feature.coords[1]) + '\t' + str(feature.coverage) + '\n'
                else:
                    s += feature.feature_type  + '\t' + str(feature.coords[0])  + '\t' +  str(feature.coords[1]) + '\n'
        return s


class ContigFeatures:
    genes = {}

    def __init__(self):
        self.genes = {}

    def init_gene(self, gene_feature, chr_name, start_codons, stop_codons):
        if gene_feature.feature_id in self.genes:
            return
        self.genes[gene_feature.feature_id] = Gene(gene_feature, chr_name, start_codons, stop_codons)

    def add_contig(self, gene_id, contig_alignment, start_codons, stop_codons, barcode_count = 0):
        if gene_id not in self.genes:
            print("Gene " + gene_id + " is not initialized")
            return
        
        gene = self.genes[gene_id]
        contig_name = contig_alignment.query_name.strip()
        blocks = sorted(contig_alignment.get_blocks())
        coords = (blocks[0][0], blocks[-1][1])
        strand = "-" if contig_alignment.is_reverse else "+"
        gene.add_contig_to_gene(Feature(contig_name, "contig", strand, coords, barcode_count), start_codons, stop_codons, contig_alignment)

        
    def dump_to_file(self, fname):
        out_f = open(fname, 'w')
        for g in self.genes.values():
            out_f.write(g.to_string())
            
        out_f.close()


def convert_sam(samfile_in, codon_info, barcode_map):
    features = ContigFeatures()
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

        if gene_index_range != (-1, -1) and gene_index_range[1] == gene_index_range[0]:
            gene_index = gene_index_range[0]
            contig_name = r.query_name.strip()
            barcode_count = len(barcode_map[contig_name])
            gene_id = codon_info.chr_dict[chr_name].get_gene_id(gene_index)
            strand = codon_info.chr_dict[chr_name].get_strand(gene_index)
            start_codons, stop_codons =  codon_info.chr_dict[chr_name].get_all_codons(gene_index)
            if strand == "-":
                coords = (min(stop_codons)[0], max(start_codons)[1])
            else:
                coords = (min(start_codons)[0], max(stop_codons)[1])

            gene_feature = Feature(gene_id, "gene", strand, coords)
            features.init_gene(gene_feature, chr_name, start_codons, stop_codons)

            covered_start_codons, covered_stop_codons = codon_info.chr_dict[chr_name].find_covered_codons(gene_index, blocks)
            if len(covered_start_codons) > 0 and len(covered_stop_codons) > 0:
                features.add_contig(gene_id, r, covered_start_codons, covered_stop_codons, barcode_count)

    return features

    
def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: " + sys.argv[0] + " <Codon coordinates file> <BAM/SAM file> <output coords file>\n")
        exit(0)

    sys.stderr.write("Reading codon structure... ")
    codon_info = GenoWideCodonStructure(sys.argv[1])
    sys.stderr.write("done.\nReading barcode file... ")
    barcode_map = get_barcode_map(sys.argv[2])
    sys.stderr.write("done.\nProcessing SAM file...\n")
    samfile_in = pysam.AlignmentFile(sys.argv[2], "rb")
    features = convert_sam(samfile_in, codon_info, barcode_map)
    features.dump_to_file(sys.argv[3])
    sys.stderr.write(" done.\n")

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
