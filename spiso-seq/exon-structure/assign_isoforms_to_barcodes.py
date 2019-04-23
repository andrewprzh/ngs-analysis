import os
import sys
import gffutils
import pysam
from common import *


class FeatureVector:
    profile = []
    strategy = ""
    strategies = {"exact" : equal_ranges, "end" : overlaps_to_right, "start" : overlaps_to_left }
    
    def __init__(self, num, strategy = "start"):
        self.profile = [0 for i in range(0, num)]  
        self.strategy = strategy    

    def add_from_blocks(self, read_features, known_features): 
        read_pos = 0
        ref_pos = 0
        
        #if self.strategy == "start":
        #print read_features
        #print known_features

        features_present = [0 for i in range(0, len(known_features))]

        while read_pos < len(read_features) and ref_pos < len(known_features):
            while read_pos < len(read_features) and left_of(read_features[read_pos], known_features[ref_pos]):
                read_pos += 1
            if read_pos == len(read_features):
                break

            while  ref_pos < len(known_features) and left_of(known_features[ref_pos], read_features[read_pos]):
                ref_pos += 1
            if ref_pos == len(known_features):
                break

            if self.strategies[self.strategy](known_features[ref_pos], read_features[read_pos]):
                features_present[ref_pos] = 1
                ref_pos += 1
                #read_pos += 1
            elif overlaps(known_features[ref_pos], read_features[read_pos]):
                features_present[ref_pos] = -1
                ref_pos += 1
            elif known_features[ref_pos] < read_features[read_pos]:
                ref_pos += 1
            else:
                read_pos +=1
               
        #if self.strategy == "exact": 
        #    features_present = self.fill_gaps(features_present)
        for i in range(0, len(self.profile)):
            self.profile[i] += features_present[i]

    def fill_gaps(self, features_present):
        i = 0
        while i < len(features_present):
            if features_present[i] != 1:
                i += 1
                continue

            j = i + 1
            while j < len(features_present) and features_present[j] == 0:
                j += 1

            if j == len(features_present):
                break

            k = i + 1
            while k < j:
                features_present[k] = -1
                k += 1
            i = j
        return features_present

class BacrodeInfo:
    barcode = ""
    junctions_counts = None
    five_utrs_counts = None
    three_utrs_counts = None

    def __init__(self, bc, junctions_num, five_utrs_num, three_utrs_num, gene_strand):
        self.barcode = bc
        self.junctions_counts = FeatureVector(junctions_num, "exact")
        self.five_utrs_counts = FeatureVector(five_utrs_num, "end" if gene_strand == '+' else "start")
        self.three_utrs_counts = FeatureVector(three_utrs_num, "start" if gene_strand == '+' else "end")

    def add_read(self, alignment, known_junctions, known_five_utrs, known_three_utrs):
        blocks = sorted(alignment.get_blocks())
        if len(blocks) >= 2:
            read_junctions = []
            for i in range(0, len(blocks) - 1):
                read_junctions.append((blocks[i][1], blocks[i+1][0]))
            self.junctions_counts.add_from_blocks(read_junctions, known_junctions)

        self.five_utrs_counts.add_from_blocks(blocks, known_five_utrs)
        self.three_utrs_counts.add_from_blocks(blocks, known_three_utrs)


class GeneBarcodeInfo:
    gene_db = None
    db = None
    barcodes = {}
    junctions = []
    five_prime_utrs = []
    three_prime_utrs = []
    exons = []
    isoform_profiles = {}
    isoform_exon_profiles = {}

    def __init__(self, gene_db, db):
        self.db = db
        self.gene_db = gene_db
        self.barcodes = {}

        self.junctions = set()
        self.exons = set()
        i_junctions = {}
        i_exons = {}
        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            cur_exon = None
            i_junctions[t.id] = set()
            i_exons[t.id] = set()
            for e in self.db.children(t, order_by='start'):
                if e.featuretype != 'exon':
                    continue
                if cur_exon is None:
                    self.exons.add((e.start, e.end))
                    i_exons[t.id].add((e.start, e.end))
                    cur_exon = e
                    continue
                self.junctions.add((cur_exon.end, e.start))
                i_junctions[t.id].add((cur_exon.end, e.start))
                self.exons.add((e.start, e.end))
                i_exons[t.id].add((e.start, e.end))
                cur_exon = e
        self.junctions = sorted(list(self.junctions))
        self.exons = sorted(list(self.exons))
        print(self.junctions)
        print(self.exons)

        self.isoform_profiles = {}
        self.isoform_exon_profiles = {}
        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            self.isoform_profiles[t.id] = [-1 for i in range(0, len(self.junctions))]
            self.isoform_exon_profiles[t.id] = [-1 for i in range(0, len(self.exons))]
            for j in i_junctions[t.id]:
                pos = self.junctions.index(j)
                self.isoform_profiles[t.id][pos] = 1
            for e in i_exons[t.id]:
                pos = self.exons.index(e)
                self.isoform_exon_profiles[t.id][pos] = 1
            print(self.isoform_profiles[t.id])
            print(self.isoform_exon_profiles[t.id])

        self.five_prime_utrs = set()
        for u in self.db.children(self.gene_db, featuretype = 'five_prime_utr', order_by='start'):
            self.five_prime_utrs.add((u.start, u.end))
        self.five_prime_utrs = sorted(list(self.five_prime_utrs))
        print(self.five_prime_utrs)

        self.three_prime_utrs = set()
        for u in self.db.children(self.gene_db, featuretype = 'three_prime_utr', order_by='start'):
            self.three_prime_utrs.add((u.start, u.end))
        self.three_prime_utrs = sorted(list(self.three_prime_utrs))
        print(self.three_prime_utrs)


    def add_read(self, alignment, barcode_id):
        chr_name = alignment.reference_name.strip()
        if chr_name !=self.gene_db.seqid:
            return
        blocks = sorted(alignment.get_blocks())
        if len(blocks) == 0:
            return
        read_start = blocks[0][0]
        read_end = blocks[-1][1]

        if not overlaps((read_start, read_end), (self.gene_db.start, self.gene_db.end)):
            return

        if barcode_id not in self.barcodes:
            self.barcodes[barcode_id] = BacrodeInfo(barcode_id, len(self.junctions), len(self.five_prime_utrs), len(self.three_prime_utrs), self.gene_db.strand)
        self.barcodes[barcode_id].add_read(alignment, self.junctions, self.five_prime_utrs, self.three_prime_utrs)

    def assign_isoform(self, barcode_id):
        barcode_info = self.barcodes[barcode_id]
        print(barcode_info.junctions_counts.profile)
        cov = sum(map(abs, barcode_info.junctions_counts.profile))
        if cov < len(barcode_info.junctions_counts.profile):
            print("Low covered")
            return ""

        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)

        f = open("matrix.txt", "a+")
        matrix_exons = [7, 8, 10, 11, 12, 20]

        print(bacrode_jprofile)
        matched_isoforms = set()
        for t in self.isoform_profiles.keys():
            isoform_jprofile = self.isoform_profiles[t]
            if diff_only_present(isoform_jprofile, bacrode_jprofile) == 0:
                print("Matched " + t)
                print(isoform_jprofile)
                matched_isoforms.add(t)

        if len(matched_isoforms) == 0:
            print("Contraditory profile")
        elif len(matched_isoforms) > 1:
            print("Ambigous match")
        else:
            print("Unique match")
            t = list(matched_isoforms)[0]
            exon_profile = self.isoform_exon_profiles[t]
            sublist = []
            for i in matrix_exons:
                sublist.append(exon_profile[i])
            f.write(barcode_id + "\t" + "\t".join(map(str, sublist)) + "\n")

        f.close()
    

def test_BIN1(db, samfile_in):
    gene_db = db["ENSG00000136717"]
    gene_info = GeneBarcodeInfo(gene_db, db)
    print(gene_db.seqid, gene_db.start, gene_db.end, gene_db.strand)
    counter = 0
    for alignment in samfile_in:
        counter += 1
        if counter % 10000 == 0:
           sys.stderr.write("\r   " + str(counter) + " lines processed")

        if alignment.reference_id == -1:
            continue
        tokens = alignment.query_name.strip().split("___")
        if len(tokens) != 2:
            continue

        barcode_id = tokens[1]
        gene_info.add_read(alignment, barcode_id)

    f = open("matrix.txt", "w")
    f.close()
    counter = 0
    for b in gene_info.barcodes.keys():
        print(b)
        gene_info.assign_isoform(b)

#        print(gene_info.barcodes[b].junctions_counts.profile)
#        print(gene_info.barcodes[b].five_utrs_counts.profile)
#        print(gene_info.barcodes[b].three_utrs_counts.profile)

    
def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: " + sys.argv[0] + " <gene DB> <BAM/SAM file> <output.tsv>\n")
        exit(0)
    db = gffutils.FeatureDB(sys.argv[1], keep_order=True)
    samfile_in = pysam.AlignmentFile(sys.argv[2], "rb")
    test_BIN1(db, samfile_in)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()

