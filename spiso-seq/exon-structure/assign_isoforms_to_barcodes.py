import os
import sys
import gffutils
import pysam
from common import *

DEDUCE_CODONS_FROM_CDS = True
KEEP_ISOFORMS_WITHOUT_CODONS = False
READS_CUTOFF = 10

def print_ints(l):
    print("\t".join(map(str,l)))


def list_by_indices(indices, l):
    sublist = []
    for i in indices:
        sublist.append(l[i])
    return sublist

class FeatureVector:
    profile = []
    reads = 0
    strategy = ""
    strategies = {"exact" : equal_ranges, "end" : overlaps_to_right, "start" : overlaps_to_left }
    
    def __init__(self, num, strategy = "start"):
        self.profile = [0 for i in range(0, num)]  
        self.reads = 0
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
        self.reads += 1
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
    total_reads = 0
    junctions_counts = None
    five_utrs_counts = None
    three_utrs_counts = None

    def __init__(self, bc, junctions_num, five_utrs_num, three_utrs_num, gene_strand):
        self.barcode = bc
        self.total_reads = 0
        self.junctions_counts = FeatureVector(junctions_num, "exact")
        self.five_utrs_counts = FeatureVector(five_utrs_num, "end" if gene_strand == '+' else "start")
        self.three_utrs_counts = FeatureVector(three_utrs_num, "start" if gene_strand == '+' else "end")

    def add_read(self, alignment, known_junctions, known_five_utrs, known_three_utrs):
        self.total_reads += 1
        blocks = sorted(alignment.get_blocks())
        if len(blocks) >= 2:
            read_junctions = []
            for i in range(0, len(blocks) - 1):
                read_junctions.append((blocks[i][1], blocks[i+1][0]))
            self.junctions_counts.add_from_blocks(read_junctions, known_junctions)

        self.five_utrs_counts.add_from_blocks(blocks, known_five_utrs)
        self.three_utrs_counts.add_from_blocks(blocks, known_three_utrs)


class BarcodeAssignmentStats:
    low_covered = 0
    contradictory = 0
    ambiguous = 0
    ambiguous_assigned = 0
    unique = 0

    def __init__(self):
        self.low_covered = 0
        self.contradictory = 0
        self.ambiguous = 0
        self.ambiguous_assigned = 0
        self.unique = 0

    def merge(self, stat):
        self.low_covered += stat.low_covered
        self.contradictory += stat.contradictory
        self.ambiguous += stat.ambiguous
        self.ambiguous_assigned += stat.ambiguous_assigned 
        self.unique += stat.unique

    def to_str(self):
        total_bc = self.low_covered + self.contradictory + self.ambiguous + self.unique
        return "low covered/contradictory/ambiguous/ambiguous_assigned/unique/total: %d / % d / % d / % d / % d / % d" % \
            (self.low_covered, self.contradictory, self.ambiguous, self.ambiguous_assigned, self.unique, total_bc)
         

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
    codon_pairs = {}

    def get_codon_pair(self, transcript, db):
        start_codon = None
        stop_codon = None
        for s in db.children(transcript, featuretype='start_codon', order_by='start'):
            start_codon = s.start
        for s in db.children(transcript, featuretype='stop_codon', order_by='start'):
            stop_codon = s.start

        if not DEDUCE_CODONS_FROM_CDS:
            return start_codon, stop_codon

        if start_codon is None:
            for s in db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    start_codon = s.start
                else:
                    start_codon = s.end
                break
        if stop_codon is None:
            for s in db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    stop_codon = s.end + 1
                else:
                    stop_codon = s.start - 2

        return start_codon, stop_codon
        

    def __init__(self, gene_db, db):
        self.codon_pairs = {}
        self.db = db
        self.gene_db = gene_db
        self.barcodes = {}
        self.junctions = set()
        self.exons = set()
        i_junctions = {}
        i_exons = {}
        i_5utr = {}
        i_3utr = {}

        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            start_codon, stop_codon = self.get_codon_pair(t, db)

            if not KEEP_ISOFORMS_WITHOUT_CODONS and (stop_codon is None or start_codon is None):
                continue

            self.codon_pairs[t.id] = (start_codon, stop_codon)

            cur_exon = None
            i_junctions[t.id] = set()
            i_exons[t.id] = set()
            for e in self.db.children(t, order_by='start'):
                if e.featuretype == 'five_prime_utr':
                    i_5utr[t.id] = (e.start, e.end)
                elif e.featuretype == 'three_prime_utr':
                    i_3utr[t.id] = (e.start, e.end)
                elif e.featuretype == 'exon':
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
        #print(self.junctions)
        #print(self.exons)

        self.five_prime_utrs = set()
        for u in self.db.children(self.gene_db, featuretype = 'five_prime_utr', order_by='start'):
            self.five_prime_utrs.add((u.start, u.end))
        self.five_prime_utrs = sorted(list(self.five_prime_utrs))
        #print(self.five_prime_utrs)

        self.three_prime_utrs = set()
        for u in self.db.children(self.gene_db, featuretype = 'three_prime_utr', order_by='start'):
            self.three_prime_utrs.add((u.start, u.end))
        self.three_prime_utrs = sorted(list(self.three_prime_utrs))
        #print(self.three_prime_utrs)

        self.isoform_profiles = {}
        self.isoform_exon_profiles = {}
        for t_id in self.codon_pairs.keys():
            self.isoform_profiles[t_id] = [-1 for i in range(0, len(self.junctions))]
            self.isoform_exon_profiles[t_id] = [-1 for i in range(0, len(self.exons))]
            for j in i_junctions[t_id]:
                pos = self.junctions.index(j)
                self.isoform_profiles[t_id][pos] = 1
            for e in i_exons[t_id]:
                pos = self.exons.index(e)
                self.isoform_exon_profiles[t_id][pos] = 1

#            index_5utr = self.three_prime_utrs.index(i_5utr[t_id])
#            self.utr_indices(
#            self.isoform_exon_profiles[t_id][pos] = 1

            #print(self.isoform_profiles[t_id])
            #print(self.isoform_exon_profiles[t_id])

        #check for identical profiles
        #isoform_profiles_tuples = map(tuple, self.isoform_profiles.values())
        #if len(isoform_profiles_tuples) > len(set(isoform_profiles_tuples)):
        #    print("Identical junction profiles for different isoforms")




    def get_gene_region(self):
        return self.gene_db.seqid, self.gene_db.start, self.gene_db.end       


    def add_read(self, alignment, barcode_id):
        chr_name = alignment.reference_name.strip()
        if chr_name != self.gene_db.seqid:
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


    def assign_isoform(self, barcode_id, stat, coverage_cutoff):
        barcode_info = self.barcodes[barcode_id]
        if barcode_info.total_reads < coverage_cutoff:
            stat.low_covered += 1
            return None, None          

        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)

        #print(bacrode_jprofile)
        matched_isoforms = set()
        for t in self.isoform_profiles.keys():
            isoform_jprofile = self.isoform_profiles[t]
            if diff_only_present(isoform_jprofile, bacrode_jprofile) == 0:
                #print("Matched " + t)
                #print(isoform_jprofile)
                matched_isoforms.add(t)

        transcript_id = None
        codon_pair = None
        if len(matched_isoforms) == 0:
            stat.contradictory += 1
            #print("Contraditory profile")
        elif len(matched_isoforms) > 1:
            codons = set()
            for t in matched_isoforms:
                codons.add(self.codon_pairs[t])
            if len(codons) == 1:
                codon_pair = list(codons)[0]
                stat.ambiguous_assigned += 1   
                transcript_id = list(matched_isoforms)[0]
            else:                
                stat.ambiguous += 1
            #print("Ambigous match")
        else:
            stat.unique += 1
            #print("Unique match")
            transcript_id = list(matched_isoforms)[0]
            codon_pair = self.codon_pairs[transcript_id]
        return transcript_id, codon_pair


def get_barcode_map(sam_file_name):
    barcode_map = {}
    contigs_name, ext = os.path.splitext(sam_file_name)
    barcode_map_file = contigs_name + "_map.txt"
    for line in open(barcode_map_file):
        tokens = line.strip().split("_barcodeIDs_")
        if len(tokens) != 2:    
            print("Wrong format, _barcodeIDs_ was not found in " + line)
            continue
        barcode_map[tokens[0]] = tokens[1].replace("_", "-").split(',')
    return barcode_map

#needs changing
def get_id(query_name, is_reads_sam = True):
    if is_reads_sam:
        tokens = query_name.strip().split("___")
        if len(tokens) != 2:
            return ""

        return tokens[1]
    else:
        return query_name.strip()


def get_gene_barcodes(db, gene_info, samfile_name, total_stats, is_reads_sam):
    samfile_in = pysam.AlignmentFile(samfile_name, "rb")
    gene_chr, gene_start, gene_end = gene_info.get_gene_region()

    counter = 0
    for alignment in samfile_in.fetch(gene_chr, gene_start, gene_end):
        counter += 1
        if counter % 10000 == 0:
           sys.stderr.write("\r   " + str(counter) + " lines processed")

        if alignment.reference_id == -1:
            continue
        seq_id = get_id(alignment.query_name, is_reads_sam)
        gene_info.add_read(alignment, seq_id)

    bc_map = get_barcode_map(samfile_name) if not is_reads_sam else None
    barcodes = {}
    stats = BarcodeAssignmentStats()
    for b in gene_info.barcodes.keys():
        cutoff = READS_CUTOFF if is_reads_sam else 0
        isoform, codons = gene_info.assign_isoform(b, stats, cutoff)
        if isoform is not None:
            if bc_map is None:
                barcodes[b] = (isoform, codons)
            else:
                for bc in bc_map[b]:
                    barcodes[bc] = (isoform, codons)

    sys.stderr.write("\nDone. Barcodes stats " + stats.to_str() + "\n")
    total_stats.merge(stats)
    return barcodes

def table_to_str(d):
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
        res += s[:-1] + '\n'
    return res


def write_gene_stats(db, gene_name, barcodes, out_tsv, out_codon_stats):
    #writing TSV with barcode -> isoform id
    outf = open(out_tsv, "a+")
    outf.write(gene_name + "\t" + str(len(barcodes)) + "\n")
    for b in barcodes.keys():
        if b != "":
            outf.write(b + "\t" + barcodes[b][0] + "\n")
    outf.close()

    #writing codon stats
    codon_count_table = {}
    for b in barcodes.keys():
        codon_pair = barcodes[b][1]
        if codon_pair[0] is None or codon_pair[1] is None:
            continue
        if codon_pair not in codon_count_table:
            codon_count_table[codon_pair] = 0
        codon_count_table[codon_pair] += 1

    outf = open(out_codon_stats, "a+")
    outf.write("====" + gene_name + "\n")
    outf.write(table_to_str(codon_count_table))
    outf.close()


def process_all_genes(db, samfile_name, outf_prefix, is_reads_sam = True):
    out_tsv = outf_prefix + ".barcodes.tsv"
    outf = open(out_tsv, "w")
    outf.close()
    out_codon_stats = outf_prefix + ".codon_stats.tsv"
    outf = open(out_codon_stats, "w")
    outf.close()
    stats = BarcodeAssignmentStats()
    for g in db.features_of_type('gene', order_by='start'):
        gene_name = g.id

        #checking codons
        gene_db = db[gene_name]
        start_codons = set()
        stop_codons = set()
        for t in db.children(gene_db, featuretype='transcript', order_by='start'):
            start_codon = -1
            stop_codon = -1
            for s in db.children(t, featuretype='start_codon', order_by='start'):
                start_codon = s.start
            for s in db.children(t, featuretype='stop_codon', order_by='start'):
                stop_codon = s.start

            if start_codon != -1 and stop_codon != -1:
                start_codons.add(start_codon)
                stop_codons.add(stop_codon)

        if len(start_codons) <= 1 or len(stop_codons) <= 1:
            continue

        gene_db = db[gene_name] 
        gene_info = GeneBarcodeInfo(gene_db, db)
        sys.stderr.write("Processing gene " + gene_name + "\n")

        barcodes = get_gene_barcodes(db, gene_info, samfile_name, stats, is_reads_sam)
        write_gene_stats(db, gene_name, barcodes, out_tsv, out_codon_stats)

    sys.stderr.write("\nFinished. Total stats " + stats.to_str() + "\n")

    
def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: " + sys.argv[0] + " <gene DB> <BAM/SAM file> <output.tsv> [ READS / contgs ]\n")
        exit(0)

    is_read_samfile = False if len(sys.argv) == 5 and sys.argv[4].upper().startswith("C") else True
    db = gffutils.FeatureDB(sys.argv[1], keep_order = True)
    process_all_genes(db, sys.argv[2], sys.argv[3], is_read_samfile)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()

