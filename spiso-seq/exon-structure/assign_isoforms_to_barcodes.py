import os
import sys
import gffutils
import pysam
from common import *

DEDUCE_CODONS_FROM_CDS = True
KEEP_ISOFORMS_WITHOUT_CODONS = False
READS_CUTOFF = 10
MIN_CODON_COUNT = 2
ASSIGN_CODONS_WHEN_AMBIGUOUS = False

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
               
        self.fill_gaps(features_present)

        self.reads += 1
        for i in range(0, len(self.profile)):
            self.profile[i] += features_present[i]

    def fill_gaps(self, features_present):
        start = 0
        while start < len(features_present) and features_present[start] == 0:
            start += 1

        end = len(features_present) - 1
        while end > 0 and features_present[end] == 0:
            end -= 1

        for i in range(start, end):
            if features_present[i] == 0:
                features_present[i] = -1


class BacrodeInfo:
    barcode = ""
    total_reads = 0
    junctions_counts = None
    #five_utrs_counts = None
    #three_utrs_counts = None

    def __init__(self, bc, junctions_num, gene_strand):
        self.barcode = bc
        self.total_reads = 0
        self.junctions_counts = FeatureVector(junctions_num, "exact")
        #self.five_utrs_counts = FeatureVector(five_utrs_num, "end" if gene_strand == '+' else "start")
        #self.three_utrs_counts = FeatureVector(three_utrs_num, "start" if gene_strand == '+' else "end")

    def add_read(self, alignment, known_junctions):
        self.total_reads += 1
        blocks = sorted(alignment.get_blocks())
        if len(blocks) >= 2:
            read_junctions = []
            for i in range(0, len(blocks) - 1):
                read_junctions.append((blocks[i][1], blocks[i+1][0]))
            self.junctions_counts.add_from_blocks(read_junctions, known_junctions)

        #self.five_utrs_counts.add_from_blocks(blocks, known_five_utrs)
        #self.three_utrs_counts.add_from_blocks(blocks, known_three_utrs)


class BarcodeAssignmentStats:
    low_covered = 0
    uniquely_assigned = 0
    assigned_to_ncrna = 0
    contradictory = 0
    empty = 0
    ambiguous = 0
    ambiguous_codon_assigned = 0
    ambiguous_subisoform_assigned = 0

    correctly_assigned = 0
    unassigned = 0
    mismapped = 0
    incorrectly_assigned_same_gene = 0
    incorrectly_assigned_other_gene = 0

    def __init__(self):
        self.low_covered = 0
        self.uniquely_assigned = 0
        self.assigned_to_ncrna = 0
        self.contradictory = 0
        self.empty = 0
        self.ambiguous = 0
        self.ambiguous_codon_assigned = 0
        self.ambiguous_subisoform_assigned = 0

        self.correctly_assigned = 0
        self.unassigned = 0
        self.mismapped = 0
        self.incorrectly_assigned_same_gene = 0
        self.incorrectly_assigned_other_gene = 0

    def merge(self, stat):
        self.low_covered += stat.low_covered
        self.uniquely_assigned += stat.uniquely_assigned
        self.assigned_to_ncrna += stat.assigned_to_ncrna
        self.contradictory += stat.contradictory
        self.empty += stat.empty
        self.ambiguous += stat.ambiguous
        self.ambiguous_codon_assigned += stat.ambiguous_codon_assigned
        self.ambiguous_subisoform_assigned += stat.ambiguous_subisoform_assigned

        self.correctly_assigned += stat.correctly_assigned
        self.unassigned += stat.unassigned
        self.mismapped += stat.mismapped
        self.incorrectly_assigned_same_gene += stat.incorrectly_assigned_same_gene
        self.incorrectly_assigned_other_gene += stat.incorrectly_assigned_other_gene

    def isoform_stats(self):
        total = self.correctly_assigned+self.unassigned+self.mismapped+self.incorrectly_assigned_same_gene+self.incorrectly_assigned_other_gene
        s = "\nTotal  correct  wrong_same  wrong_other  unassigned  mismapped\n"
        return s + "%d\t%d\t%d\t%d\t%d\t%d" % \
            (total, self.correctly_assigned, self.incorrectly_assigned_same_gene, self.incorrectly_assigned_other_gene, self.unassigned, self.mismapped)

    def to_str(self):
        total_bc = self.low_covered + self.uniquely_assigned + self.assigned_to_ncrna + self.contradictory + self.empty + self.ambiguous + self.ambiguous_codon_assigned + self.ambiguous_subisoform_assigned
        s = "\nTotal  low_covered  unique  ncrna  contradictory  empty  ambiguous  ambiguous_codon  ambiguous_assigned:\n"
        return s + "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t" % \
            (total_bc, self.low_covered, self.uniquely_assigned, self.assigned_to_ncrna, self.contradictory, self.empty, self.ambiguous, self.ambiguous_codon_assigned, self.ambiguous_subisoform_assigned)

class ProfileStorage:
    isoform_profiles = {}
    isoform_exon_profiles = {}

    def __init__(self):
        self.isoform_profiles = {}
        self.isoform_exon_profiles = {}
         

class GeneBarcodeInfo:
    gene_db = None
    db = None
    coding_rna_profiles = ProfileStorage()
    all_rna_profiles = ProfileStorage()
    junctions = []
    exons = []
    codon_pairs = {}
    barcodes = {}

    def get_codon_pair(self, transcript):
        start_codon = None
        stop_codon = None
        for s in self.db.children(transcript, featuretype='start_codon', order_by='start'):
            start_codon = s.start
        for s in self.db.children(transcript, featuretype='stop_codon', order_by='start'):
            stop_codon = s.start

        if not DEDUCE_CODONS_FROM_CDS:
            return start_codon, stop_codon

        if start_codon is None:
            for s in self.db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    start_codon = s.start
                else:
                    start_codon = s.end
                break
        if stop_codon is None:
            for s in self.db.children(transcript, featuretype='CDS', order_by='start'):
                if s.strand == "+":
                    stop_codon = s.end + 1
                else:
                    stop_codon = s.start - 2
        return start_codon, stop_codon


    def isoform_is_coding(self, t):
        start_codon, stop_codon = self.get_codon_pair(t)
        return stop_codon is not None and start_codon is not None


    def get_junctions_and_exons(self, keep_isoforms_without_codons = False):
        i_junctions = {}
        i_exons = {}
        codon_pairs_map = {}

        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                continue

            cur_exon = None
            i_junctions[t.id] = set()
            i_exons[t.id] = set()
            for e in self.db.children(t, order_by='start'):
                if e.featuretype == 'exon':
                    if cur_exon is None:
                        i_exons[t.id].add((e.start, e.end))
                        cur_exon = e
                        continue
                    i_junctions[t.id].add((cur_exon.end, e.start))
                    i_exons[t.id].add((e.start, e.end))
                    cur_exon = e

        self.junctions = set()
        self.exons = set()
        for i in i_junctions.keys():
            for j in i_junctions[i]:
                self.junctions.add(j)
        for i in i_exons.keys():
            for e in i_exons[i]:
                self.exons.add(e)
            
        self.junctions = sorted(list(self.junctions))
        self.exons = sorted(list(self.exons))
        print(self.junctions)
        print(self.exons)

        return i_junctions, i_exons


    def set_junction_profiles(self, profile_storage, i_junctions, i_exons, keep_isoforms_without_codons = False):
        profile_storage.isoform_profiles = {}
        profile_storage.isoform_exon_profiles = {}

        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            if not keep_isoforms_without_codons and not self.isoform_is_coding(t):
                continue

            profile_storage.isoform_profiles[t.id] = [-1 for i in range(0, len(self.junctions))]
            profile_storage.isoform_exon_profiles[t.id] = [-1 for i in range(0, len(self.exons))]
            for j in i_junctions[t.id]:
                pos = self.junctions.index(j)
                profile_storage.isoform_profiles[t.id][pos] = 1
            for e in i_exons[t.id]:
                pos = self.exons.index(e)
                profile_storage.isoform_exon_profiles[t.id][pos] = 1

            print("Isoform " + t.id)
            print(profile_storage.isoform_profiles[t.id])
            print(profile_storage.isoform_exon_profiles[t.id])

            if all(x == -1 for x in profile_storage.isoform_profiles):
                del profile_storage.isoform_profiles[t.id]
                del profile_storage.isoform_exon_profiles[t.id]


    def set_codon_pairs(self):
        for t in self.db.children(self.gene_db, featuretype = 'transcript', order_by='start'):
            start_codon, stop_codon = self.get_codon_pair(t)
            if stop_codon is not None and start_codon is not None:
                self.codon_pairs[t.id] = (start_codon, stop_codon)            


    def __init__(self, gene_db, db):
        self.codon_pairs = {}
        self.db = db
        self.gene_db = gene_db
        self.barcodes = {}
        self.junctions = []
        self.exons = []
        self.coding_rna_profiles = ProfileStorage()
        self.all_rna_profiles = ProfileStorage()

        self.set_codon_pairs()
        i_junctions, i_exons = self.get_junctions_and_exons(True)
        print("Coding isoforms")
        self.set_junction_profiles(self.coding_rna_profiles, i_junctions, i_exons)
        print("NC isoforms")
        self.set_junction_profiles(self.all_rna_profiles, i_junctions, i_exons, True)

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
            self.barcodes[barcode_id] = BacrodeInfo(barcode_id, len(self.junctions), self.gene_db.strand)
        self.barcodes[barcode_id].add_read(alignment, self.junctions)


    def find_matches(self, barcode_info, profile_storage):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
        matched_isoforms = set()
        for t in profile_storage.isoform_profiles.keys():
            isoform_jprofile = profile_storage.isoform_profiles[t]
            if diff_only_present(isoform_jprofile, bacrode_jprofile) == 0:
                print("Matched " + t)
                matched_isoforms.add(t)
        return matched_isoforms


    def is_empty_alignment(self, barcode_info):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
        print(barcode_info.barcode)
        print(bacrode_jprofile)

        if all(el != 1 for el in bacrode_jprofile):
            return True
        return False


    def assign_isoform(self, barcode_id, stat, coverage_cutoff):
        barcode_info = self.barcodes[barcode_id]
        if barcode_info.total_reads < coverage_cutoff:
            stat.low_covered += 1
            return None, None          

        if self.is_empty_alignment(barcode_info):
            stat.empty += 1
            return None, None          

        matched_isoforms = self.find_matches(barcode_info, self.coding_rna_profiles)

        transcript_id = None
        codon_pair = None
        if len(matched_isoforms) == 0:
            transcript_id = self.resolve_contradictory(barcode_info, stat)

        elif len(matched_isoforms) > 1:
            matched_isoforms = self.resolve_ambiguous(barcode_info, matched_isoforms, self.coding_rna_profiles)
            
            if len(matched_isoforms) == 1:
                stat.ambiguous_subisoform_assigned += 1
                print("Unique match after resolution")
                transcript_id = list(matched_isoforms)[0]
                codon_pair = self.codon_pairs[transcript_id]
            else:
                codons = set()
                for t in matched_isoforms:
                    codons.add(self.codon_pairs[t])
                if ASSIGN_CODONS_WHEN_AMBIGUOUS and len(codons) == 1:
                    codon_pair = list(codons)[0]
                    stat.ambiguous_codon_assigned += 1   
                    transcript_id = list(matched_isoforms)[0]
                else:                
                    stat.ambiguous += 1
                print("Ambigous match")
        else:
            stat.uniquely_assigned += 1
            print("Unique match")
            transcript_id = list(matched_isoforms)[0]
            codon_pair = self.codon_pairs[transcript_id]
        return transcript_id, codon_pair


    def resolve_contradictory(self, barcode_info, stat):
        matched_isoforms = self.find_matches(barcode_info, self.all_rna_profiles)
        transcript_id = None
        if len(matched_isoforms) == 0:
            stat.contradictory += 1
            print("Contraditory profile")
        elif len(matched_isoforms) == 1:
            stat.assigned_to_ncrna += 1
            print("Non-coding assigned")
            transcript_id = list(matched_isoforms)[0]
        else:
            matched_isoforms = self.resolve_ambiguous(barcode_info, matched_isoforms, self.all_rna_profile)
            if len(matched_isoforms) == 1:
                stat.assigned_to_ncrna += 1
                transcript_id = list(matched_isoforms)[0]
            else:
                stat.ambiguous += 1

        return transcript_id


    def resolve_ambiguous(self, barcode_info, matched_isoforms, profile_storage):
        bacrode_jprofile = map(sign, barcode_info.junctions_counts.profile)
        for t in matched_isoforms:
            matched_positions = find_matching_positions(profile_storage.isoform_profiles[t], bacrode_jprofile)
            
            all_junctions_detected = True
            for i in range(len(matched_positions)):
                if matched_positions == 0 and profile_storage.isoform_profiles[t] != -1:
                    all_junctions_detected = False
                    break

            if all_junctions_detected:
                print("Ambiguity resolved")
                print(profile_storage.isoform_profiles[t])
                print(matched_positions)
                print(bacrode_jprofile)
                return set([t])

        return matched_isoforms


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


def get_gene_barcodes(db, gene_info, samfile_name, total_stats, is_reads_sam, bc_map):
    samfile_in = pysam.AlignmentFile(samfile_name, "rb")
    gene_chr, gene_start, gene_end = gene_info.get_gene_region()

    counter = 0
    for alignment in samfile_in.fetch(gene_chr, gene_start, gene_end):
        counter += 1
        if counter % 100 == 0:
           sys.stderr.write("\r   " + str(counter) + " lines processed")

        if alignment.reference_id == -1:
            continue
        seq_id = get_id(alignment.query_name, is_reads_sam)
        gene_info.add_read(alignment, seq_id)

    #bc_map = get_barcode_map(samfile_name) if not is_reads_sam else None
    gene_isoform_ids = set(gene_info.all_rna_profiles.isoform_profiles.keys())
    barcodes = {}
    stats = BarcodeAssignmentStats()
    for b in gene_info.barcodes.keys():
        print(" ===== ")
        cutoff = READS_CUTOFF if is_reads_sam else 0
        isoform, codons = gene_info.assign_isoform(b, stats, cutoff)
        if isoform is not None:
            if b == isoform:
                stats.correctly_assigned += 1
            elif b in gene_isoform_ids:
                stats.incorrectly_assigned_same_gene += 1
            else:
                stats.incorrectly_assigned_other_gene += 1
                
            if bc_map is None:
                barcodes[b] = (isoform, codons)
            else:
                for bc in bc_map[b]:
                    barcodes[bc] = (isoform, codons)
        else:
            if b in gene_isoform_ids:
                stats.unassigned += 1
            else:
                stats.mismapped += 1
            

    sys.stderr.write("\nDone. Barcodes stats " + stats.to_str() + "\n")
    sys.stderr.write("\nDone. Isoform stats " + stats.isoform_stats() + "\n")
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
    start_codons = set()
    stop_codons = set()
    for b in barcodes.keys():
        codon_pair = barcodes[b][1]
        if codon_pair is None or codon_pair[0] is None or codon_pair[1] is None:
            continue
        start_codons.add(codon_pair[0])
        stop_codons.add(codon_pair[1])
        if codon_pair not in codon_count_table:
            codon_count_table[codon_pair] = 0
        codon_count_table[codon_pair] += 1

    if len(start_codons) >= MIN_CODON_COUNT and len(stop_codons) >= MIN_CODON_COUNT:
        outf = open(out_codon_stats, "a+")
        outf.write("====" + gene_name + "\n")
        outf.write(table_to_str(codon_count_table))
        outf.close()


def need_to_process(gene_db, db):
    #checking codons
    start_codons = set()
    stop_codons = set()
    for t in db.children(gene_db, featuretype='transcript', order_by='start'):
        start_codon = None
        stop_codon = None
        for s in db.children(t, featuretype='start_codon', order_by='start'):
            start_codon = s.start
        for s in db.children(t, featuretype='stop_codon', order_by='start'):
            stop_codon = s.start

        if start_codon is not None and stop_codon is not None:
            start_codons.add(start_codon)
            stop_codons.add(stop_codon)

    return len(start_codons) >= MIN_CODON_COUNT and len(stop_codons) >= MIN_CODON_COUNT


def process_all_genes(db, samfile_name, outf_prefix, is_reads_sam = True):
    out_tsv = outf_prefix + ".barcodes.tsv"
    outf = open(out_tsv, "w")
    outf.close()
    out_codon_stats = outf_prefix + ".codon_stats.tsv"
    outf = open(out_codon_stats, "w")
    outf.close()
    stats = BarcodeAssignmentStats()
    bc_map = get_barcode_map(samfile_name) if not is_reads_sam else None

    for g in db.features_of_type('gene', order_by='start'):
        gene_name = g.id

        gene_db = db[gene_name]
        if not need_to_process(gene_db, db):
            continue

        gene_info = GeneBarcodeInfo(gene_db, db)
        sys.stderr.write("Processing gene " + gene_name + "\n")

        barcodes = get_gene_barcodes(db, gene_info, samfile_name, stats, is_reads_sam, bc_map)
        write_gene_stats(db, gene_name, barcodes, out_tsv, out_codon_stats)

    sys.stderr.write("\nFinished. Total stats " + stats.to_str() + "\n")
    sys.stderr.write("\nFinished. Isoform stats " + stats.isoform_stats() + "\n")

    
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

