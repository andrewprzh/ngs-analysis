import os
import sys
import pysam
import common
from common import *

DELIM = ','

dbg_barcodes = set() #['BX:Z:CTAGGTGTATAACG-1', 'BX:Z:AACGCTTCTTGAGC-1', 'BX:Z:TGCCAAGTTTACAT-1', 'BX:Z:CTAGTGAGATTCTT-1', 'BX:Z:GAGCGAGTCTAAAG-1', 'BX:Z:AGTACCAGCTGGAC-1', 'BX:Z:AGTACCAGCTGGAC-1', 'BX:Z:CTTAGGTCGTGGCC-1', 'BX:Z:ATTTCCCATCGCAT-1', 'BX:Z:AGCTGTAGGTTGAT-1', 'BX:Z:CATCCCAGTGTCCT-1'])
#['BX:Z:ACGGGAAGGCCATG-1'])


def debug(bc, s):
    if bc in dbg_barcodes:
        print(s)

def process_read(alignment, exons, barcode):
    exons_present = [0 for i in range(2 + len(exons))]

    blocks = sorted(alignment.get_blocks())
    block_pos = 0
    exon_pos = 0

    if barcode in dbg_barcodes:
        print blocks
        print exons

    
    while block_pos < len(blocks) and exon_pos < len(exons):
        while block_pos < len(blocks) and left_of(blocks[block_pos], exons[exon_pos]):
            block_pos += 1
        if block_pos == len(blocks):
            break

        while exon_pos < len(exons) and left_of(exons[exon_pos], blocks[block_pos]):
            if block_pos == 0:
                exons_present[exon_pos + 1] = 0
            else:
                exons_present[exon_pos + 1] = -1
            exon_pos += 1
        if exon_pos == len(exons):
            break

        if overlaps(exons[exon_pos], blocks[block_pos]):
            exons_present[exon_pos + 1] = 1
            exon_pos += 1
            block_pos += 1

    rr = ReadRecord(blocks, exons_present)
#    rr.FillGaps()

    debug(barcode, str(barcode) + ": " + str(exons_present))
    return barcode, rr


def merge_exon_info(bc_record):
    if len(exon_info) != len(new_info):
        print("Cannot merge exon information")
        return exon_info
    result = [0 for i in range(len(exon_info))]

    for i in range(len(exon_info)):
        if exon_info[i] == 0 and new_info[i] != 0:
            result[i] = new_info[i]
        elif new_info[i] == 0:
            result[i] = exon_info[i]
        elif exon_info[i] != 0 and exon_info[i] != new_info[i]:
            print("Information contradicts: " + str(exon_info) + ", " + str(new_info))
            return None

    return result


def read_matrix(matrix_file):
    hagen_dict = {}
    exons_number = len(initial_coords)
    with open(maxtrix_file, "r") as hagen_file:
        for line in hagen_file.readlines():
            tokens = line.split("\t")
            barcode = tokens[0]
            hagen_dict[barcode] = [int(tokens[i]) for i in range(1,exons_number + 1)]
    return hagen_dict

def process_sam(samfile_in, exon_coords, barcode_map, matrix_file):
    hagen_dict = read_matrix(matrix_file)
    barcode_info = {}
    alignments = []
    contigs = set()
    used_bc = {}
    for r in samfile_in:
        contig_name = r.query_name.strip()
        if contig_name in contigs:
            pass
            #print contig_name
        else:
            contigs.add(contig_name)

        alignments.append(contig_name)
        if contig_name not in barcode_map:
            print("No info on " + contig_name)
            continue

        barcodes = barcode_map[contig_name]
        
        for barcode in barcodes:
            barcode, read_record = process_read(r, exon_coords, barcode)
            if barcode not in barcode_info:
                barcode_info[barcode] = BarcodeRecord(barcode, [read_record])
                used_bc[barcode] = [contig_name]
            else:
                barcode_info[barcode].AddAlignment(read_record)
                used_bc[barcode].append(contig_name)

    for bc, ctgs in used_bc.iteritems():
        if len(ctgs) > 1 and bc in hagen_dict:
            pass#print bc, " ".join(ctgs)

    return barcode_info, alignments



def get_barcode_map(sam_file_name):
    barcode_map = {}
    contigs_name, ext = os.path.splitext(sam_file_name)
    barcode_map_file = contigs_name + "_map.txt"
    for line in open(barcode_map_file):
        tokens = line.strip().split("_barcodeIDs_")
        if len(tokens) != 2:
            print("Wrong format, _barcodeIDs_ was not found in " + line)
            continue
        barcode_map[tokens[0]] = tokens[1].replace("_", "-").split(DELIM)
    return barcode_map


def assess_alignment(barcode, alignment, hagen_dict, stats, unique_isoforms):
    if not alignment.Complete():
        stats.incomplete.append(barcode)
        if alignment.exon_present[1] == 0 and alignment.exon_present[-2] != 0:
            stats.incomplete5.append(barcode)
        elif alignment.exon_present[1] != 0 and alignment.exon_present[-2] == 0:
            stats.incomplete3.append(barcode)

        if barcode not in hagen_dict.keys():
            stats.incomplete_missing.append(barcode)
            return
        
        stats.incomplete_in_dict.append(barcode)
        d = diff_only_present(alignment.exon_present[1:-1], hagen_dict[barcode])
        if d == 0:
            stats.incomplete_correct.append(barcode)
        else:
            stats.incomplete_incorrect.append(barcode)
        return

    unique_isoforms.add("".join(str(x) for x in alignment.exon_present[1:-1]))

    if barcode not in hagen_dict.keys():
        stats.contigs_only.append(barcode)
        return
    
    d = diff_only_present(alignment.exon_present[1:-1], hagen_dict[barcode])
    if d == 0:
        stats.correct.append(barcode)
    else:
        if d == 1:
            stats.incorrect1.append(barcode)
            #print(str(barcode) + " with D = " + str(d) + ":\n " + str(alignment.exon_present[1:-1])+ "\n " + str(hagen_dict[barcode]))
            #print(str(initial_coords))
            #print(str(alignment.block_coordinates))
        stats.incorrect.append(barcode)


def get_stats(maxtrix_file, initial_coords, barcode_info, barcode_map):
    stats = BarcodeStats()
    secondary_stats = BarcodeStats()
    hagen_dict = read_matrix(maxtrix_file)

    print("Total barcodes in maxtrix " + str(len(hagen_dict)))
    stats.total_matrix_barcodes = len(hagen_dict)
    stats.contig_alignment_count = len(barcode_map)
    stats.distinct_barcodes_in_contigs = len(barcode_info.items())

    unique_assembled_isoforms = set()

    for barcode in barcode_info.keys():
        bc_info = barcode_info[barcode]
        stats.used.add(barcode)
        if len(bc_info.read_infos) > 1:
            stats.with_secondary.append(barcode)

        contig_alignments = sorted(bc_info.read_infos, reverse = True, key = lambda x: x.InformativeExons())

        complete_alignments = 0
        for r in contig_alignments:
            if r.Complete():
                complete_alignments += 1
        if complete_alignments > 1:
            stats.ambiguous.append(barcode)

        assess_alignment(barcode, contig_alignments[0], hagen_dict, stats, unique_assembled_isoforms)

        for r in contig_alignments[1:]:
            assess_alignment(barcode, r, hagen_dict, secondary_stats, unique_assembled_isoforms)

    unique_isoforms = set()     
    for bc in hagen_dict.keys():
        isoform_str = "".join(str(x) for x in hagen_dict[bc])
        unique_isoforms.add(isoform_str)
        if bc not in stats.used:
             stats.matrix_only.append(bc)

    print("Unique assebmled isoformsx " + str(len(unique_assembled_isoforms)))
    print(unique_assembled_isoforms)
    print("Unique isoforms in matrix " + str(len(unique_isoforms)))
    print(unique_isoforms)
    return stats, secondary_stats


if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " <Maxtrix file starting with BIN1 | MAPT> <BAM/SAM file> [Use generated contig->barcodes map: y|n]")
    exit(0)


maxtrix_file = sys.argv[1]
if (maxtrix_file.upper().startswith('BIN1')):
    remapped_coords = [(1160, 1249), (1350, 1457), (1558,1581),(1682,1789),(1890,2132),(3423,3515)]
    initial_coords = [(127051154, 127051243),(127052255,127052362),(127053422, 127053445),(127053905, 127054012),(127057473, 127057601),(127068163, 127068255)]
elif (maxtrix_file.upper().startswith('MAPT')):
    remapped_coords = [(1252, 1338), (1439, 1525), (2987, 3184)]
    #initial_coords = [ (45971859, 45971945), (45974385, 45974471), (45989878, 45990075)]
    initial_coords = [(45969165, 45969299), (45971859, 45971945), (45974385, 45974471), (45989878, 45990075), (45993924,45993977), (46010310, 46010402)]
else:
    print("Gene " + sys.argv[1] + " is not supported")
    exit(1)

samfile_in = pysam.AlignmentFile(sys.argv[2], "rb")

barcode_map = {}
if len(sys.argv) == 3 or (len(sys.argv) == 4 and not sys.argv[3].upper().startswith("N")):
    barcode_map = get_barcode_map(sys.argv[2])


barcode_info, alignments  = process_sam(samfile_in, initial_coords, barcode_map, maxtrix_file)

alignments_set = set(alignments)
print(len(set(alignments)))
print(len(barcode_map.keys()))

count = 0
for bc in barcode_map.keys():
    if bc not in alignments_set:
        count += 1
        #print(bc)
print count


stats, seondary_stats = get_stats(maxtrix_file, initial_coords, barcode_info, barcode_map)
stats.contig_alignment_count = len(alignments)

print("Total contigs: " + str(stats.contig_count) + " alignments: " + str(stats.contig_alignment_count) + ", distinct barcodes " + str(stats.distinct_barcodes_in_contigs) )
print("Total barcodes in matix: " + str(stats.total_matrix_barcodes) )
print("Total barcodes used: " + str(len(stats.used) ))

print("Ambiguous contigs: " + str(len(stats.ambiguous)) + ", distinct barcodes " + str(len(set(stats.ambiguous))))
print("Contigs with secondary alignments: " + str(len(stats.with_secondary)) + ", distinct barcodes " + str(len(set(stats.with_secondary))))

print("== Primary stats ==")
print("Correct contigs: " + str(len(stats.correct)) + ", distinct barcodes " + str(len(set(stats.correct))) )
print("Incorrect contigs: " + str(len(stats.incorrect)) + ", distinct barcodes " + str(len(set(stats.incorrect))))
print("Incorrect1 contigs: " + str(len(stats.incorrect1)) + ", distinct barcodes " + str(len(set(stats.incorrect1))))
print("Incomplete contigs: " + str(len(stats.incomplete)) + ", distinct barcodes " + str(len(set(stats.incomplete))))
print("Incomplete contigs in the matrix: " + str(len(stats.incomplete_in_dict)) + ", distinct barcodes " + str(len(set(stats.incomplete_in_dict))))
print("Incomplete contigs NOT in the matrix: " + str(len(stats.incomplete_missing)) + ", distinct barcodes " + str(len(set(stats.incomplete_missing))))

print("5' " + str(len(stats.incomplete5)) + ", 3': " + str(len(stats.incomplete3)));
print("5' " + str(len(set(stats.incomplete5))) + ", 3': " + str(len(set(stats.incomplete3))));
print("Incomplete corr " + str(len(stats.incomplete_correct)) + ", Incomplete incorr " + str(len(stats.incomplete_incorrect)));
print("Incomplete corr " + str(len(set(stats.incomplete_correct))) + ", Incomplete incorr " + str(len(set(stats.incomplete_incorrect))));

print("Contigs with barcodes that are not in the matrix " + str(len(stats.contigs_only)) + ", distinct barcodes " + str(len(set(stats.contigs_only))))
print("Found only in matrix, but not found in contigs " + str(len(stats.matrix_only)) + ", distinct barcodes " + str(len(set(stats.matrix_only))))


print("== Secondary stats ==")
print("Correct contigs: " + str(len(seondary_stats.correct)) + ", distinct barcodes " + str(len(set(seondary_stats.correct))) )
print("Incorrect contigs: " + str(len(seondary_stats.incorrect)) + ", distinct barcodes " + str(len(set(seondary_stats.incorrect))))
print("Incomplete contigs: " + str(len(seondary_stats.incomplete)) + ", distinct barcodes " + str(len(set(seondary_stats.incomplete))))
print("Incomplete contigs in the matrix: " + str(len(seondary_stats.incomplete_in_dict)) + ", distinct barcodes " + str(len(set(seondary_stats.incomplete_in_dict))))
print("Incomplete contigs NOT in the matrix: " + str(len(seondary_stats.incomplete_missing)) + ", distinct barcodes " + str(len(set(seondary_stats.incomplete_missing))))

print("5' " + str(len(seondary_stats.incomplete5)) + ", 3': " + str(len(seondary_stats.incomplete3)));
print("Incomplete corr " + str(len(seondary_stats.incomplete_correct)) + ", Incomplete incorr " + str(len(seondary_stats.incomplete_incorrect)));


#print("\n".join(stats.matrix_only))

