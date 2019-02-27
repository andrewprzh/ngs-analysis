#!/usr/bin/python

import os
import sys
import pysam
import common
from common import *

dbg_barcodes = set()

#['BX:Z:AGCTGTTCGGGAGT-1', 'BX:Z:TAGACCTCGTTGGT-1', 'BX:Z:AGTAAGTCTCCTAT-1', 'BX:Z:CTGCAATCTCAAGC-1', 'BX:Z:GTACCAAGCCCTGC-1', 'BX:Z:GTAGTCAGACTATG-1', 'BX:Z:GTGGAGGTAGTCGG-1', 'BX:Z:TTCCCGCATCAGGG-1', 'BX:Z:AATAGGTCTATCCT-1', 'BX:Z:CACCTCCATCGTGT-1', 'BX:Z:GTCGATCAAGCTTC-1', 'BX:Z:AGACCTAGTTCTCG-1', 'BX:Z:GTAGAATCCTATGG-1', 'BX:Z:CCTAGATCGTCCTC-1', 'BX:Z:CGATTGCAGGTCGA-1', 'BX:Z:GGTACTCAGTCCGC-1', 'BX:Z:GTGAGGTCCTCATT-1'])
#['BX:Z:ACGGGAAGGCCATG-1'])
#set(['BX:Z:TAGACCTCGTTGGT-1'])

def debug(bc, s):
    if bc in dbg_barcodes:
        print(s)


def process_read(alignment, exons):
    tokens = alignment.query_name.strip().split("___")
    if len(tokens) != 2:
        return None, None

    barcode = tokens[1]
    exons_present = [0 for i in range(2 + len(exons))]

    blocks = sorted(alignment.get_blocks())
    block_pos = 0
    exon_pos = 0
    
    #if barcode in dbg_barcodes:
    #    print blocks
    #    print exons

    
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

    if left_of(blocks[0], exons[0]):
        exons_present[0] = 1

    if left_of(exons[-1], blocks[-1]):
        exons_present[-1] = 1

    rr = ReadRecord(blocks, exons_present)

    #debug(barcode, str(barcode) + ": " + str(exons_present))
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


def process_sam(samfile_in, exon_coords):
    barcode_info = {}
    count = 0
    for r in samfile_in:
        count += 1
        if count % 10000 == 0: 
            print count
        barcode, read_record = process_read(r, exon_coords)
        if barcode is None:
            continue
        if barcode not in barcode_info:
            barcode_info[barcode] = BarcodeRecord(barcode, [read_record])
        else:
            barcode_info[barcode].AddAlignment(read_record)

    return barcode_info


def merge_barcodes(barcode_info):
    new_bc_info = {}
    single_am_count = 0
    multiple_am_count = 0

    for barcode, bc_records in barcode_info.iteritems():
        exon_scores = [{-1:0,0:0,1:0} for i in range(0, len(bc_records.read_infos[0].exon_present))]
        
        coverage = 0
        for r in bc_records.read_infos:
            coverage += 1
            for i in range(0, len(r.exon_present)):
                exon_scores[i][r.exon_present[i]] += 1
        #debug(barcode, barcode + ": " + str(exon_scores))

        support = 1000000
        consensus_exons = [0 for i in range(0, len(exon_scores))]
        ambigous = []
        for i in range(0, len(exon_scores)):
            if (exon_scores[i][-1] == 0 and exon_scores[i][1] > 0) or (exon_scores[i][-1] * RELATIVE_COUNT_COEFF < exon_scores[i][1]):
                if barcode in dbg_barcodes and exon_scores[i][-1] != 0:
                    print("Ambiguous at " + str(i))
                    print(barcode + ": " + str(exon_scores))
                consensus_exons[i] = 1
                support = min(support, exon_scores[i][1])
            elif (exon_scores[i][1] == 0 and exon_scores[i][-1] > 0) or (exon_scores[i][1] * RELATIVE_COUNT_COEFF < exon_scores[i][-1]):
                if barcode in dbg_barcodes and exon_scores[i][1] != 0:
                    print("Ambiguous at " + str(i))
                    print(barcode + ": " + str(exon_scores))
                consensus_exons[i] = -1
                support = min(support, exon_scores[i][-1])
            elif exon_scores[i][-1] >= 1 and exon_scores[i][1] >= 1:
                #ambiguous exons
                if barcode in dbg_barcodes:
                    print("Ambiguous at " + str(i))
                    print(barcode + ": " + str(exon_scores))
                ambigous.append(i)

        #debug(barcode, barcode + ": " + str(ambigous) + " : " + str(exon_scores))

        if len(ambigous) == 0:
#            for r in bc_records.read_infos:
#                if diff_only_present(r.exon_present, consensus_exons) == 0:

            rr = ReadRecord([], consensus_exons, support, coverage)
            new_bc_info[barcode] = BarcodeRecord(barcode, [rr])
        elif len(ambigous) == 1:
            print("Single ambiguity: " + str(ambigous))
            print(barcode + ": " + str(exon_scores))
            single_am_count += 1
        else:
            multiple_am_count += 1

    print ("Single exon: " + str(single_am_count))
    print ("Multiple exons: " + str(multiple_am_count))
    return new_bc_info


def filter_barcodes(barcode_info, include_incomplete = False):
    new_bc_info = {}
    stats = BarcodeStats()
    for barcode, bc_records in barcode_info.iteritems():
        debug(barcode, str(barcode) + ": " + str(bc_records.read_infos[0].exon_present))
        stats.contig_count += 1
        ri = []
        low_cov_ri = []
        for r in bc_records.read_infos:
            if r.Supported():
                debug(barcode, "Supported")
                ri.append(r)
            else:
                debug(barcode, "Not supported")
                low_cov_ri.append(r)

        if len(ri) == 0 and len(low_cov_ri) > 0:
            stats.low_covered.append(barcode)
            if len(low_cov_ri) == 1 and low_cov_ri[0].Complete():
                stats.missed.append(barcode)

        if len(ri) > 1:
            print(barcode + " has multiple isoforms")
            stats.ambiguous.append(barcode)
            for r in ri:
                print(str(r.exon_present) + "; " + str(r.block_coordinates) + "; " + str(r.support))

            ri_sorted = sorted(ri, key = lambda x: x.support, reverse=True)
            if ri_sorted[0].Complete() and not ri_sorted[1].Complete() and ri_sorted[0].support / ri_sorted[1].support >= RELATIVE_COUNT_COEFF:
                print("Resolved")
                new_bc_info[barcode] = BarcodeRecord(barcode, [ri_sorted[0]])
                stats.incomplete_correct.append(barcode)
                debug(barcode, "Resolved")

        if len(ri) == 1:
            if ri[0].Complete():
                new_bc_info[barcode] = BarcodeRecord(barcode, ri)
                stats.correct.append(barcode)
                debug(barcode, "Complete")
            else:
                if reduce(lambda x, y: x and y, map(lambda x: x == 0, ri[0].exon_present[1:-1])):
                    stats.incorrect.append(barcode)
                    continue

                stats.incomplete.append(barcode)
                debug(barcode, "Incomplete")

                if include_incomplete:
                    new_bc_info[barcode] = BarcodeRecord(barcode, ri)

                if ri[0].exon_present[1] == 0:
                    stats.incomplete5.append(barcode)
                if ri[0].exon_present[-2] == 0:
                    stats.incomplete3.append(barcode)
                
                if ri[0].exon_present[1] != 0 and ri[0].exon_present[-2] != 0:
                    stats.incomplete_middle.append(barcode)

    return new_bc_info, stats



def print_barcodes(barcode_info):
    for barcode, bc_records in barcode_info.iteritems():
        print(barcode)
        for r in bc_records.read_infos:
            print(str(r.exon_present) + "; " + str(r.block_coordinates) + "; " + str(r.support))
        raw_input("Press Enter to continue...")


def print_maxtrix(barcode_info, file_name):
    print("Printing " + str(len(barcode_info)) + " barcodes")
    outf = open(file_name, "w")
    for barcode, br in barcode_info.iteritems():
        if len(br.read_infos) == 0:
            print("Empty exon info")
        #print(br.read_infos[0].exon_present)
        outf.write(barcode + "\t" + "\t".join(map(lambda x: str(x), br.read_infos[0].exon_present[1:-1])) + "\n")
    outf.close()


if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " [BIN1 | MAPT] <SAM file> <output matrix file> [File with debug barcodes]")
    exit(0)


if len(sys.argv) > 4:
    bcf = open(sys.argv[4])
    for l in bcf:
        dbg_barcodes.add(l.strip())


if (sys.argv[1] == 'BIN1'):
    remapped_coords = [(1160, 1249), (1350, 1457), (1558,1581),(1682,1789),(1890,2132),(3423,3515)]
    initial_coords = [(127051154, 127051243),(127052255,127052362),(127053422, 127053445),(127053905, 127054012),(127057359, 127057601),(127068163, 127068255)]
    maxtrix_file = "BIN1.txt"
    exons_number = 6
elif  (sys.argv[1] == 'MAPT'):
    remapped_coords = [(1252, 1338), (1439, 1525), (2987, 3184)]
    initial_coords = [(45969165, 45969299), (45971859, 45971945), (45974385, 45974471), (45989878, 45990075), (45993924,45993977), (46010310, 46010402)]
    maxtrix_file = "MAPT.txt"
    exons_number = 6
else:
    print("Gene " + sys.argv[1] + " is not supported")
    exit(1)

samfile_in = pysam.AlignmentFile(sys.argv[2], mode = "r")
barcode_info = process_sam(samfile_in, initial_coords)
print("Total barcodes read from SAM file " + str(len(barcode_info)))
bc_info, stats = filter_barcodes(merge_barcodes(barcode_info), False)
print_maxtrix(bc_info, sys.argv[3])
print("Barcodes processed ", stats.contig_count)
print("Low covered ", len(stats.low_covered))
print("Low covered but complete ", len(stats.missed))
print("Ambiguous", len(stats.ambiguous), " of them resolved ", len(stats.incomplete_correct) )
print("Complete ", len(stats.correct))
print("Empty ", len(stats.incorrect))
print("Incomplete ", len(stats.incomplete))
print("Incomplete start ", len(stats.incomplete5))
print("Incomplete end ", len(stats.incomplete3))
print("Incomplete middle ", len(stats.incomplete_middle))





