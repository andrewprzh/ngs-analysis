############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
import gzip
from Bio import SeqIO
from traceback import print_exc


BARCODE_LENGTHS = [16, 14]

def read_barcode_map(inf):
    barcode_map = {}
    for l in open(inf):
        tokens = l.strip().split()
        if len(tokens) != 5:
            continue
        barcode_map[tokens[0].split('_')[0]] = tokens[4]

    return barcode_map


def process_single_read_sample(sample_id, barcode_file, read_file, barcode_map, outf_map):
    if barcode_file.endswith("gz"):
        bc_handle = SeqIO.parse(gzip.open(barcode_file, "rt"), "fastq")
    else:
        bc_handle = SeqIO.parse(open(barcode_file, "rt"), "fastq")
    if read_file.endswith("gz"):
        read_handle = SeqIO.parse(gzip.open(read_file, "rt"), "fastq")
    else:
        read_handle = SeqIO.parse(open(read_file, "rt"), "fastq")

    for bc_record in bc_handle:
        read_rec = next(read_handle)
        if read_rec.id.split()[0] != bc_record.id.split()[0]:
            sys.stderr.write("Unequal read ids " + read_rec.id.split()[0] + ",  " + bc_record.id.split()[0] + "\n")
            continue

        for l in BARCODE_LENGTHS:
            barcode = str(bc_record.seq[:l])
            if barcode in barcode_map:
                cluster_id = barcode_map[barcode]
                read_rec.description += " sample=" + sample_id + " barcode=" + barcode + " cluster=" + cluster_id
                outf = outf_map[cluster_id]
                SeqIO.write(read_rec, outf, "fastq")
                break


def process_paired_read_sample(sample_id, barcode_file, read_file1, read_file2, barcode_map, outf_map1, outf_map2):
    if barcode_file.endswith("gz"):
        bc_handle = SeqIO.parse(gzip.open(barcode_file, "rt"), "fastq")
    else:
        bc_handle = SeqIO.parse(open(barcode_file, "rt"), "fastq")
    if read_file1.endswith("gz"):
        read_handle1 = SeqIO.parse(gzip.open(read_file1, "rt"), "fastq")
    else:
        read_handle1 = SeqIO.parse(open(read_file1, "rt"), "fastq")
    if read_file2.endswith("gz"):
        read_handle2 = SeqIO.parse(gzip.open(read_file2, "rt"), "fastq")
    else:
        read_handle2 = SeqIO.parse(open(read_file2, "rt"), "fastq")

    for bc_record in bc_handle:
        read_rec1 = next(read_handle1)
        read_rec2 = next(read_handle2)
        if read_rec1.id.split()[0] != bc_record.id.split()[0] or read_rec2.id.split()[0] != bc_record.id.split()[0]:
            sys.stderr.write("Unequal read ids " + read_rec1.id.split()[0] + ",  " + read_rec1.id.split()[0] + ",  " + bc_record.id.split()[0] + "\n")
            continue

        for l in BARCODE_LENGTHS:
            barcode = str(bc_record.seq[:l])
            if barcode in barcode_map:
                cluster_id = barcode_map[barcode]
                read_rec1.description += " sample=" + sample_id + " barcode=" + barcode + " cluster=" + cluster_id
                read_rec2.description += " sample=" + sample_id + " barcode=" + barcode + " cluster=" + cluster_id
                outf1 = outf_map1[cluster_id]
                outf2 = outf_map2[cluster_id]
                SeqIO.write(read_rec1, outf1, "fastq")
                SeqIO.write(read_rec2, outf2, "fastq")
                break


def get_sample_name(fname):
    for r in ['R1', 'R2', 'R3']:
        if fname.find(r) != -1:
            return fname.split("/")[-1].split(r)[0][:-1]

    return fname.split('.')[0]


def read_samples(sample_file):
    sample_tuples = []
    for l in open(sample_file, "rt"):
        if l.strip() == "":
            continue

        tokens = l.strip().split()
        if len(tokens) == 2:
            sample_tuples.append((get_sample_name(tokens[0]), tokens[0], tokens[1]))
        elif len(tokens) == 3:
            sample_tuples.append((get_sample_name(tokens[0]), tokens[0], tokens[1], tokens[2]))
        else:
            sys.stderr.write("Malformed input, line " + l + "\n")
    return sample_tuples


def create_output_file_map(barcode_map, output_prefix):
    name_map = {}
    cluster_ids = set(barcode_map.values())
    for cluster_id in cluster_ids:
        name_map[cluster_id] = output_prefix + ".cluster_" + cluster_id

    read_type_suffix = { "single" : ".S.fastq", "left" : ".R1.fastq", "right" : ".R2.fastq"}
    outf_maps = {}
    for read_type in read_type_suffix.keys():
        outf_maps[read_type] = {}
        for cluster_id in name_map.keys():
            outf_maps[read_type][cluster_id] = open(name_map[cluster_id] + read_type_suffix[read_type], "wt")

    return outf_maps


def close_files(outf_maps):
    for read_type in outf_maps.keys():
        for cluster_id in outf_maps[read_type]:
            outf_maps[read_type][cluster_id].close()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", "-s", help="files with tab-separated sample info, each line has barcoded file and file with reads", type=str, default = "10x")
    parser.add_argument("--barcode_map", "-b", help="tsne.data extracted table with barcodes and cluster numbers", type=str)
    parser.add_argument("--output_prefix", "-o", help="output prefix", type=str)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    barcode_map = read_barcode_map(args.barcode_map)
    print("Loaded " + str(len(barcode_map.keys())) + " barcodes with " + str(len(set(barcode_map.values()))) + " clusters")
    samples = read_samples(args.samples)
    outf_maps = create_output_file_map(barcode_map, args.output_prefix)

    for sample in samples:
        print("Processing sample " + sample[0])
        if len(sample) == 3:
            process_single_read_sample(sample[0], sample[1], sample[2], barcode_map, outf_maps["single"])
        elif len(sample) == 4:
            process_paired_read_sample(sample[0], sample[1], sample[2], sample[3], barcode_map, outf_maps["left"], outf_maps["right"])

    close_files(outf_maps)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



