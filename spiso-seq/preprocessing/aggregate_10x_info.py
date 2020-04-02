import os
import sys
import argparse
import pysam
from traceback import print_exc
import glob


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dirs", nargs='+', help="folders to process", type=str)
    parser.add_argument("--output_file", "-o", help="output file", type=str)

    args = parser.parse_args()
    return args


def flush_read_map(output_file, read_map):
    outf = open(output_file, 'w')
    for read_id in read_map.keys():
        info = read_map[read_id]
        outf.write(read_id + '\t' + info[0] + '\t' + info[1] + '\n')
    outf.close()


def process_bam(bam_file, group_id, barcode_set):
    read_map = {}
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    for alignment in bamfile:
        read_id = alignment.query_name
        barcode = alignment.get_tag('CB').split(":")[-1][:-2]
        if barcode not in barcode_set:
            print("Barcode " + barcode + " is not in set for group " + group_id)
        read_map[read_id] = (barcode, group_id)

    return read_map


def process_dir(dir, output_file):
    tsv_files = glob.glob(os.path.join(dir, '*.tsv'))
    region_name = os.path.basename(dir)
    for f in tsv_files:
        read_group = os.path.splitext(os.path.basename(f))[0]
        print("Processing " + read_group)
        barcode_set = set()
        for l in open(f):
            barcode_set.add(l.strip()[:-2])
        read_map = process_bam(os.path.join(dir, 'OutputBAM/' + region_name + '_' + read_group + '.bam'), read_group, barcode_set)
        flush_read_map(output_file, read_map)
        exit(0)


def main():
    args = parse_args()

    for dir in args.dirs:
        print("Processing " + dir)
        process_dir(dir, args.output_file)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
