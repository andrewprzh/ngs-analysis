############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc


def read_abundances(inf):
    values = []
    for l in open(inf):
        t = l.strip().split()
        values.append(int(t[0]))
    return values


def read_genomes(inf):
    genomes = []
    for l in open(inf):
        t = l.strip().split()
        genomes.append(t[0])
    return genomes


def simulate_metag(genomes, counts, output_prefix):
    assert len(genomes) == len(counts)
    tmp_output_fastas = []
    for i, genome in enumerate(genomes):
        tmp_out = output_prefix + "_" + str(i) + ".fastq"
        tmp_output_fastas.append(tmp_out)
        print("Simulating %d reads for %s" % (counts[i], genome))
        os.system('~/.local/bin/iss generate -p 16  --genomes ' + genome + ' --n_reads ' + str(counts[i]) + ' --model HiSeq  -a uniform -o ' + tmp_out + ' 2> /dev/null')

    print("Combining reads from %d files" % len(tmp_output_fastas))
    main_output = output_prefix + ".fastq"
    open(main_output, "w").close()
    for tmp_out in tmp_output_fastas:
        os.system('cat ' + tmp_out + ' >> ' + main_output)
        os.remove(tmp_out)

    print("All reads written to " + main_output)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", "-g", help="initial file with refrenece sequences", type=str, required=True)
    parser.add_argument("--abundances", "-a", help="abundances in the same order as genomes", type=str, required=True)
    parser.add_argument("--output", "-o", help="output folder", type=str, required=True)
    parser.add_argument("--count", "-c", help="total reads to generate", type=int, default=5000000)

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    return args


def main():
    args = parse_args()
    genomes = read_genomes(args.gnomes)
    abundances = read_abundances(args.abundances)
    total_abundance = sum(abundances)
    read_counts = [args.count * abundance / total_abundance for abundance in abundances]
    simulate_metag(genomes, read_counts, args.output + "")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
