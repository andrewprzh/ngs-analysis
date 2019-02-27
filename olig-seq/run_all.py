#!/usr/bin/python

import os
import sys
 
if len(sys.argv) < 6:
    print("Usage: " + sys.argv[0] + " <oligs in FASTA> <file list> <output folder> <full | restricted> <output stat file> ")
    exit(0)


DB = sys.argv[1]
FILES = []
OUT_DIR = sys.argv[3]

restricted = sys.argv[4].startswith('r')
stat_file = sys.argv[5]

for l in open(sys.argv[2]):
    files = l.strip().split()
    if len(files) != 2:
        print("Wrong input list, must contain 2 files per line, ignoring entry")
        continue
    FILES.append((files[0], files[1]))

if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)
index_dir = os.path.join(OUT_DIR, "index")
trimmed_dir = os.path.join(OUT_DIR, "trimmed")
if not os.path.exists(index_dir):
    os.mkdir(index_dir)
if not os.path.exists(trimmed_dir):
    os.mkdir(trimmed_dir)

statf = open(stat_file, "w")
statf.write("Sample\tInput read pairs\tRead pairs with adapters\tMapped\tOligs covered\n")
statf.close()

for f in FILES:
    fbase = f[0].split('/')[-1].split('R1')[0]
    log_f = os.path.join(OUT_DIR, fbase + ".log")
    print(fbase + "\n")
    statf = open(stat_file, "a+")
    statf.write(fbase + "\t")
    statf.close()
    trimmed_l = os.path.join(trimmed_dir, fbase + "R1.trimmed.fastq")
    trimmed_r = os.path.join(trimmed_dir, fbase + "R2.trimmed.fastq")

    print("   Running cutadapt")
    if not restricted:
        os.system("cutadapt -e 0.25 -g CCCAGCCGGCCATGGCC -a GCTAGCAGTGGTGGAGGCGG -g CCGCCTCCACCACTGCTAGC -a GGCCATGGCCGGCTGGG  -G CCCAGCCGGCCATGGCC -A GCTAGCAGTGGTGGAGGCGG -G  CCGCCTCCACCACTGCTAGC -A GGCCATGGCCGGCTGGG  -o " + trimmed_l + " -p " + trimmed_r + " " + f[0] + " " + f[1] + " >> " + log_f)
    else:
        os.system("cutadapt -e 0.1 -g CATGGCC -a GCTAG -g CTAGC -a GGCCATG  -G CATGGCC -A GCTAG -G  CTAGC -A GGCCATG  -o " + trimmed_l + " -p " + trimmed_r + " " + f[0] + " " + f[1] + " >> " + log_f)    

    merged_f = os.path.join(OUT_DIR, fbase + "merged.fastq")
    merge_stat_f = os.path.join(OUT_DIR, fbase + "merge_stats.tsv")

    print("   Merging reads")
    os.system("python2 merge_and_filter.py "  + f[0] + " " + f[1] + " " + trimmed_l + " " + trimmed_r + " " + merged_f + " FASTQ " + ("r" if restricted else "f") + " " + stat_file + " > " + merge_stat_f)
    map_stat_f = os.path.join(OUT_DIR, fbase + "map_stats.tsv")
    sam_f = os.path.join(OUT_DIR, fbase + ".sam")

    print("   Running minimap")
    os.system("./minimap2-2.14_x64-linux/minimap2 " + DB + " " + merged_f + " -a -x sr >  " + sam_f + " 2> " + os.path.join(OUT_DIR, fbase + "minimap2.log"))
    merged_fasta = os.path.join(OUT_DIR, fbase + "merged.fasta")

    print("   Converting merged reads to FASTA")
    os.system("seqtk seq -a " + merged_f + " > " + merged_fasta)

    print("   Mapping merged reads to olig DB " + DB)
    os.system("python2 map_all_with_sam.py " + DB + " " + merged_fasta + " " + sam_f + " " + stat_file + " > " + map_stat_f)


