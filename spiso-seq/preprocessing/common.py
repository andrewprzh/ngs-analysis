############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from Bio import SeqIO


def get_barcodes(contig_id, delim1, delim2):
    tokens = contig_id.split(delim1)
    if len(tokens) != 2:
        print("Wrong fromat " + contig_id)

    return tokens[0], tokens[1].strip().split(delim2) if len(tokens) > 1 else []


def convert_fasta_with_barcodes(inf, args):
    contigs_file = inf
    contigs_name, ext = os.path.splitext(contigs_file.split('/')[-1])
    short_id_contigs_name = os.path.join(args.output_prefix, contigs_name + "_short_ids.fasta")
    map_file_name = os.path.join(args.output_prefix, contigs_name + "_map.txt")
    mapf = open(map_file_name, "w")

    new_fasta = []
    for record in SeqIO.parse(contigs_file, "fasta"):
        if len(record.seq) > args.max_len:
            continue
        name = record.id
        num = int(name.split('_')[1])
        if num % 2 == 1:
            continue
        record.id, bc = get_barcodes(name, args.delim, args.delim2)
        new_fasta.append(record)
        mapf.write(record.id + "_barcodeIDs_" + ",".join(bc) + "\n")

    mapf.close()
    SeqIO.write(new_fasta, short_id_contigs_name, "fasta")
    return contigs_name, short_id_contigs_name



def align_fasta(contigs_infile, args, out_name = ""):
    method = args.aligner
    contigs_name, ext = os.path.splitext(contigs_infile.split('/')[-1])
    if out_name == "":
        out_name = contigs_name

    alignment_name = os.path.join(args.output_prefix, out_name)
    alignment_sam_path = alignment_name + '.sam'
    alignment_bam_path = alignment_name + '.bam'

    if method.lower() == "gmap":
        gmap_path = '/Bmo/prjbel/tools/gmap-2019-06-10/src/gmap'
        command = gmap_path + ' -D ./  -d {ref_index_name} {transcripts} --format=samse -t {threads} -O ' \
                              '> {alignment_out} '.format(ref_index_name=args.index,
                                                          transcripts=contigs_infile,
                                                          threads=str(args.threads),
                                                          alignment_out=alignment_sam_path)
        exit_code = os.system(command)

        if exit_code != 0:
            print("GMAP finished with errors")

        os.system('samtools sort -o ' + alignment_bam_path + ' ' + alignment_sam_path)

    elif method.lower() == "star":
        star_path = '/Bmo/prjbel/tools/STAR-2.7.0f/bin/Linux_x86_64/STARlong'
        zcat_option = " --readFilesCommand zcat " if ext.endswith('gz') else ""

        # Simple
        # command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM
        #  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        # Hagen options
        command = '{star} {zcat} --runThreadN {threads} --genomeDir {ref_index_name}  --readFilesIn {transcripts}  ' \
                  '--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 \
                   --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen ' \
                  '-1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 \
                   --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 1000000 --seedPerWindowNmax 1000 ' \
                  '--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 \
                   --outFileNamePrefix {alignment_out}'.format(star=star_path,
                                                               zcat=zcat_option,
                                                               threads=str(args.threads),
                                                               ref_index_name=args.index,
                                                               transcripts=contigs_infile,
                                                               alignment_out=alignment_name)

        exit_code = os.system(command)
        os.system('mv ' + alignment_name + 'Aligned.out.sam ' + alignment_bam_path)

        if exit_code != 0:
            print("STAR finished with errors")

    elif method.lower == "minimap":
        minimap_path = '/Bmo/prjbel/tools/minimap2-2.8_x64-linux/minimap2'
        command =  minimap_path+ ' {ref_index_name} {transcripts} -a -x splice -t {threads} ' \
                              '> {alignment_out} '.format(ref_index_name=args.index,
                                                          transcripts=contigs_infile,
                                                          threads=str(args.threads),
                                                          alignment_out=alignment_sam_path)
        exit_code = os.system(command)
        if exit_code != 0:
            print("minimap finished with errors")

    else:
        print("Method " + method + " is not supported")
        return 

    os.system('samtools index ' + alignment_bam_path)

