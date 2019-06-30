from Bio import SeqIO


def get_spades_barcode(contig_id):
    delimeter = "_barcodeIDs_"
    tokens = contig_id.split(delimeter)
    if len(tokens) != 2:
        print("Wrong fromat " + contig_id)

    return tokens[0], tokens[1].strip() if len(tokens) > 1 else ""
    

def convert_fasta_with_barcodes(contigs_file, output_dir, barcode_from_id_function):
    contigs_name, ext = os.path.splitext(contigs_file.split('/')[-1])
    short_id_contigs_name = os.path.join(output_dir, contigs_name + "_short_ids.fasta")
    map_file_name = os.path.join(output_dir, contigs_name + "_map.txt")
    mapf = open(map_file_name, "w")

    contigs = read_fasta(contigs_file)
    new_fasta = []
    for record in SeqIO.index(contigs_file, "fasta"):
        record.id, bc = barcode_from_id_function(name)
        new_fasta.append(record)
        mapf.write(record.id + "_barcodeIDs_" + bc + "\n")

    mapf.close()
    SeqIO.write(new_fasta, short_id_contigs_name, "fasta")
    return contigs_name, short_id_contigs_name


def align_fasta(contigs_name, short_id_contigs_name, index, method = "gmap"):
    alignment_name = os.path.join(output_dir, contigs_name)
    alignment_sam_path = alignment_name + '.sam'
    if method.lower() == "gmap":
        gmap_path = '/Bmo/prjbel/tools/gmap-2019-06-10/src/gmap'
        command = gmap_path + ' -D ./  -d {ref_index_name} {transcripts} --format=samse -t 16 -O > {alignment_out} '.format(ref_index_name=index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        exit_code = os.system(command)

        if exit_code != 0:
            print("GMAP finished with errors")

    elif method.lower() == "star":
        star_path = '/Bmo/prjbel/tools/STAR-2.7.0f/bin/Linux_x86_64/STARlong'
        #Simple 
        #command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
        #Hagen options
        command = '{star} --runThreadN 16 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 \
                   --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 \
                   --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 \
                   --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=index, transcripts=short_id_contigs_name, alignment_out=alignment_name)

        exit_code = os.system(command)
        os.system('mv ' + alignment_name + 'Aligned.out.sam ' + alignment_sam_path)

        if exit_code != 0:
            print("STAR finished with errors")

    else:
        print("Method " + method + " is not supported")
        return 

    alignment_bam_path = os.path.join(output_dir, contigs_name + '.bam')
    os.system('samtools sort -o ' + alignment_bam_path + ' ' + alignment_sam_path)
    os.system('samtools index ' + alignment_bam_path)

