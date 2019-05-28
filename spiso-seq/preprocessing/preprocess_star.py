import os
import sys

def read_fasta(fpath):
    """
        Returns list of FASTA entries (in tuples: name, seq)
    """
    first = True
    seq = ''
    name = ''

    fasta_file = open(fpath)

    for raw_line in fasta_file:
        if raw_line.find('\r') != -1:
            lines = raw_line.split('\r')
        else:
            lines = [raw_line]
        for line in lines:
            if not line:
                continue
            if line[0] == '>':
                if not first:
                    yield name, seq

                first = False
                name = line.strip()[1:]
                seq = ''
            else:
                seq += line.strip()

    if name or seq:
        yield name, seq

    fasta_file.close()


def print_fasta(fasta):
    for name, seq in fasta:
        print ('>%s' % name)
        for i in xrange(0, len(seq), 60):
            print (seq[i:i + 60])


def write_fasta(fpath, fasta, mode='w'):
    outfile = open(fpath, mode)

    for name, seq in fasta:
        outfile.write('>%s\n' % name)
        for i in xrange(0, len(seq), 60):
            outfile.write(seq[i:i + 60] + '\n')
    outfile.close()



if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <FASTA file> <STAR index> <output dir>")
    exit(0)

contigs_file = sys.argv[1]
star_index = sys.argv[2]
output_dir = sys.argv[3]

contigs_name, ext = os.path.splitext(contigs_file.split('/')[-1])
short_id_contigs_name = os.path.join(output_dir, contigs_name + "_short_ids.fasta")
map_file_name = os.path.join(output_dir, contigs_name + "_map.txt")
mapf = open(map_file_name, "w")

contigs = read_fasta(contigs_file)
new_fasta = []
for name,seq in contigs:
    tokens = name.split("_barcodeIDs_")
    if len(tokens) != 2:
        print("Wrong fromat " + name)
    contig_id = tokens[0]
    new_fasta.append((contig_id, seq))
    mapf.write(contig_id + "_barcodeIDs_" + tokens[1].strip() + "\n")

mapf.close()
write_fasta(short_id_contigs_name, new_fasta)

alignment_sam_path = os.path.join(output_dir, contigs_name)
star_path = '/Bmo/prjbel/tools/STAR-2.7.0f/bin/Linux_x86_64/STARlong'
#Simple 
#command = '{star} --runThreadN 8 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM  --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
#Hagen options
command = '{star} --runThreadN 8 --genomeDir {ref_index_name}  --readFilesIn {transcripts}  --outSAMtype SAM --outSAMattributes NH HI NM MD --outFilterMultimapScoreRange 1 \
           --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 \
           --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 \
           --outFileNamePrefix {alignment_out}'.format(star=star_path, ref_index_name=star_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)

exit_code = os.system(command)
os.system('mv ' + alignment_sam_path + 'Aligned.out.sam ' + alignment_sam_path + '.sam')

if exit_code != 0:
    print("STAR finished with errors")


