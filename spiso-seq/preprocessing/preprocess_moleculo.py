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
    print("Usage: " + sys.argv[0] + " <FASTA file> <GMAP index> <output dir>")
    exit(0)

contigs_file = sys.argv[1]
gmap_index = sys.argv[2]
output_dir = sys.argv[3]

contigs_name, ext = os.path.splitext(contigs_file.split('/')[-1])
short_id_contigs_name = os.path.join(output_dir, contigs_name + "_short_ids.fasta")
map_file_name = os.path.join(output_dir, contigs_name + "_map.txt")
mapf = open(map_file_name, "w")

contigs = read_fasta(contigs_file)
new_fasta = []
for name,seq in contigs:
    tokens = name.split("-Barcode=")
    if len(tokens) != 2:
        print("Wrong fromat " + name)
    contig_id = tokens[0]
    bc =  tokens[1].strip().split(':')[0] if len(tokens) > 1 else ""
    new_fasta.append((contig_id, seq))
    mapf.write(contig_id + "_barcodeIDs_" + bc + "\n")

mapf.close()
write_fasta(short_id_contigs_name, new_fasta)

alignment_sam_path = os.path.join(output_dir, contigs_name + '.sam')
command = 'gmap -D ./  -d {ref_index_name} {transcripts} --format=samse -t 16 -O > {alignment_out} '.format(ref_index_name=gmap_index, transcripts=short_id_contigs_name, alignment_out=alignment_sam_path)
exit_code = os.system(command)

if exit_code != 0:
    print("GMAP finished with errors")


