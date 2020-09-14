import sys
import os
from Bio import SeqIO
from collections import defaultdict


def get_fasta_ids(fasta):
    ids = set()
    for seq_record in SeqIO.parse(fasta, "fasta"):
        ids.add(seq_record.id)
    return ids


def create_dataset_dict(inf):
    contig_dict = defaultdict(str)
    for l in open(inf):
        l = l.strip()
        if not l or l.find("assembly_SRR") == -1:
            continue
        dataset_name = l[l.find("assembly_SRR") + 9:]
        fasta = get_fasta_ids(os.path.join(l, "transcripts.fasta"))
        for id in fasta:
            if id in contig_dict:
                print("Duplicate id: " + id)
                continue
            contig_dict[id] = dataset_name
    return contig_dict


def process_fasta(infasta, contig_dict, outfasta):
    with open(outfasta, "w") as output_handle:
        new_contigs = []
        for seq_record in SeqIO.parse(infasta, "fasta"):
            id = seq_record.id
            if id not in contig_dict:
                print("Not found: " + id)
            else:
                print("Contig " + id + " will be assigned to " + contig_dict[id])
                id += "_" + contig_dict[id]
                seq_record.id = id
                seq_record.description=""
            new_contigs.append(seq_record)

        SeqIO.write(new_contigs, output_handle, "fasta")


contig_dict = create_dataset_dict(sys.argv[2])
process_fasta(sys.argv[1], contig_dict, sys.argv[3])
