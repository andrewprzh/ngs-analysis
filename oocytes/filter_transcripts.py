import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


MOTIF_FWD='ATTCC'
MOTIF_BCK='GGAAT'


def check_sequence(seq, max_len, motif_percentage):
    seq_len = len(seq)
    if seq_len > max_len:
        return False, False
    fwd_count = seq.count(MOTIF_FWD)
    bck_count = seq.count(MOTIF_BCK)
    return max(fwd_count, bck_count) * 5.0 >= motif_percentage * seq_len, fwd_count >= bck_count


def filter_fasta(infasta, outfasta, max_len, motif_percentage):
    all_datasets = set()
    filtered_datasets = set()
    with open(outfasta, "w") as output_handle:
        new_contigs = []
        for seq_record in SeqIO.parse(infasta, "fasta"):
            dataset_id = seq_record.id.split('_')[-1]
            all_datasets.add(dataset_id)
            passed, rc = check_sequence(seq_record.seq, max_len, motif_percentage)
            if not passed:
                print("Contig " + seq_record.id + " seems unreliable")
                continue
            if rc:
                seq_record.id += "_RC"
                seq_record.description = ""
                seq_record.seq = seq_record.seq.reverse_complement()
            new_contigs.append(seq_record)
            filtered_datasets.add(dataset_id)

        SeqIO.write(new_contigs, output_handle, "fasta")
    print(all_datasets)
    print("Total datasets: "+str(len(all_datasets))+", fileterd: "+str(len(filtered_datasets))+", missed: "+str(all_datasets-filtered_datasets))



filter_fasta(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
