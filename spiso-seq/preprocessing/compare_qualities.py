import pysam
import sys
import argparse
import os
from traceback import print_exc


def process_bam(bam, split_id=False):
    q_dict = {}
    for a in pysam.AlignmentFile(bam, "rb"):
        if not a.is_primary:
            continue
        read_id = a.query_name
        if split_id:
            read_id = read_id.split("_")[0]
        mapq = a.mapping_quality
        try:
            dpas = a.get_tag("AS")
        except KeyError:
            dpas = -1
        q_dict[read_id] = [mapq, dpas]
    return q_dict


def process_bam_pair(raw_bam, corr_bam):
    raw_qdict = process_bam(raw_bam)
    corr_qdict = process_bam(corr_bam, split_id=True)
    corrected_reads = 0
    dropped_reads = 0
    for read_id in corr_qdict:
        if read_id in raw_qdict:
            corr_qdict[read_id] += raw_qdict[read_id]
            corrected_reads +=1
        else:
            dropped_reads += 1
    return corr_qdict, corrected_reads, dropped_reads


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rawbam", type=str, help="initial bam")
    parser.add_argument("--corrbam", type=str, help="corrected bam")
    parser.add_argument("--list", "-l", type=str, help="bam files list, 2 per row")
    parser.add_argument("--output", "-o", type=str, help="output folder")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    file_pairs = []
    if args.list:
        for l in open(args.list):
            t = l.strip().split()
            if len(t) == 2:
                file_pairs.append((t[0], t[1]))
    else:
        file_pairs = [(args.rawbam, args.corrbam)]

    for pair in file_pairs:
        print("Processing pair[0]")
        qdict, corrected_reads, dropped_reads = process_bam_pair(pair[0], pair[1])
        print("Corrected reads %d, dropped %d" % (corrected_reads, dropped_reads))

        with open(os.path.join(args.output, os.path.basename(pair[0]) + ".qualities.tsv")) as outf:
            for read_id in qdict:
                outf.write("\t".join([read_id] + list(map(str, qdict[read_id]))) + "\n")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

