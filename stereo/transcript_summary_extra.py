import argparse
import gzip
import sys
from collections import defaultdict, Counter

def open_maybe_gzip(file_path):
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')

def load_transcript_ids(file_path):
    with open(file_path) as f:
        return set(line.strip() for line in f if line.strip())

def parse_read_files(files, valid_transcripts):
    read_info = dict()
    transcript_to_reads = defaultdict(lambda: defaultdict(list))  # transcript_id -> sample -> list of read_ids

    for sample_idx, file_path in enumerate(files):
        with open_maybe_gzip(file_path) as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 7:
                    continue
                read_id, transcript_id, assign_type, splicing = cols[0], cols[3], cols[5], cols[6]
                splicing_events = list(map(lambda x: x.split(':')[0], splicing.split(','))) if splicing != '.' else []

                read_info[read_id] = {
                    'transcript_id': transcript_id,
                    'assignment_type': assign_type,
                    'splicing_events': splicing_events
                }

                if transcript_id in valid_transcripts:
                    is_unique = assign_type == 'unique'
                    is_fsm = 'fsm' in splicing_events
                    platform = 'pacbio' if read_id.startswith('m') else 'nanopore'

                    transcript_to_reads[transcript_id][sample_idx].append((read_id, is_unique, is_fsm, platform))

    return read_info, transcript_to_reads

def summarize_transcripts(transcript_to_reads, read_info, valid_transcripts, num_samples):
    summary = dict()

    for transcript_id, sample_reads in transcript_to_reads.items():
        all_reads = [r for sample in sample_reads.values() for r in sample]
        total_reads = len(all_reads)
        unique_reads = [r for r in all_reads if r[1]]
        unique_fsm_reads = [r for r in unique_reads if r[2]]

        pct_unique = len(unique_reads) / total_reads if total_reads else 0
        pct_fsm_unique = (len(unique_fsm_reads) / len(unique_reads)) if unique_reads else 0.0

        def count_samples(threshold, fsm_only=False):
            count = 0
            for reads in sample_reads.values():
                filtered = [r for r in reads if r[1] and (r[2] if fsm_only else True)]
                if len(filtered) >= threshold:
                    count += 1
            return count

        has_pacbio = any(r[3] == 'pacbio' for r in all_reads)
        has_nanopore = any(r[3] == 'nanopore' for r in all_reads)

        summary[transcript_id] = (
            total_reads,
            len(unique_reads),
            len(unique_fsm_reads),
            pct_fsm_unique * 100,
            count_samples(1),
            count_samples(3),
            count_samples(10),
            count_samples(1, fsm_only=True),
            count_samples(3, fsm_only=True),
            count_samples(10, fsm_only=True),
            has_pacbio and has_nanopore
        )

        # Perform detailed analysis if <50% unique OR <50% fsm of unique
        if pct_unique < 0.5 or pct_fsm_unique < 0.5:
            assignment_counts = Counter()
            splicing_counts = Counter()
            other_transcript_counts = Counter()

            for sample_read_list in sample_reads.values():
                for read_id, *_ in sample_read_list:
                    read = read_info.get(read_id)
                    if not read:
                        continue

                    assignment_counts[read['assignment_type']] += 1
                    splicing_counts.update(read['splicing_events'])
                    if read['transcript_id'] != transcript_id:
                        other_transcript_counts[read['transcript_id']] += 1

            print(f"\nTranscript: {transcript_id} failed threshold checks")
            print(f"  Total reads: {total_reads}, Unique: {len(unique_reads)} ({pct_unique:.2%}), FSM unique: {len(unique_fsm_reads)} ({pct_fsm_unique:.2%})")
            print(f"  Assignment types: {dict(assignment_counts)}")
            print(f"  Splicing events: {dict(splicing_counts)}")
            print(f"  Other assigned transcripts: {dict(other_transcript_counts)}\n")

    return summary

def main():
    parser = argparse.ArgumentParser(description="Summarize transcript read assignments across experiments.")
    parser.add_argument("--transcripts", required=True, help="File with transcript IDs (one per line).")
    parser.add_argument("--read_files", required=True, nargs='+', help="TSV or TSV.GZ files with read assignments.")

    args = parser.parse_args()

    transcript_ids = load_transcript_ids(args.transcripts)
    read_info, transcript_to_reads = parse_read_files(args.read_files, transcript_ids)

    summary = summarize_transcripts(transcript_to_reads, read_info, transcript_ids, len(args.read_files))

    # Header
    print("Transcript_ID\tTotal_Reads\tUnique_Reads\tFSM_Unique_Reads\tFSM_Unique_%\t"
          "Samples_≥1_Unique\tSamples_≥3_Unique\tSamples_≥10_Unique\t"
          "Samples_≥1_Unique_FSM\tSamples_≥3_Unique_FSM\tSamples_≥10_Unique_FSM\t"
          "Has_Both_PacBio_and_Nanopore")

    # Sort by number of FSM unique reads descending
    for transcript_id, values in sorted(summary.items(), key=lambda x: x[1][2], reverse=True):
        print(f"{transcript_id}\t" + "\t".join(str(v) for v in values))

if __name__ == "__main__":
    main()
