#!/usr/bin/env python3

import os
import re
import sys
import subprocess
from pathlib import Path

def count_lines_fast(filepath):
    """Fast line counting using wc -l."""
    try:
        result = subprocess.run(['wc', '-l', filepath],
                              capture_output=True, text=True, check=True)
        return int(result.stdout.split()[0])
    except:
        return 0

def count_barcoded_reads(log_dir, log_path):
    """Count lines in barcoded_reads files."""
    # First try to find barcoded_reads files in the output directory
    total = 0
    for f in Path(log_dir).glob("**/*.barcoded_reads*.tsv"):
        total += count_lines_fast(str(f))

    if total > 0:
        return total

    # If not found, try to parse the command line from log to find the barcoded_reads path
    with open(log_path) as f:
        for line in f:
            if "Command line:" in line and "--barcoded_reads" in line:
                # Extract barcoded_reads path from command line
                match = re.search(r'--barcoded_reads\s+(\S+)', line)
                if match:
                    barcoded_path = match.group(1)
                    # Try as absolute path first
                    if os.path.exists(barcoded_path):
                        return count_lines_fast(barcoded_path)
                    # Try relative to log directory
                    rel_path = os.path.join(os.path.dirname(log_dir), barcoded_path)
                    if os.path.exists(rel_path):
                        return count_lines_fast(rel_path)

    return total

def extract_from_log(log_path):
    """Extract statistics from isoquant.log file."""
    stats = {}

    with open(log_path) as f:
        for line in f:
            # Barcode detected
            if "Barcode detected:" in line:
                match = re.search(r"Barcode detected:\s+(\d+)", line)
                if match:
                    stats['barcode_detected'] = int(match.group(1))

            # PCR duplicates filtering stats - get the "and barcoded" versions
            elif "Assigned to any gene and barcoded:" in line:
                match = re.search(r"Assigned to any gene and barcoded:\s+(\d+)", line)
                if match:
                    stats['assigned_to_any_gene_barcoded'] = int(match.group(1))
            elif "Spliced and barcoded:" in line:
                match = re.search(r"Spliced and barcoded:\s+(\d+)", line)
                if match:
                    stats['spliced_barcoded'] = int(match.group(1))
            elif "Uniquely assigned and barcoded:" in line and "spliced" not in line:
                match = re.search(r"Uniquely assigned and barcoded:\s+(\d+)", line)
                if match:
                    stats['uniquely_assigned_barcoded'] = int(match.group(1))
            elif "Uniquely assigned and spliced and barcoded:" in line:
                match = re.search(r"Uniquely assigned and spliced and barcoded:\s+(\d+)", line)
                if match:
                    stats['uniquely_assigned_spliced_barcoded'] = int(match.group(1))
            elif "Unique gene-barcodes pairs:" in line:
                match = re.search(r"Unique gene-barcodes pairs:\s+(\d+)", line)
                if match:
                    stats['unique_gene_barcode_pairs'] = int(match.group(1))
            elif "Total reads saved:" in line:
                match = re.search(r"Total reads saved:\s+(\d+)", line)
                if match:
                    stats['total_reads_saved'] = int(match.group(1))
            elif "Spliced reads saved:" in line:
                match = re.search(r"Spliced reads saved:\s+(\d+)", line)
                if match:
                    stats['spliced_reads_saved'] = int(match.group(1))

    return stats

def extract_from_stats_file(log_dir):
    """Extract statistics from .stats.tsv file if log parsing incomplete."""
    stats = {}

    # Find stats.tsv file
    stats_files = list(Path(log_dir).glob("**/*.stats.tsv"))
    if not stats_files:
        return stats

    # Use the first stats file found
    stats_file = stats_files[0]

    with open(stats_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue
            key, val = parts[0], parts[1]

            if key == "Unique gene-barcodes pairs":
                stats['unique_gene_barcode_pairs'] = int(val)
            elif key == "Total reads saved":
                stats['total_reads_saved'] = int(val)
            elif key == "Spliced reads saved":
                stats['spliced_reads_saved'] = int(val)
            elif key == "Assigned to any gene and barcoded":
                stats['assigned_to_any_gene_barcoded'] = int(val)
            elif key == "Spliced and barcoded":
                stats['spliced_barcoded'] = int(val)
            elif key == "Uniquely assigned and barcoded":
                stats['uniquely_assigned_barcoded'] = int(val)
            elif key == "Uniquely assigned and spliced and barcoded":
                stats['uniquely_assigned_spliced_barcoded'] = int(val)

    return stats

def process_log(log_path):
    """Process a single log file and extract all stats."""
    log_dir = os.path.dirname(log_path)
    sample_name = os.path.basename(log_dir)

    # Count barcoded reads
    total_reads = count_barcoded_reads(log_dir, log_path)

    # Extract from log
    stats = extract_from_log(log_path)

    # Try to fill missing stats from .stats.tsv file
    stats_file_data = extract_from_stats_file(log_dir)
    for key, val in stats_file_data.items():
        if key not in stats:
            stats[key] = val

    # If barcode_detected not found in log but we have total_reads,
    # it means barcoded_reads was provided externally
    barcode_detected = stats.get('barcode_detected', 0)
    if barcode_detected == 0 and total_reads > 0:
        barcode_detected = total_reads

    # Build row
    row = [
        sample_name,
        total_reads,
        barcode_detected,
        stats.get('assigned_to_any_gene_barcoded', 0),
        stats.get('spliced_barcoded', 0),
        stats.get('uniquely_assigned_barcoded', 0),
        stats.get('uniquely_assigned_spliced_barcoded', 0),
        stats.get('unique_gene_barcode_pairs', 0),
        stats.get('total_reads_saved', 0),
        stats.get('spliced_reads_saved', 0),
    ]

    return row

def main():
    log_paths = [
        "/abga/work/andreyp/stereo/isoquant/jan26/IsoQuant_stereo_split.PacBio_S1/isoquant.log",
        "/abga/work/andreyp/stereo/isoquant/jan26/IsoQuant_stereo_split.PacBio_S2/isoquant.log",
        "/abga/work/andreyp/stereo/isoquant/jan26/IsoQuant_stereo_split.S1_4kGenesJunctions/isoquant.log",
        "/abga/work/andreyp/stereo/isoquant/jan26/IsoQuant_stereo_split.S2_PC1-CID_L_Exome_LW/isoquant.log",
        "/abga/work/andreyp/stereo/isoquant/subsampling/Stereo.S1.combined_full/isoquant.log",
        "/abga/work/andreyp/visium/isoquant/jan26/Visium.PacBio.All.2um/isoquant.log",
        "/abga/work/andreyp/visium/isoquant/jan26/Visium.PacBio.All.8um/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.f1.extended_gtf/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.f2.extended_gtf/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.f3.extended_gtf/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.m1.extended_gtf/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.m2.extended_gtf/isoquant.log",
        "/abga/work/andreyp/analysis/sc_anoushka/IsoQuant.SC_Anoushka.m3.extended_gtf/isoquant.log",
    ]

    # Header
    header = [
        "Sample",
        "Total_reads",
        "Barcode_detected",
        "Assigned_to_any_gene_and_barcoded",
        "Spliced_and_barcoded",
        "Uniquely_assigned_and_barcoded",
        "Uniquely_assigned_and_spliced_and_barcoded",
        "Unique_gene_barcode_pairs",
        "Reads_after_dedup",
        "Spliced_after_dedup"
    ]

    # Collect all rows
    rows = [header]

    for log_path in log_paths:
        if not os.path.exists(log_path):
            print(f"Warning: {log_path} not found", file=sys.stderr)
            continue

        try:
            row = process_log(log_path)
            rows.append(row)
        except Exception as e:
            print(f"Error processing {log_path}: {e}", file=sys.stderr)

    # Transpose the table
    transposed = list(zip(*rows))

    # Print transposed table
    for row in transposed:
        print("\t".join(map(str, row)))

if __name__ == "__main__":
    main()
