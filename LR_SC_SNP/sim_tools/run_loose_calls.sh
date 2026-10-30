#!/bin/bash
# Phase 1 calls at p <= 1e-4 (the default since 2026-10-31) with the unreliable column, plus evaluation,
# written as lr_sc_snp/loose.* next to the older p <= 1e-9 calls (lr_sc_snp/calls.*).
# Usage: bash run_loose_calls.sh DATASET THREADS
set -euo pipefail
D=$1; T=$2
B=/abga/work/andreyp/simulation/data/rna_variants
GENOME=/abga/work/andreyp/data/reference/mouse/GRCm39.primary_assembly.genome.fa
LR=/home/andreyp/ngs-analysis/LR_SC_SNP
cd $B/$D
O=lr_sc_snp/loose
/usr/bin/python3 $LR/call_variants.py --read_info $(ls isoquant/*/*.read_info.tsv.gz) \
    --umi_groups $(ls isoquant/*/*.UMI_filtered.ED*.umi_groups.tsv.gz) --bam $D.bam --reference $GENOME \
    --max_error_pvalue 1e-4 -o $O -t $T 2> $O.log
/usr/bin/python3 $LR/eval_calls.py --truth $D.variants.vcf --truth_cell_counts $D.variant_cell_counts.tsv.gz \
    --sites $O.sites.tsv --barcode_counts $O.barcode_counts.tsv.gz --reference $GENOME -o $O.eval > $O.eval.summary.txt
touch $O.done
