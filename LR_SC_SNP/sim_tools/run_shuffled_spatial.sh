#!/bin/bash
# Negative control for the spatial filter: Visium HD calls with 8 um bin coordinates permuted among bins.
set -euo pipefail
BASE=/abga/work/andreyp/simulation/data/rna_variants
D=Mouse.VisiumHD.full.ONT_R10.brain_expr.50M
S=$(dirname $(readlink -f $0))
LR=/home/andreyp/ngs-analysis/LR_SC_SNP
BIN=$BASE/ref/visium_hd.slide1.bins8um
cd $BASE/$D
OUT=lr_sc_snp
echo "[$(date '+%F %T')] shuffling bin coordinates"
/usr/bin/python3 $S/shuffle_coords.py $BIN.cell_coords.tsv $OUT/shuffled_bins.cell_coords.tsv 1
echo "[$(date '+%F %T')] spatial filter on shuffled coordinates"
PYTHONNOUSERSITE=1 /home/andreyp/miniconda3/envs/splisofind/bin/python $LR/spatial_filter.py --sites $OUT/calls.sites.tsv \
    --barcode_counts $OUT/calls.barcode_counts.tsv.gz --barcode_to_cell $BIN.barcode_to_cell.tsv \
    --cell_coords $OUT/shuffled_bins.cell_coords.tsv -t 16 -o $OUT/calls.shuffled > $OUT/spatial.shuffled.log 2>&1
/usr/bin/python3 $LR/eval_filters.py --eval_sites $OUT/calls.eval.sites.tsv --spatial_sites $OUT/calls.shuffled.spatial_sites.tsv \
    > $OUT/spatial.shuffled.eval.summary.txt
touch $OUT/spatial.shuffled.done
echo "[$(date '+%F %T')] done"
cat $OUT/spatial.shuffled.eval.summary.txt
