# LR_SC_SNP: somatic variants in long-read single-cell and spatial RNA

A proof-of-concept caller built on top of IsoQuant. It works in two phases:

1. **Raw calls** (`call_variants.py`): simple UMI-level counts with hard cut-offs, tuned for recall.
2. **Pattern filters**: does a variant follow a cell type (`celltype_filter.py`) or a spatial pattern (`spatial_filter.py`)?
   These remove most false raw calls.

The idea: bulk callers (Clair3-RNA, ClairS-TO) miss everything below about 5–10% pseudobulk VAF. In single-cell or spatial
data a rare variant becomes believable when separate cells carry it, and more so when the carriers form a pattern.
Calls are made on the whole-sample pseudobulk; per-group pseudobulks are not used, since spatial patterns need not
follow tissue regions.

Status (2026-10-31): benchmarked on simulated ONT data (10x, 10x with cell types, Visium HD); not yet run on real data.
Nothing is committed yet (IsoQuant branch `umi_groups_output`, these scripts).

## 0. IsoQuant
Run IsoQuant in single-cell or spatial mode with `--large_output read_info umi_groups`.
`umi_groups` lists the reads of every UMI deduplication group, and is available on the IsoQuant branch `umi_groups_output`.

Example for simulated ONT 10x data that is already aligned:
```
isoquant.py --reference genome.fa --genedb genes.db --complete_genedb --bam reads.bam --data_type nanopore \
    --mode tenX_v3 --split_molecules false --barcode_whitelist whitelist.tsv \
    --large_output read_info umi_groups -t 16 -o OUT --prefix SAMPLE
```
- **Outputs used:** `SAMPLE.read_info.tsv.gz`, `SAMPLE.UMI_filtered.ED3.umi_groups.tsv.gz`, and the input BAM.
- **Aligning in IsoQuant:** if IsoQuant aligns FASTQ itself, use its BAM in `OUT/SAMPLE/aux/`, not the `_dedup_` or `_tagged_` ones.
- **Spatial modes:** IsoQuant writes a dense RNA velocity loom for every barcode grouping. For Visium HD (11.2M barcodes)
  that is about 3.5 TB. It is being fixed in IsoQuant separately; use an IsoQuant version with the fix for spatial data.

## 1. Phase 1: `call_variants.py`
```
call_variants.py --read_info SAMPLE.read_info.tsv.gz --umi_groups SAMPLE.UMI_filtered.ED3.umi_groups.tsv.gz \
    --bam reads.bam --reference genome.fa -o PREFIX -t 16
```
- **Reads:**
  - Only reads assigned to a single gene are used (`--gene_assignment_types`).
  - The BAM must be the original one with all reads, not the deduplicated BAM.
  - Only primary alignments are used. Exons reached only by supplementary alignments are not seen; see the open items.
- **Variants:**
  - SNVs, and indels up to `--max_indel_len` (10) bp.
  - Indels are left-aligned. Indels next to introns, clips or alignment ends are ignored.
- **One UMI is one piece of evidence.** At each site, consider the reads of a molecule that cover it:
  - the molecule is **alt** if strictly more than `--consensus` (0.5) of them carry the allele;
  - it is **ref** if strictly more than 0.5 carry the reference;
  - otherwise it is a **no call**, counted neither in alt nor in depth.
- **Reporting:** a variant is reported when all of these hold:
  - at least `--min_alt_umis` (10) alt molecules, from at least `--min_alt_barcodes` (10) barcodes;
  - VAF at least `--min_vaf` (0.005);
  - more alt molecules than sequencing errors explain: the binomial p-value of alt out of depth molecules, at the
    variant class's error rate, is at most `--max_error_pvalue` (1e-4).
- **Error classes:** SNVs by transcript-strand change; indels by type, length (1, 2, 3+) and homopolymer length (1–5, 6+).
- **Error rates** are estimated before calling, on random loci holding `--error_sample_fraction` (5%) of molecules:
  - they are per-read errors at every covered position, from single-read molecules;
  - alleles in at least 10% (and at least 3) of the reads at a position are left out as likely variants;
  - a class with fewer than 100 errors takes the rate of its closest pooled parent (e.g. `DEL:3+:hp6+` → `DEL:3+` → `DEL`);
  - `--error_rates` reuses a table from an earlier run instead;
  - multi-read molecules need a read majority, so their error rate is lower and the test is conservative.
- **Why the error test:**
  - Without it, ONT errors (about 0.3–0.4% per specific substitution) give more than 10 alt UMIs from more than 10
    cells at any site deeper than about 3,000 UMIs. On the full 10x simulation that was 4.4M calls for about 6K true SNVs.
  - About 1M sites are tested per 20M-read dataset, so at p ≤ 1e-4 the error model alone predicts about 100 false calls.
  - The default was 1e-9 until 2026-10-31. It was relaxed because 1e-4 adds about 20% true calls on the cell-type
    simulation, and the cell-type filter removes nearly all the added false ones (see the benchmark).
- **Unreliable flags.** Nothing is removed; the column `unreliable` lists reasons to doubt a site, or `.`:
  - **`read_end`:** in at least `--end_fraction` (0.5) of the alt reads, the site is within `--end_distance` (20)
    aligned bases of an alignment end. Introns are skipped and soft clips don't count.
    - This catches adapter bases aligned past transcript ends: the TSO (ending in GGG) at the 5' end and the polyA at
      the 3' end. Nearly every full-length read then shows the "variant", so VAF approaches 100%.
  - **`editing`:** A>G on the transcript strand in at least half of the stranded alt molecules, i.e. genomic A>G in
    + strand molecules or T>C in − strand ones.
    - The strand is that of each molecule's assigned isoform (read_info), so overlapping genes on opposite strands
      are handled per molecule. The same rule sets `editing_like`.
  - **`homopolymer`:** an indel overlapping or next to (1 bp) a reference run of at least `--homopolymer_flag_len` (5)
    identical bases.
  - **Not included:** a "cluster" flag (other calls within 50 bp). It marked 60–90% of false calls in the simulation,
    but real somatic SNVs and editing sites can cluster.
- **Outputs:**
  - **`PREFIX.sites.tsv`:** one line per allele. Columns:
    - gene, alt genes, strand, transcript-strand change;
    - `editing_like`, `homopolymer_len` (the reference run at the site);
    - alt / depth / no-call UMIs, VAF, alt barcodes, covered barcodes;
    - error class, error rate and the error test's p-value;
    - `alt_end_fraction` (share of alt reads with the site within `--end_distance` of an alignment end);
    - `unreliable`.
  - **`PREFIX.barcode_counts.tsv.gz`:** variant, barcode, alt UMIs and depth UMIs, for every barcode with depth > 0.
  - **`PREFIX.error_rates.tsv`:** per class: errors, depth, excluded sites and rate, with pooled rates per type and length.
- **Algorithm.** Overlapping genes are merged into loci, and each locus is read up to three times:
  1. Count molecules with at least one alt read, an upper bound that selects the candidates.
  2. Compute the molecule consensus at the candidates: totals and alt molecules per barcode. Apply the cut-offs and the error test.
  3. Only for loci with reported sites: count depth per barcode, alt read positions and alt molecule strands at those sites.
  - Depth per barcode is left out of pass 2 on purpose. At deep loci, pass 1 keeps nearly every covered position
    (errors alone pass it), and almost every (candidate, barcode) pair is covered. Counting them all took most of the
    time and up to 50 GB per worker.
  - A full 20M-read 10x dataset takes 35–50 min on 30 threads, and the 50M-read Visium HD dataset about 2 h.

## 2. Phase 2 filters
- **`celltype_filter.py`:** `--sites --barcode_counts --cell_types barcode_types.tsv -o PREFIX`.
  - For each variant and type: a one-sided Fisher's exact test of alt vs ref UMIs in the type against all other typed
    barcodes. BH correction is across all variant × type tests.
  - Outputs: `PREFIX.celltype_tests.tsv` and `PREFIX.celltype_sites.tsv` (best type and the significant types).
- **`spatial_filter.py`:** `--sites --barcode_counts --barcode_to_cell map.tsv (--cell_coords | --barcode_coords) -o PREFIX`.
  - Barcodes are summed per cell, and every covering cell gets alt / depth.
  - The test is Spl-IsoFind's `moransI_sparse`: kNN Moran's I, a permutation p-value and a BY q-value.
  - The Spl-IsoFind parameters `--k --nperm --mincells --mincellspergroup --imb` keep the package defaults unless set.
    The default `mincellspergroup=20` means a variant needs at least 20 cells with alt fraction above 0.5 to be tested.
  - Output: `PREFIX.spatial_sites.tsv`.
- **`make_bin_cells.py`:** groups grid barcodes into square bins as pseudo-cells, e.g. Visium HD 2 µm into 8 µm when
  there is no segmentation.

## 3. Running on real data (to be completed when the data is provided)
The first target is PacBio Kinnex long reads with Visium HD, from one mouse brain section (GRCm39, GENCODE vM36 basic).
The steps and inputs, as far as they are known now:
1. **IsoQuant:** branch `umi_groups_output` plus the loom fix, `--data_type pacbio_ccs --mode visium_hd`, the
   Visium HD barcode whitelists and `--large_output read_info umi_groups`.
   - Decide `--split_molecules`. Kinnex reads are already segmented by skera, and the simulations used `false`.
2. **Phase 1:** `call_variants.py` with defaults. The error rates are re-estimated from the data, since PacBio errors
   differ from ONT.
3. **Spatial filter** needs IsoQuant barcode → cell and cell → coordinates. For the existing real section:
   - IsoQuant reports Visium HD barcodes as `BC1|BC2` sequences. The spaceranger table
     `/abga/work/andreyp/data/visium/spaceranger/barcode_seq_to_2um.tsv` maps them to 2 µm bins (`s_002um_RRRRR_CCCCC`).
   - The segmentation `FromJustine/InDel_matrix_gen/Bin_to_Cell.csv`, in `/abga/work/andreyp/data/spatial_indel/`, maps
     2 µm bins (column `barcode`) to `cell_id`, with `in_nucleus` and `in_cell`. It has 77,759 cells.
   - `Cell_Labels.csv`, in the same folder, gives `cell_id → x, y, celltype`. About 47K cells are typed "other".
   - The two tables chain into `--barcode_to_cell` (IsoQuant barcode → cell_id) and `--cell_coords` (cell_id, x, y).
4. **Cell-type filter:** barcode → cell type, through the same chain.
5. **What to expect, and how to read the flags:**
   - RNA editing: A>G on the transcript strand, often clustered; flagged `editing`.
   - Adapter artifacts at transcript ends: `read_end`.
   - Homopolymer indels, which are PacBio's main error mode: `homopolymer`.
   - Germline differences from the reference: high VAF in all cells, so they should not pass the pattern filters.
   - Paralog mismapping: clusters of calls with identical alt counts. There is no flag for it yet.

## 4. Benchmark on simulated data
- **Data:**
  - Datasets are in `/abga/work/andreyp/simulation/data/rna_variants`, made with
    `~/ngs-analysis/barcode_detection/simulate_barcoded_mutations.py` and NanoSim ONT R10 reads; see the README there.
  - The driver is `run_lr_sc_snp.sh DATASET [THREADS]` in that folder: IsoQuant, phase 1, evaluation, filters.
  - Outputs are in `DATASET/lr_sc_snp/`:
    - `calls.*`: the run at p ≤ 1e-9, before the flags existed;
    - `loose.*`: the 10x runs at p ≤ 1e-4 with flags, the current defaults (`sim_tools/run_loose_calls.sh`);
    - `loose.cutoffs.txt`: results per cut-off.
- **Results at the current defaults** (Visium HD was run at p ≤ 1e-9):

| dataset | SNV precision | SNV recall | indel precision / recall | phase 2 filter |
|---|---|---|---|---|
| 10x 5K, 20M reads, random carriers | 0.41 | 0.60 | no indels simulated | random cell-type labels: 0 pass |
| 10x HippSNf1, 20M reads, cell-type carriers | 0.24 | 0.28 | INS 0.35 / 0.26, DEL 0.09 / 0.27 | cell type: SNV precision 0.98, recall 0.26 |
| Visium HD, 50M reads, stripe carriers (p ≤ 1e-9) | 0.37 | 0.69 | no indels simulated | spatial: precision 0.996, 7,502 of 9,366 true calls pass; shuffled bins: 0 of 13,671 pass |

- **Recall by support:**
  - At p ≤ 1e-4, recall is at least 98.7% at ≥50 alt templates, 90–96% at 20–49 and 26–31% at 10–19 (10x and
    cell-type). Visium HD at p ≤ 1e-9: 98–99.7% at ≥50, 82% at 20–49.
  - Variants with fewer than 10 alt templates (62% of the cell-type truth, 32% of the 10x truth) are below the
    10-UMI cut-off by design.
  - About 37% of templates are never sequenced at the site (the simulator samples about one read per template), so a
    variant with 10–19 alt templates has only 6–11 alt molecules.
- **Error-test cut-off** (cell-type simulation; SNV recall in phase 1 and after the cell-type filter):

| p cut-off | SNV recall | after filter: recall / false SNVs | false DEL in phase 1 |
|---|---|---|---|
| 1e-9 | 0.225 | 0.219 / 70 | 1,903 |
| 1e-6 | 0.260 | 0.248 / 74 | 2,498 |
| 1e-4 | 0.279 | 0.263 / 75 | 3,631 |

- **False raw calls**, traced to their simulated templates through the read names (`sim_tools/screen_*.py`):
  - **Transcript ends:** about 45% of false SNVs. Adapter bases are aligned 1–15 bp past the transcript end (`read_end`).
  - **Paralogs:** about 43% of false SNVs and most false insertions. Reads from a paralog, e.g. in the Pcdha cluster,
    align inside an intron with about 12% mismatches (against about 3%) and give clusters of calls.
  - **Sequencing errors at deep sites:** mostly deletions at VAF around 1.5%, many in homopolymers (`homopolymer`).
  - **Templates match the genome** except for the adapters, so there is no reference-bias artifact.
  - **Flags at p ≤ 1e-4:** they mark 38% of false SNVs, 34% of false deletions and 19% of false insertions. Among true
    calls, 9% of SNVs are flagged (nearly all `editing`: random A>G is 1/12 of simulated SNVs), 5% of deletions, 3% of insertions.
- **Missed variants with ≥5 alt templates:**
  - Most fail the 10-UMI cut-off. The rest failed the old 1e-9 error test (VAF about 2% at depth 600–800) or the barcode cut-off.
  - About 2% of true variants have almost no reads: exons reached only by supplementary alignments in very long genes
    (Auts2, Cadm2, Nkain2), a paralog cluster at chr2:175 Mb, and the first or last 10 bp of transcripts.

## 5. Simulation tools (`sim_tools/`)
These are debugging and benchmark scripts for the simulated datasets, with hard-coded paths.
- **`screen_sites.py DATASET_DIR OUT_PREFIX [THREADS] [SITES_PER_GROUP]`** (1–2 h): traces the reads at sampled
  TP / FP / FN sites to their templates. It writes per-site and per-read features, the counts the caller would see,
  and the template funnel for true variants.
- **Summaries of the screen:**
  - `screen_summary.py PREFIX`: the recall funnel, the first failing cut-off and FP vs TP features;
  - `screen_mechanism.py PREFIX`: where false alt reads lie relative to their source transcript;
  - `screen_flags.py PREFIX`: shares of FP and TP calls marked by candidate flags.
- **`run_loose_calls.sh DATASET THREADS`:** reruns phase 1 and evaluation into `lr_sc_snp/loose.*`.
- **`eval_cutoffs.py DATASET_DIR CELL_TYPES`:** precision and recall per p cut-off, with and without flags and after
  the cell-type filter.
- **`hp_flag_check.py EVAL_SITES`:** homopolymer flag variants on labelled indel calls.
- **Negative control for the spatial filter:** `shuffle_coords.py` and `run_shuffled_spatial.sh` permute the 8 µm bin
  coordinates.

## 6. Open items
- **Duplicate rows:** a few sites are reported twice in `sites.tsv`, once from each of two adjacent loci whose reads
  cover the same position: 6–38 rows per dataset, mostly false calls.
  - The phase 2 filters sum the barcode rows per variant, which is correct because the loci's molecules are disjoint.
    Only the sites table has the duplicate rows.
- **Supplementary alignments** are ignored, both here and in IsoQuant, which loses exons of very long genes.
- **Paralog mismapping** has no flag. Candidate signals from the simulation: alt reads noisier than ref reads, and
  low MAPQ.
- **Pass 1** could use a stricter bound. At deep loci nearly every covered position becomes a candidate.

## Environments
- **`call_variants.py`, `celltype_filter.py`, `eval_calls.py`:** pysam, numpy and scipy (`/usr/bin/python3` with the user site).
- **`spatial_filter.py`:** the conda env `splisofind` (`pip install SplIsoFind`; 0.1.8). Run it with `PYTHONNOUSERSITE=1`.
- **`sim_tools/screen_*.py`:** pandas in addition.
