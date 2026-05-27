# PolyA / TSS site prediction — implementer reference

Implementation: `isoquant_lib/terminal_counter.py`. Ported on top of
master's `CompositeCounter` / `string_pools` rewrite from the
`start_position` branch (preserved as the `start_position_pre_merge`
tag).

> Companion doc: `.claude/POLYA_TSS_TRAINING.md` covers how to collect
> training data and retrain the XGBoost models.

## What it does

For each reference transcript, accumulate a per-base histogram of read
terminal genomic coordinates (read end for polyA, read start for TSS),
detect peaks with `scipy.signal.find_peaks`, drop the ones an XGBoost
classifier rejects, and emit one row per accepted peak with a
`Known` / `Novel` flag computed against the annotated transcript end.

Gating:

| Flag           | PolyA | TSS                  |
| -------------- | ----- | -------------------- |
| `--genedb`     | on    | on                   |
| `--fl_data`    | —     | also required        |

TSS additionally requires `--fl_data` because read start coordinates
without full-length evidence are unreliable (5′ degradation, RT
fall-off).

## Pipeline walkthrough

`add_read_info(read_assignment)` is called once per assignment from the
per-chromosome worker (`parallel_workers.collect_reads_in_parallel`,
line ~276). It applies the filter, extracts the position, and pushes
into `self.transcripts[transcript_id]['data']` (plus the per-group
bucket).

`dump()` is called once per chromosome (line ~300) after every
assignment for that chromosome has been processed. The whole peak
detection / classification / write pipeline runs there:

1. Build a `pandas.DataFrame` row per transcript with `data_min`,
   `data_max`, the bin-1 histogram of all positions, the histogram
   summary statistics (`var`, `skew`, `entropy`, `mean_height`).
   Histogram is padded with `HISTOGRAM_PAD` zeros on each side so
   `find_peaks` has room to detect edge peaks and `peak_widths` does
   not run off the array.
2. Run `find_peaks(histogram, distance=PEAK_DISTANCE)`, then
   `peak_prominences` and `peak_widths(..., rel_height=PEAK_REL_HEIGHT)`
   on the same peak indices.
3. Transcripts with zero detected peaks: fall back to the
   `scipy.stats.mode` of the position list. One row emitted, peak left
   = peak right = mode.
4. Transcripts with ≥1 peak:
   - `_rank_peaks` sorts peaks by descending height, computes
     `relative_height = h / h_top`, and assigns `rank` starting at 1.
     One-peak rows get `rank=0, relative_height=1.0`.
   - `df.explode(...)` flattens the per-transcript arrays into one row
     per (transcript, peak).
   - Feed `FEATURE_COLUMNS` into the loaded XGBoost model, drop rows
     where `predict()` returns 0.
   - `prediction` column is set to the peak_location (start-relative);
     it becomes a genomic coordinate after adding `start`.
5. Survivors from steps 3 and 4 are concatenated, classified with
   `flag = 'Novel' if abs(prediction - annotated) > ANNOTATION_TOLERANCE`,
   and the per-peak count window is computed by `_counts_for_peak`.
6. Write to the per-chr TSV (always `mode='w'` — each worker has its
   own per-chr file). The sample-level counter never calls `dump()`;
   it only provides `output_counts_file_name` for the
   `merge_counts → merge_files` concatenation in
   `dataset_processor.merge_assignments`.

## Tunable constants (must match training)

Defined at module top in `isoquant_lib/terminal_counter.py`:

| Name                   | Value | Used for                                      |
| ---------------------- | ----- | --------------------------------------------- |
| `PEAK_DISTANCE`        | 10    | minimum spacing between peaks (bp)            |
| `PEAK_REL_HEIGHT`      | 0.98  | rel-height for `peak_widths`                  |
| `HISTOGRAM_PAD`        | 10    | zero-padding on each side of the histogram    |
| `ANNOTATION_TOLERANCE` | 10    | bp window for `Known` / `Novel` classification |

`FEATURE_COLUMNS = ['var', 'skew', 'peak_count', 'peak_width', 'entropy',
'mean_height', 'peak_heights', 'relative_height']`. If you change any
of these or the feature names, retrain the model in lockstep — the
classifier has no schema check at inference time and will silently
mis-score peaks if columns drift.

## Filter logic

`_ACCEPTED_ASSIGNMENT_TYPES` =
`{unique, unique_minor_difference, inconsistent, inconsistent_non_intronic}`.

Counter-specific:
- **PolyA**: also requires `read_assignment.polyA_found` and a
  non-null `polya_info`. Position is `external_polya_pos` for `+`
  strand, `external_polyt_pos` for `-` strand; the read is dropped if
  the matching field is -1.
- **TSS**: requires strand `+` or `-` (rejects `.`); position is
  `corrected_exons[0][0]` for `+`, `corrected_exons[-1][1]` for `-`.

Annotated position is read from
`read_assignment.gene_info.all_isoforms_exons[transcript_id]`:

| Counter | `+` strand                 | `-` strand                 |
| ------- | -------------------------- | -------------------------- |
| PolyA   | `exons[-1][1] + 1`         | `exons[0][0] - 1`          |
| TSS     | `exons[0][0] - 1`          | `exons[-1][1] + 1`         |

## Multi-group integration

Each counter holds a `group_index` and reads from
`read_assignment.read_group_ids[self.group_index]` (a list of
**integer** pool indices since master's string-interning rewrite —
not the original branch's bare strings).

At dump time, integer IDs are translated back to display names via
`string_pools.get_read_group_pool(self.group_index).get_str(id)`.

`ReadAssignmentAggregator.__init__` (`isoquant_lib/assignment_aggregator.py`)
constructs counters in two passes:

1. Always: one ungrouped `PolyACounter` (and, with `--fl_data`, one
   ungrouped `TSSCounter`), constructed with `string_pools=None`.
   Writes to `sample.out_polya_prediction_tsv` /
   `sample.out_tss_prediction_tsv`.
2. If `--read_group`: one grouped counter per grouping strategy with
   `string_pools=self.string_pools, group_index=group_idx`. Writes to
   `f"{sample.out_polya_prediction_grouped_tsv}_{strategy_name}"`.

Both passes append to `self.global_counter`; everything else (per-chr
parallel execution, file merge, resume) is handled by the existing
counter infrastructure unchanged.

## Output format

Ungrouped TSV columns:
```
chromosome  transcript_id  gene_id  prediction  counts  flag
```

Grouped TSV (one row per peak × group with non-zero count):
```
chromosome  transcript_id  gene_id  prediction  counts  flag  counts_byGroup  group_id
```

- `prediction`: genomic coordinate of the predicted site.
- `counts`: total reads supporting the peak window (`peak_left ≤ pos ≤
  peak_right` in the padded histogram).
- `counts_byGroup` (grouped only): reads from this group in the
  genomic window `[start + peak_left, start + peak_right]`.
- `flag`: `Known` if `|prediction − annotated| ≤ ANNOTATION_TOLERANCE`,
  else `Novel`.

Per-chr files (`SAMPLE_chrid.polyA_prediction.tsv` etc.) live alongside
the merged ones during execution and are deleted by `merge_files` after
concatenation.

## Counter / `AbstractCounter` integration contract

`TerminalCounter` does **not** call `AbstractCounter.__init__`. The
base class's `counts_file_name` suffix machinery would clobber our
output path (which is already a fully-formed `.tsv`). Instead we set
the four fields the rest of the pipeline reads — `output_file`,
`output_counts_file_name`, `output_stats_file_name = None`,
`usable_file_name = None` — and clear the file in our own `__init__`.
With both stats/usable fields `None`, `merge_counts()` skips the
stats and usable handling and just concatenates the per-chr TSVs.

The remaining `AbstractCounter` methods (`add_read_info_raw`,
`add_confirmed_features`, `add_unassigned`, `add_unaligned`,
`finalize`) are explicit no-ops — they exist because the
`CompositeCounter` iterates over every counter for each of these
calls, but they have no meaning for terminal-position prediction.

## Lazy model loading (fork-safety)

`XGBClassifier` initialises the OpenMP runtime when a model is loaded.
If the parent process loads the model before forking workers, the
inherited OpenMP semaphores can deadlock inside the workers and the
pipeline hangs silently in `process_assigned_reads`. The counter
defers model loading to a `model` property that triggers on first
`predict()` call — i.e. only inside the per-chr worker, after fork.
The parent's sample-level counter never reaches `dump()`, so its
`_model` stays `None`.

Symptom if this regresses: `-t 1` works fine, `-t 2` hangs partway
through "Processing assigned reads" with empty per-chr output files.

## Model files

- `isoquant_lib/data/model_polya.json` (~261 KB)
- `isoquant_lib/data/model_tss.json` (~310 KB)

Shipped with the package via `[tool.setuptools.package-data]` in
`pyproject.toml` and `recursive-include isoquant_lib/data *.json` in
`MANIFEST.in`. Paths are resolved with
`Path(__file__).parent / "data" / ...` — independent of cwd, survives
`pip install` and the `src/ → isoquant_lib/` rename.

## Edge cases handled

- `_counts_for_peak` clamps `[peak_left, peak_right]` to histogram
  bounds (covers edge peaks where `peak_left < 0` or
  `peak_right ≥ len(hist) − HISTOGRAM_PAD`).
- `_counts_for_group` returns 0 if the group never contributed to the
  transcript.
- `_write_empty` always writes the header (column-only TSV) so
  downstream pipelines can rely on the file existing with a header
  row even when no transcripts had eligible reads.
- Empty input (`self.transcripts == {}`) short-circuits to
  `_write_empty` without touching the model — useful for chromosomes
  with no polyA-tailed reads.
- Concat of zero-peak fallback + model-survivor frames filters out
  empty frames first, avoiding the pandas `FutureWarning` about
  all-NA columns in `pd.concat`.

## Performance notes

- Memory: `self.transcripts` holds every contributing read's terminal
  position as a Python int in lists keyed by transcript and group.
  Worst case = O(reads × (1 + n_groups)). For a 1 B-read run with
  ~200 K transcripts and 100 groups this is the dominant memory cost
  for the counter — but it's per-chr and per-worker, so it scales
  down by chromosome count and clears on `dump()`.
- CPU: `dump()` is dominated by the per-row `.apply(lambda …)` calls
  building the peak DataFrame. Acceptable today but the obvious target
  if profiling shows the prediction step becoming the long pole.
  Vectorizing `find_peaks`/`peak_widths`/`peak_prominences` across the
  whole histogram array is non-trivial (each transcript has a
  different histogram length); a more realistic optimization is to
  batch the XGBoost `predict()` calls across transcripts (already
  done after `explode`).
- I/O: per-chr TSVs use `mode='w'`, so each worker writes its own
  file exactly once; the merge step is straight `shutil.copyfileobj`.

## Training-data collection (developer mode)

Hidden flags `--collect_polya_training PATH` and
`--collect_tss_training PATH` switch the matching counter from
inference to training-feature emission. In that mode `dump()` runs
the same peak-detection + feature-extraction pipeline but skips the
XGBoost `predict()` call and appends per-peak rows
(`FEATURE_COLUMNS + ['chromosome', 'true_peak']`) to a per-chr CSV
next to the per-chr prediction TSV. `dataset_processor.merge_assignments`
concatenates fragments into the user-supplied PATH and removes the
per-chr files. Grouped counters are skipped entirely in training
mode. Predictions emit a header-only TSV.

`isoquant.py:main` prints a multi-line WARNING when either flag is
present so the dev mode is obvious in the log.

End-to-end workflow:

```
python isoquant.py … --genedb GENEDB --collect_polya_training peaks.csv
python misc/train_polya_tss_model.py --features peaks.csv \
    --output isoquant_lib/data/model_polya.json
```

Full procedure, schema, hyperparameters, validation checklist, and
model-replacement workflow live in `.claude/POLYA_TSS_TRAINING.md`.

## Tests

`isoquant_tests/test_polya_prediction.py` — 22 unit tests covering
filtering, accumulation per group, annotation lookups for both
strands and both counters, dump structure (empty / ungrouped /
grouped / model-rejected / multi-transcript / zero-peak fallback),
peak ranking, count-window clamping, Known/Novel boundary,
no-op-method safety, output-file truncation on init. XGBoost is
mocked via a `_StubModel` so the tests do not depend on the shipped
classifier.

`isoquant_tests/console_test.py::test_with_bam_and_polya` asserts
the ungrouped TSV is produced and non-empty.
`test_with_bam_polya_and_fl_data` adds the `--fl_data` path and
asserts both polyA and TSS TSVs.

End-to-end accuracy is gated by the weekly
`.github/workflows/Mouse.ONT_simulated.polyA.yml` workflow. It runs
IsoQuant against the GENCODE vM36 annotation on simulated mouse reads
generated from a synthetic GTF with alternative transcript ends,
then `misc/assess_polya_prediction.py` compares the resulting
`*.polyA_prediction.tsv` to the synthetic GTF's annotated polyA
coordinates. CI YAML lives at
`/abga/work/andreyp/ci_isoquant/data/Mouse.ONT_simulated.polyA.yaml`
(on the self-hosted runner). New launcher run-type:
`polya_prediction` (`isoquant_tests/github/run_pipeline.py`).

### Assessment script (`misc/assess_polya_prediction.py`)

**Site-level evaluation.** Each truth polyA site (known or novel) is an
independent target. A prediction is TP if it matches any truth site
within `--tolerance` (default 10 bp, mirrors `ANNOTATION_TOLERANCE`).
A truth site is recovered if any prediction matches it. Precision is
prediction-centric, recall is site-centric.

**Novel ID handling.** The synthetic GTF contains novel transcript
entries with IDs like `ENSMUST00000000033.12_142204503` (base GENCODE
ID + `_` + novel polyA coordinate). `load_truth` merges these into
the base transcript's truth list so a single prediction for the base
ID can match either the known or the novel site.

**Expression filter (`--transcript_counts`).** When provided (e.g. a
simulated read count TSV), both the truth set and predictions are
filtered to only include expressed transcripts (count > 0). This
prevents unexpressed transcripts in the reference GTF from inflating
the recall denominator. Novel IDs are mapped to base IDs before
filtering. The YAML config key is `simulated_counts`; the launcher
(`run_polya_prediction`) passes it as `--transcript_counts`.

**Known / novel split.** The report includes overall precision /
recall / F1 plus per-category breakdowns:
- `precision_known`, `recall_known`, `f1_known` — sites from standard
  GENCODE transcript entries.
- `precision_novel`, `recall_novel`, `f1_novel` — sites from synthetic
  `base_<digits>` entries.

For precision, each prediction's closest truth site determines its
known/novel category. For recall, each truth site is checked
independently against all predictions for that transcript.

## Port notes (one-time, for archaeology)

The feature was implemented on the `start_position` branch (forked
from 3.8.0, 4acb5d5) in `src/terminal_counter.py` and
`src/polya_position_model.py`. Master had moved on substantially
(rename `src/` → `isoquant_lib/`, multi-group `CompositeCounter`,
string-interning rewrite, barcode-calling pipeline, fusion detection,
pip packaging). Rather than `git merge master`, a fresh
`polya_tss_on_master` branch was cut off master and the polyA/TSS
feature was ported on top with the following cleanups:

- `polya_position_model.py` was a near-duplicate of `terminal_counter.py`
  with drifted class signatures; tests imported the duplicate while
  production imported the real module. Dropped entirely.
- `collect_data()` and `result()` were 95 %-identical to `create_df()`
  and never called. Dropped.
- `print(matrix)` and other debug prints replaced with `logger.debug`.
- Hardcoded `src/model_polya.json` / `src/model_tss.json` paths replaced
  with `Path(__file__).parent / "data" / ...`.
- Inline training-data writes (`self.peaks.to_csv("src/model_df.csv")`)
  removed from the hot inference path. The original branch's
  toggle-style mechanism (comment out `create_df()` / uncomment
  `collect_data()` in `dump()`, uncomment `train_model()` in
  `finalize()`) was preserved in spirit but moved behind the hidden
  `--collect_polya_training` / `--collect_tss_training` flags, with
  the trainer extracted to `misc/train_polya_tss_model.py`. The
  46 MB `model_df.csv` artifact was not ported; gitignored.
- Unused imports (`pickle`, `matplotlib.pyplot`, `Enum`, `unique`,
  `random`) removed.
- Magic numbers extracted to module-level constants.
- `read_assignment.read_group` (single string) replaced with
  `read_assignment.read_group_ids[self.group_index]` (integer pool
  index) — matches master's CompositeCounter contract.
- Type hints added; module-level / class / method docstrings added.
- Model loading made lazy to keep XGBoost's OpenMP out of the parent
  process before fork.
- Counter constructor takes `string_pools` and `group_index` like the
  other master-side counters; `AbstractCounter.__init__` is bypassed
  (its `counts_file_name` suffix machinery would mangle the
  `.polyA_prediction.tsv` path).
- Wired into `ReadAssignmentAggregator.__init__` rather than directly
  into `DatasetProcessor` — that's where master constructs counters.
- New tests rewritten against the production module (the old tests
  exercised the dead duplicate).
- `xgboost` and `scikit-learn` moved into `requirements.txt` and
  `pyproject.toml` runtime deps (were test-only on the original
  branch, which meant any production install hit `ImportError` as
  soon as `--genedb` was supplied).
