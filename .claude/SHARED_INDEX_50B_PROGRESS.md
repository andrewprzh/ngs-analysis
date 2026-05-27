# 50B Barcode Index Memory Reduction — Progress and Resume

This file captures the state of the multi-step memory-reduction work on
branch `shared_index_10b` so a fresh Claude Code session can resume cleanly.
Authoritative design doc:
`/home/andreyp/.claude/plans/ancient-wishing-tulip.md`.

## Goal

Scale the sparse 3-anchor shared k-mer index
(`isoquant_lib/barcode_calling/indexers/shared_memory.py:SharedMemorySparseAnchorIndexer`)
from the current 10B-barcode / ~0.5 TB RAM ceiling to a 50B-barcode target,
**without further recall loss**. Time and disk are tradeable; recall is not.

## Implementation Plan (User-Approved Order)

| Step | Status | What |
|------|--------|------|
| **D** | ✅ committed `d5c08c6` | New `MappedSparseAnchorIndexer` — disk-backed mmap variant alongside the SHM class |
| **B** | ✅ committed `076a30e` | Stable-sort by first anchor → collapse first-anchor postings to `(start, length)` table (`LAYOUT_VERSION = 2`, new `first_anchor_starts.bin`, inverted list shrunk to 2N×5) |
| **A** | conditional, next | 50-bit packing of `known_bin_seq` — only if profiling after D+B shows it as the binding constraint |
| **C** | conditional | Elias-Fano / PForDelta on remaining 2 posting lists — only if 50B still doesn't fit |
| **E** | reserve | Sharded query batching — only if mmap page-cache thrash proves unacceptable |

User-stated workflow constraint: **commit each step (D, B, ...) separately
with passing unit tests** before moving to the next.

## Where We Are: Steps D and B Committed; Step A Pending Decision

Step D and Step B are both on `shared_index_10b` (commits `d5c08c6` and
`076a30e`). All 24 mapped-index unit tests pass. The integration into
`isoquant_lib/barcode_calling/callers/stereo.py` to actually use the
mapped indexer is still deferred — SHM is still the caller's choice.

**Next decision:** profile a 100M-or-so build with the post-B layout to
see if `known_bin_seq` (uint64, 8 B/barcode) becomes the binding RSS
constraint. If yes → Step A (50-bit packing). If not → consider wiring
the caller to use the mapped indexer at 10B-50B scale (a separate small
commit).

## Historical Notes from Step D (kept for context)

### Files created / modified for Step D

- **New:** `isoquant_lib/barcode_calling/indexers/mapped.py`
  - Classes: `MappedIndexInfo`, `MappedSparseAnchorIndexer`
  - File layout under `cache_dir`:
    - `meta.json` (versioned, `LAYOUT_VERSION = 1`)
    - `known_bin_seq.bin` — N × 8 bytes uint64
    - `index_ranges.bin` — (4^k + 1) × 8 bytes uint64
    - `inverted_list.bin` — `index_size × 5` bytes uint40-packed
  - Reuses the existing Numba kernels from `shared_memory.py`
    (`_compute_and_count_anchors`, `_compute_and_populate_anchors`, etc.)
    by writing into `numpy.memmap`-backed arrays.
  - Best-effort `posix_madvise` (`MADV_RANDOM` / `MADV_SEQUENTIAL`) via
    `_madvise_random` / `_madvise_sequential` helpers.
  - Default `persistent=False` → main instance unlinks `cache_dir` on
    `__del__`. Workers (constructed via `from_sharable_info`) never
    unlink. `persistent=True` skips cleanup for build-once / query-many
    flows.
  - `cache_dir=None` falls back to `tempfile.mkdtemp(prefix="barcode_index_")`
    but callers are expected to pass an explicit fast-disk path for real
    runs.

- **Modified:** `isoquant_lib/barcode_calling/indexers/__init__.py`
  - Added `from .mapped import MappedIndexInfo, MappedSparseAnchorIndexer`
    and extended `__all__`.

- **New tests:** `isoquant_tests/test_mapped_index.py`
  - 18 tests covering: init, meta file contents + sizes, empty index,
    `seq_len` validation, exact match, front/back/overlap error recovery,
    pickling sharable info, layout-version mismatch error,
    `max_hits`/`min_kmers`/`ignore_equal`, score ordering, low-mem
    parity, parity with `SharedMemorySparseAnchorIndexer`, `persistent`
    behavior, `MappedIndexInfo.__getstate__/__setstate__`.
  - All 18 PASS (last verified run logged at
    `local_tmp/mapped_test.log` — 431.02s, the slow part is allocating
    the 4^14-entry ranges table per test; Numba JIT warmup amortizes).

### How to verify mapped-index tests quickly

```bash
mkdir -p local_tmp/pytest
TMPDIR=/home/andreyp/IsoQuant2/local_tmp/pytest \
  /usr/bin/python3 -m pytest isoquant_tests/test_mapped_index.py -v \
  > local_tmp/mapped_test.log 2>&1
# inspect the log: expect "24 passed" (18 from D + 6 from B)
```

**Important:** `TMPDIR` must point under the repo. Without it, pytest's
`tmp_path` lives in `/tmp` (≤ 20 GB free) and gigabyte-scale memmap files
created by the indexer fill the disk → bus error mid-run.

### Commit message to use (per CLAUDE.md style — single-line, lowercase for additions
that "add"; tense: present)

```
Add mmap-backed sparse anchor indexer for 50B-scale barcodes
```

Files to stage:
- `isoquant_lib/barcode_calling/indexers/mapped.py` (new)
- `isoquant_lib/barcode_calling/indexers/__init__.py` (modified)
- `isoquant_tests/test_mapped_index.py` (new)
- `.claude/SHARED_INDEX_50B_PROGRESS.md` (this file, optional —
  decide whether to ship the progress note in-repo or only locally)

The integration into `isoquant_lib/barcode_calling/callers/stereo.py` to
actually use the mapped indexer is **intentionally deferred** — the SHM
class is still wired in. Switching the caller will be a separate small
commit once D and B are both proven.

## Next: Step B

Inside `MappedSparseAnchorIndexer.__init__`, after the anchor-count pass
and *before* writing `known_bin_seq.bin`:

1. Compute the first anchor for every input barcode (cheap — already
   computed in pass 1, but currently discarded). Materialize as a uint32
   array of length N.
2. Build a stable permutation `perm = numpy.argsort(first_anchors, kind="stable")`.
3. Reorder both `seqs` and (implicitly) the per-bucket counts for the
   first anchor: bucket `a` now corresponds to the contiguous slice
   `seqs[start[a] : start[a] + count[a]]`.
4. Write a new file `first_anchor_table.bin` carrying
   `(start: uint64, count: uint40-or-uint64) × 4^k` entries (2.1 GB at
   k=14, constant size).
5. Strip first-anchor entries from `inverted_list.bin`: only last and
   minimizer postings remain → file shrinks by ~1/3.
6. Bump `LAYOUT_VERSION` to 2; meta.json gains a `has_first_anchor_table`
   key.
7. In `get_occurrences`, the first-anchor probe becomes:
   ```python
   start = first_anchor_table[a, 0]
   count = first_anchor_table[a, 1]
   # mark slot_masks[barcode] |= 0b001 for each id in start..start+count
   ```
   No barcode-ID lookup needed — `start..start+count` are the IDs.

### Step B tests to add (alongside the existing ones)

- Permutation invariant: query results equal those of the pre-B mapped
  indexer (built in `low_mem` for both) for the same `BARCODES` list.
- Edge case: many barcodes sharing the same first anchor (force-collide
  with a constructed input).
- Edge case: all barcodes have distinct first anchors.
- Meta file: `layout_version == 2`, `first_anchor_table.bin` exists,
  `inverted_list.bin` is ~2/3 the size of the Step-D version.

### Step B commit message

```
Add first-anchor permutation table to mapped index
```

## Environment Notes for Future Session

- `python` is not on PATH; use `/usr/bin/python3`.
- Bash tool calls writing to `/tmp` are unreliable on this machine and
  the user has explicitly forbidden it: **never** write logs / scratch
  to `/tmp`. Use `./local_tmp/` (already exists). Pytest's `tmp_path`
  fixture is fine because pytest manages it.
- Tests are slow (4^14-bucket allocation per test). 7–8 minutes is
  normal for the full file. Don't over-shorten on retries.
- During this session the Bash tool entered a bad state where every
  command returned spurious exit code 1/2/128, even though the
  underlying subprocesses ran. Restarting the harness should clear it.
  When commands fail mysteriously again, run with output redirected to
  a file under `local_tmp/` and inspect the file with the Read tool.

## Critical Files Map

```
isoquant_lib/barcode_calling/indexers/
├── __init__.py                  # exports — updated for D
├── base.py                      # KmerIndexer, ArrayKmerIndexer (untouched)
├── two_bit.py                   # Dict2BitKmerIndexer, Array2BitKmerIndexer (untouched)
├── shared_memory.py             # existing SHM indexers + Numba kernels (untouched)
└── mapped.py                    # NEW: MappedSparseAnchorIndexer (Step D)

isoquant_tests/
├── test_shared_mem_index.py     # existing SHM tests (untouched, should still pass)
└── test_mapped_index.py         # NEW: 18 tests for D, all passing

isoquant_lib/barcode_calling/callers/
└── stereo.py                    # caller — NOT updated to use mapped class yet
```

## Resume Checklist

1. Read `/home/andreyp/.claude/plans/ancient-wishing-tulip.md` for the
   approved design.
2. Read this file for current state.
3. Verify Step D still passes:
   `/usr/bin/python3 -m pytest isoquant_tests/test_mapped_index.py -v > local_tmp/mapped_test.log 2>&1`
4. Commit Step D as a single-line lowercase / capitalized commit per
   CLAUDE.md style (message above).
5. Implement Step B inside `mapped.py`; add tests; commit.
6. Decide whether Step A is needed based on profiling.
