# SharedMemoryArray2BitKmerIndexer Optimization Plan

## Problem Summary

| Version | Time (10M barcodes) | Peak RAM | Issue |
|---------|---------------------|----------|-------|
| Original (list of lists) | ~4 min | High | Memory overhead from Python lists |
| Two-pass (80bd928) | ~8 min | Lower | 2x slower due to sorting/loop overhead |
| Vectorized (018222d) | ~8 min | Higher | Optimized wrong part, added memory |

## Root Cause Analysis

1. **Pass 1 (counting) was never the bottleneck** - the vectorization benchmark showed speedup here, but this was already fast
2. **Pass 2 bottleneck is the Python loop** (lines 153-160) that iterates over unique k-mers per chunk
3. **Sorting overhead:** `numpy.argsort` is O(n log n), called per chunk
4. **Memory overhead:** Multiple temporary arrays created per chunk

## Proposed Solutions (in order of preference)

### Option 1: Numba JIT for Pass 2 (Recommended)

Use Numba to JIT-compile the index population loop:

```python
from numba import njit

@njit
def populate_index_numba(seqs, shifts, mask, index_ranges, index, num_kmers_per_seq):
    """Single-pass direct write with position tracking."""
    total_kmers = len(index_ranges) - 1
    write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

    for seq_idx in range(len(seqs)):
        seq = seqs[seq_idx]
        for k_pos in range(num_kmers_per_seq):
            kmer_idx = (seq >> shifts[k_pos]) & mask
            write_pos = index_ranges[kmer_idx] + write_positions[kmer_idx]
            index[write_pos] = seq_idx
            write_positions[kmer_idx] += 1
```

**Pros:**
- Simple, readable code
- Native speed (100-1000x faster than Python loops)
- No sorting needed - direct O(1) writes
- Low memory overhead

**Cons:**
- Adds Numba dependency
- First-call compilation overhead (~1-2 seconds)

**Expected improvement:** Should match or beat original 4-minute time while keeping low memory.

### Option 2: Cython Extension

Similar to Numba but compiled ahead of time:

```cython
# shared_mem_index_cy.pyx
def populate_index_cython(uint64_t[:] seqs, uint64_t[:] shifts, uint64_t mask,
                          uint64_t[:] index_ranges, uint64_t[:] index, int num_kmers):
    cdef uint64_t[:] write_positions = numpy.zeros(len(index_ranges) - 1, dtype=numpy.uint64)
    cdef Py_ssize_t seq_idx, k_pos
    cdef uint64_t seq, kmer_idx, write_pos

    for seq_idx in range(len(seqs)):
        seq = seqs[seq_idx]
        for k_pos in range(num_kmers):
            kmer_idx = (seq >> shifts[k_pos]) & mask
            write_pos = index_ranges[kmer_idx] + write_positions[kmer_idx]
            index[write_pos] = seq_idx
            write_positions[kmer_idx] += 1
```

**Pros:**
- No runtime compilation
- Predictable performance

**Cons:**
- Build complexity (requires Cython compilation)
- Less portable

### Option 3: Optimized Pure NumPy (No New Dependencies)

Use scatter/gather operations to avoid sorting:

```python
def populate_index_numpy_optimized(seqs, shifts, mask, index_ranges, index, num_kmers_per_seq):
    """Direct scatter without sorting."""
    total_kmers = len(index_ranges) - 1
    write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

    # Process smaller chunks to reduce memory
    CHUNK_SIZE = 100_000
    for chunk_start in range(0, len(seqs), CHUNK_SIZE):
        chunk_end = min(chunk_start + CHUNK_SIZE, len(seqs))
        chunk_seqs = seqs[chunk_start:chunk_end]

        # Vectorized k-mer extraction
        chunk_kmers = (chunk_seqs[:, None] >> shifts) & mask  # shape: (chunk, num_kmers)

        # Flatten and process sequentially (unavoidable for correct ordering)
        flat_kmers = chunk_kmers.ravel()
        seq_indices = numpy.repeat(
            numpy.arange(chunk_start, chunk_end, dtype=numpy.uint64),
            num_kmers_per_seq
        )

        # Direct write loop (still Python, but simpler)
        for i in range(len(flat_kmers)):
            kmer_idx = flat_kmers[i]
            write_pos = index_ranges[kmer_idx] + write_positions[kmer_idx]
            index[write_pos] = seq_indices[i]
            write_positions[kmer_idx] += 1
```

**Pros:**
- No new dependencies
- Simpler than current sorting approach

**Cons:**
- Still has Python loop (can't fully vectorize scatter with variable-length buckets)
- Will be slower than Numba/Cython

### Option 4: Hybrid - Keep Original for Small, Two-Pass for Large

```python
def __init__(self, known_bin_seq, kmer_size=12, seq_len=25):
    THRESHOLD = 1_000_000  # Use list-based for < 1M barcodes

    if len(known_bin_seq) < THRESHOLD:
        self._init_list_based(known_bin_seq, kmer_size, seq_len)
    else:
        self._init_two_pass(known_bin_seq, kmer_size, seq_len)
```

**Pros:**
- Best of both worlds
- No new dependencies

**Cons:**
- Code duplication
- Doesn't solve the large barcode case

## Recommendation

**Go with Option 1 (Numba)** because:

1. Numba is already commonly used in bioinformatics
2. Single simple function to JIT compile
3. Expected to match original speed (~4 min) with lower memory
4. Easy to fall back to pure Python if Numba unavailable

## Implementation Steps

1. Add Numba to optional dependencies
2. Create `@njit` decorated function for Pass 2
3. Keep Pass 1 as-is (vectorized counting is fine)
4. Add fallback to pure Python loop if Numba not available
5. Benchmark on 10M and 500M barcodes
6. Clean up temporary arrays after construction

## Memory Analysis

For 10M barcodes, k=12, seq_len=25:
- `seqs`: 10M * 8 bytes = 80 MB
- `index_ranges`: 16M * 8 bytes = 128 MB (shared memory)
- `index`: ~140M * 8 bytes = 1.1 GB (shared memory, this is the k-mer index itself)
- `write_positions`: 16M * 8 bytes = 128 MB (temporary)
- **Total temporary during construction:** ~340 MB
- **Final shared memory:** ~1.3 GB

This is much better than the current approach which creates additional sorting arrays.
