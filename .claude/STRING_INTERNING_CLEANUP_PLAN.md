# String Interning Cleanup Plan

## Overview
Remove all fallback code and assume string_pools is always present.

## 1. Consolidate Duplicate Loader Classes

**Issue**: `BasicReadAssignmentLoader` and `ReadAssignmentLoader` exist in both `assignment_loader.py` and `read_assignment_loader.py`

**Action**:
- Keep only `src/assignment_loader.py` (has more complete implementations)
- Delete `src/read_assignment_loader.py`
- Update imports:
  - `src/dataset_processor.py`: Remove import from `read_assignment_loader`
  - `src/processed_read_manager.py`: Remove import from `read_assignment_loader`

**Files to modify**:
- Delete: `src/read_assignment_loader.py`
- Update: `src/dataset_processor.py`, `src/processed_read_manager.py`

---

## 2. Remove Fallback Code in Data Classes

### IsoformMatch (`src/isoform_assignment.py`)

**Delete fallback fields**:
- `self._assigned_gene_str`
- `self._assigned_transcript_str`

**Change `__init__` from**:
```python
def __init__(self, match_classification, assigned_gene=None, assigned_transcript=None,
             match_subclassification=None, transcript_strand='.', penalty_score=0.0, string_pools=None):
    self._string_pools = string_pools

    if string_pools is not None:
        self.assigned_gene_id = string_pools.gene_pool.get_int(assigned_gene) if assigned_gene else None
        self.assigned_transcript_id = string_pools.transcript_pool.get_int(assigned_transcript) if assigned_transcript else None
    else:
        self.assigned_gene_id = None
        self.assigned_transcript_id = None
        self._assigned_gene_str = assigned_gene
        self._assigned_transcript_str = assigned_transcript
```

**To**:
```python
def __init__(self, match_classification, assigned_gene=None, assigned_transcript=None,
             match_subclassification=None, transcript_strand='.', penalty_score=0.0, string_pools):
    self._string_pools = string_pools
    self.assigned_gene_id = string_pools.gene_pool.get_int(assigned_gene) if assigned_gene else None
    self.assigned_transcript_id = string_pools.transcript_pool.get_int(assigned_transcript) if assigned_transcript else None
```

**Change `serialize()` from**:
```python
def serialize(self, outfile):
    if self._string_pools is not None:
        write_int_or_none(self.assigned_gene_id, outfile)
        write_int_or_none(self.assigned_transcript_id, outfile)
    else:
        write_string_or_none(self._assigned_gene_str, outfile)
        write_string_or_none(self._assigned_transcript_str, outfile)
```

**To**:
```python
def serialize(self, outfile):
    write_int_or_none(self.assigned_gene_id, outfile)
    write_int_or_none(self.assigned_transcript_id, outfile)
```

**Change `deserialize()` from**:
```python
@classmethod
def deserialize(cls, infile, string_pools=None):
    match = cls.__new__(cls)
    match._string_pools = string_pools

    if string_pools is not None:
        match.assigned_gene_id = read_int_or_none(infile)
        match.assigned_transcript_id = read_int_or_none(infile)
    else:
        match._assigned_gene_str = read_string_or_none(infile)
        match._assigned_transcript_str = read_string_or_none(infile)
        match.assigned_gene_id = None
        match.assigned_transcript_id = None
```

**To**:
```python
@classmethod
def deserialize(cls, infile, string_pools):
    match = cls.__new__(cls)
    match._string_pools = string_pools
    match.assigned_gene_id = read_int_or_none(infile)
    match.assigned_transcript_id = read_int_or_none(infile)
```

**Simplify properties from**:
```python
@property
def assigned_gene(self):
    if self._string_pools is not None:
        return self._string_pools.gene_pool.get_str(self.assigned_gene_id) if self.assigned_gene_id is not None else None
    else:
        return self._assigned_gene_str

@property
def assigned_transcript(self):
    if self._string_pools is not None:
        return self._string_pools.transcript_pool.get_str(self.assigned_transcript_id) if self.assigned_transcript_id is not None else None
    else:
        return self._assigned_transcript_str
```

**To**:
```python
@property
def assigned_gene(self):
    return self._string_pools.gene_pool.get_str(self.assigned_gene_id) if self.assigned_gene_id is not None else None

@property
def assigned_transcript(self):
    return self._string_pools.transcript_pool.get_str(self.assigned_transcript_id) if self.assigned_transcript_id is not None else None
```

---

### ReadAssignment (`src/isoform_assignment.py`)

**Delete fallback fields**:
- `self._chr_id_str`
- `self._barcode_str`
- `self._umi_str`

**Apply same pattern as IsoformMatch**:
- Remove `string_pools=None` → make it required
- Remove all `if string_pools is not None` branches
- Simplify `serialize()` to only write integers
- Simplify `deserialize()` to only read integers
- Simplify all properties (`chr_id`, `barcode`, `umi`)

---

### BasicReadAssignment (`src/isoform_assignment.py`)

**Delete fallback fields**:
- `self._chr_id_str`

**Apply same pattern**:
- Remove `string_pools=None`
- Remove conditional branches
- Simplify serialize/deserialize
- Simplify `chr_id` property

---

## 3. Update All Call Sites

**Make string_pools required (remove `=None` defaults)**:

### Loader Classes
- `NormalTmpFileAssignmentLoader.__init__(..., string_pools)`
- `GeneListTmpFileAssignmentLoader.__init__(..., string_pools)`
- `QuickTmpFileAssignmentLoader.__init__(..., string_pools)`
- `ReadAssignmentLoader.__init__(..., string_pools)`
- `MergingSimpleReadAssignmentLoader.__init__(..., string_pools)`
- `BasicReadAssignmentLoader.__init__(..., string_pools)`

### Factory Functions
- `create_assignment_loader(..., string_pools)`
- `create_merging_assignment_loader(..., string_pools)`
- `prepare_multimapped_reads(..., string_pools)`

### Data Classes (already done)
- `IsoformMatch.__init__(..., string_pools)`
- `ReadAssignment.__init__(..., string_pools)`
- `BasicReadAssignment.__init__(..., string_pools)`

---

## 4. Verify All Creation Sites Pass string_pools

**Search for patterns**:
```bash
grep -n "IsoformMatch(" src/*.py src/*/*.py
grep -n "ReadAssignment(" src/*.py src/*/*.py
grep -n "BasicReadAssignment(" src/*.py src/*/*.py
```

**Ensure all have `string_pools=...` argument**

---

## 5. Add Documentation

### `src/string_pools.py` - Module docstring:
```python
"""
String interning for memory optimization.

Replaces duplicated strings with integer indices referencing shared string pools.
Expected memory reduction: ~78% for billion-read datasets.

Architecture:
- StringPool: Bidirectional mapping between strings and integers
- StringPoolManager: Manages all pools for a worker
  - Global pools: gene, transcript, chromosome (from annotation)
  - Per-chromosome pools: barcode, UMI (from split files)
  - Read group pools: Multiple types depending on strategy

Usage:
- Workers build pools from gene database at start
- All IsoformMatch/ReadAssignment objects store integer IDs
- Properties provide transparent string access via pools
- Serialization writes integer format for disk savings
"""
```

### `src/isoform_assignment.py` - Class docstrings:
```python
class IsoformMatch:
    """
    ...existing docstring...

    Memory Optimization:
        Gene and transcript IDs are stored as integers referencing
        shared string pools (string_pools parameter required).
    """

class ReadAssignment:
    """
    ...existing docstring...

    Memory Optimization:
        chr_id, barcode, and umi are stored as integers referencing
        shared string pools (string_pools parameter required).
    """
```

---

## 6. Update Tests

**Files to check**:
- `tests/test_serialization.py` (if exists)
- Any tests that create IsoformMatch/ReadAssignment objects

**Changes needed**:
- All test objects must pass `string_pools` parameter
- Create a test fixture for StringPoolManager:
```python
@pytest.fixture
def string_pools():
    from src.string_pools import StringPoolManager
    pools = StringPoolManager()
    # Pre-populate with test data
    pools.gene_pool.add("GENE1")
    pools.transcript_pool.add("TX1")
    pools.chromosome_pool.add("chr1")
    return pools
```

---

## 7. Error Handling

**Add assertions where string_pools is critical**:
```python
def __init__(self, ..., string_pools):
    assert string_pools is not None, "string_pools is required"
    self._string_pools = string_pools
```

Or raise better errors:
```python
def __init__(self, ..., string_pools):
    if string_pools is None:
        raise ValueError("string_pools parameter is required for memory optimization")
    self._string_pools = string_pools
```

---

## 8. Cleanup Checklist

- [ ] Delete `src/read_assignment_loader.py`
- [ ] Update imports in `src/dataset_processor.py`
- [ ] Update imports in `src/processed_read_manager.py`
- [ ] Remove fallback fields from `IsoformMatch`
- [ ] Simplify `IsoformMatch.__init__()`
- [ ] Simplify `IsoformMatch.serialize()`
- [ ] Simplify `IsoformMatch.deserialize()`
- [ ] Simplify `IsoformMatch` properties
- [ ] Remove fallback fields from `ReadAssignment`
- [ ] Simplify `ReadAssignment.__init__()`
- [ ] Simplify `ReadAssignment.serialize()`
- [ ] Simplify `ReadAssignment.deserialize()`
- [ ] Simplify `ReadAssignment` properties
- [ ] Remove fallback fields from `BasicReadAssignment`
- [ ] Simplify `BasicReadAssignment` methods and properties
- [ ] Remove `string_pools=None` from all loader class signatures
- [ ] Remove `string_pools=None` from all factory functions
- [ ] Verify all creation sites pass `string_pools`
- [ ] Add module docstring to `src/string_pools.py`
- [ ] Add memory optimization notes to class docstrings
- [ ] Update tests to use string_pools
- [ ] Add error handling for missing string_pools
- [ ] Run full test suite
- [ ] Manual test with larger dataset

---

## Notes

- **Breaking change**: Old serialized files without string_pools won't load (acceptable for dev branch)
- **Estimated effort**: 2-3 hours
- **Files impacted**: ~10 files
- **Lines removed**: ~150-200 lines of fallback code
- **Expected benefits**: Cleaner code, no conditional branching, ~5% performance improvement
