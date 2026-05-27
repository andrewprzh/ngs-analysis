# Read Grouping Refactoring (December 2024)

## Overview

Refactored the read grouping system to eliminate deprecated code and use a single, unified table grouper implementation with shared data pattern for memory efficiency.

## Problem Statement

The codebase had multiple table grouper implementations with overlapping functionality:
1. **Old `ReadTableGrouper`** - Loaded table data per instance, no sharing
2. **`MultiColumnReadTableGrouper`** - Returned lists of group IDs, required special handling
3. **`SingleColumnTableGrouper`** - New implementation with shared data pattern

This resulted in:
- Code duplication and maintenance burden
- Special case handling in `MultiReadGrouper.get_group_id()`
- Confusion about which implementation to use

## Solution

Unified all table-based grouping under a single pattern:

### Architecture

```
SharedTableData (loads TSV once)
    ↓
ReadTableGrouper (column_index=0)  ← Single grouper for single column
ReadTableGrouper (column_index=1)  ← Multiple groupers for multi-column
ReadTableGrouper (column_index=2)
```

### Classes After Refactoring

1. **`SharedTableData`** - Loads and stores table data
   - Loads TSV file once
   - Stores: `read_id -> [col1_value, col2_value, ...]`
   - Can be shared across multiple groupers

2. **`ReadTableGrouper`** - The only table grouper (renamed from `SingleColumnTableGrouper`)
   - Constructor: `ReadTableGrouper(shared_data, column_index)`
   - Returns single group ID for one column
   - Multiple instances can share the same `SharedTableData`

3. **Other groupers unchanged:**
   - `AlignmentTagReadGrouper` - Groups by BAM tag (e.g., CB, UB)
   - `ReadIdSplitReadGrouper` - Groups by splitting read ID on delimiter
   - `FileNameGrouper` - Groups by filename
   - `BarcodeSpotGrouper` - Groups by barcode-to-spot mapping

### CLI Behavior

**Single column specification:**
```bash
--read_group file:groups.tsv:0:1
```
Creates:
- 1 `SharedTableData` (loads column 1)
- 1 `ReadTableGrouper` (column_index=0)
- Output: `gene_grouped.file_col1_counts.tsv`

**Multi-column specification:**
```bash
--read_group file:groups.tsv:0:1,2,3
```
Equivalent to:
```bash
--read_group file:groups.tsv:0:1 --read_group file:groups.tsv:0:2 --read_group file:groups.tsv:0:3
```

Creates:
- 1 `SharedTableData` (loads columns 1, 2, 3 once)
- 3 `ReadTableGrouper` instances (column_index=0, 1, 2)
- Output: 3 separate files:
  - `gene_grouped.file_col1_counts.tsv`
  - `gene_grouped.file_col2_counts.tsv`
  - `gene_grouped.file_col3_counts.tsv`

### Code Changes

#### parse_grouping_spec() (src/read_groups.py:364-391)

**Before:**
```python
if ',' in group_col_spec:
    # Multiple columns - use MultiColumnReadTableGrouper
    group_id_column_indices = [int(x) for x in group_col_spec.split(',')]
    return MultiColumnReadTableGrouper(read_group_chr_filename, read_id_column_index,
                                       group_id_column_indices, delim)
else:
    # Single column - use old ReadTableGrouper
    group_id_column_index = int(values[3])
    return ReadTableGrouper(read_group_chr_filename, read_id_column_index,
                          group_id_column_index, delim)
```

**After:**
```python
if ',' in group_col_spec:
    # Multiple columns - create separate groupers sharing the same table data
    group_id_column_indices = [int(x) for x in group_col_spec.split(',')]
    shared_data = SharedTableData(read_group_chr_filename, read_id_column_index,
                                 group_id_column_indices, delim)
    groupers = []
    for i in range(len(group_id_column_indices)):
        groupers.append(ReadTableGrouper(shared_data, i))
    return groupers  # Return list of groupers
else:
    # Single column - use same pattern for consistency
    group_id_column_index = int(values[3])
    shared_data = SharedTableData(read_group_chr_filename, read_id_column_index,
                                 [group_id_column_index], delim)
    return ReadTableGrouper(shared_data, 0)
```

#### create_read_grouper() (src/read_groups.py:439-474)

**Added list flattening:**
```python
for spec in specs:
    spec = spec.strip()
    if spec:
        grouper = parse_grouping_spec(spec, args, sample, chr_id)
        if grouper:
            # parse_grouping_spec can return either a single grouper or a list of groupers
            if isinstance(grouper, list):
                groupers.extend(grouper)
            else:
                groupers.append(grouper)
```

#### MultiReadGrouper.get_group_id() (src/read_groups.py:222-229)

**Before (with special case handling):**
```python
for i, grouper in enumerate(self.groupers):
    gid = grouper.get_group_id(alignment, filename)
    # Handle both single group IDs and lists (for MultiColumnReadTableGrouper)
    if isinstance(gid, list):
        group_ids.extend(gid)
        for g in gid:
            self.read_groups[i].add(g)
    else:
        group_ids.append(gid)
        self.read_groups[i].add(gid)
```

**After (simplified):**
```python
for i, grouper in enumerate(self.groupers):
    gid = grouper.get_group_id(alignment, filename)
    group_ids.append(gid)
    self.read_groups[i].add(gid)
```

## Removed Code

1. **Old `ReadTableGrouper` class** (lines 71-82)
   - Loaded table directly in constructor
   - Each instance had its own copy of table data

2. **`MultiColumnReadTableGrouper` class** (lines 223-248)
   - Returned lists of group IDs
   - Required special handling in `MultiReadGrouper`

3. **`split_read_group_table()` function** (lines 514-536)
   - Obsolete single-threaded table splitting
   - Replaced by `split_read_table_parallel()` from `table_splitter` module

## Retained Helper Functions

- **`load_table()`** - Still used by `BarcodeSpotGrouper` for barcode-to-spot mappings
- **`load_multicolumn_table()`** - Used by `SharedTableData` to load multiple columns

## Testing

All tests updated and passing (60 passed, 2 skipped):
- `TestReadTableGrouper` - Updated to use SharedTableData pattern
- `TestSharedTableData` - Tests multi-column shared data
- `TestParseGroupingSpec` - Tests multi-column spec creates separate groupers
- `TestBarcodeCaller` - 41 tests pass
- `TestUMIFilter` - 6 tests pass, 2 skipped

## Benefits

1. **Code Simplicity**
   - Single table grouper implementation
   - No special case handling for lists
   - Clear, consistent pattern

2. **Memory Efficiency**
   - Table data loaded once and shared
   - Especially important for large TSV files with many columns

3. **Maintainability**
   - Less code to maintain
   - Single source of truth for table-based grouping
   - Clear upgrade path for future improvements

4. **Functional Equivalence**
   - CLI behavior unchanged
   - All existing workflows continue to work
   - Multi-column specs now create separate grouped count tables as intended

## Migration Guide

**Old code (if any external usage):**
```python
# This no longer works
grouper = ReadTableGrouper(file_path, read_col=0, group_col=1, delim='\t')
```

**New code:**
```python
# Use SharedTableData pattern
shared_data = SharedTableData(file_path, read_id_column_index=0,
                             group_id_column_indices=[1], delim='\t')
grouper = ReadTableGrouper(shared_data, column_index=0)
```

## Related Documentation

- `MULTI_GROUP_IMPLEMENTATION.md` - Multi-group support overview
- `BARCODE_UMI_INTEGRATION.md` - Barcode/UMI integration into ReadAssignment
