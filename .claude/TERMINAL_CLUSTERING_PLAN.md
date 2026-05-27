# Terminal Position Clustering Improvement Plan

## Problem Statement

Current terminal position clustering in `src/intron_graph.py` has two asymmetric approaches:

1. **3' ends (polyA sites)**: `cluster_polya_positions()` uses greedy clustering with `apa_delta` window
2. **5' ends (TSS)**: `cluster_terminal_positions()` simply takes the extreme (min/max) position

For full-length reads (with both polyA and TSS evidence), we should apply equally robust clustering to both ends.

---

## Current Code Locations

### intron_graph.py

**`cluster_polya_positions()`** (lines 518-550):
- Greedy: picks highest-count position, absorbs neighbors within `apa_delta`
- Returns list of `(position, count, representative_read_ids)` tuples
- Called from `thread_ends()` for 3' positions

**`cluster_terminal_positions()`** (lines 552-561):
- Simple: returns single extreme position (leftmost for 5', rightmost for 3')
- No clustering, no count weighting
- Called from `thread_starts()` for 5' positions

**`thread_ends()`** (lines 563-620):
- Handles 3' terminal positions
- Uses `cluster_polya_positions()` when `polya_read_ids` available
- Falls back to `cluster_terminal_positions()` otherwise

**`thread_starts()`** (lines 622-680):
- Handles 5' terminal positions
- Only uses `cluster_terminal_positions()` (no peak detection)
- Should use equivalent clustering to `thread_ends()`

---

## New Algorithm in polya_position Branch

Branch: `origin/polya_position`
File: `src/polya_position_model.py`

### Key Components:

1. **Peak Detection** (`PositionClusterer.find_peaks()`):
```python
from scipy.signal import find_peaks
peaks, properties = find_peaks(
    smoothed_counts,
    prominence=min_prominence,
    width=1
)
```

2. **Feature Extraction** (`extract_peak_features()`):
- peak_height, prominence, width_half, width_full
- left/right_base positions
- entropy, variance of surrounding region

3. **XGBoost Classification** (`PolyAPositionModel`):
- Binary classifier: real polyA site vs artifact
- Trained on labeled data
- Pre-trained model can be loaded

4. **Integration Point** (`cluster_with_model()`):
- Takes position histogram
- Returns filtered peaks with confidence scores

---

## Implementation Plan

### Phase 1: Unify Terminal Clustering Logic

**Goal**: Make `thread_starts()` use same quality clustering as `thread_ends()`

**Changes to `intron_graph.py`**:

1. Rename `cluster_polya_positions()` → `cluster_terminal_peaks()` (generic)
2. Add parameter `is_5prime: bool` to handle strand-specific logic
3. Update `thread_starts()` to call `cluster_terminal_peaks()` when FL read data available
4. Keep fallback to `cluster_terminal_positions()` for non-FL reads

**Code sketch**:
```python
def cluster_terminal_peaks(self, position_dict, polya_read_ids, intron,
                           is_5prime, apa_delta):
    # Same greedy clustering, works for both ends
    # is_5prime affects assertion direction only
    ...

def thread_starts(self, intron_path, fl_positions, ...):
    if fl_positions:  # Full-length read data available
        clustered = self.cluster_terminal_peaks(
            fl_positions, fl_read_ids, intron,
            is_5prime=True, apa_delta=self.apa_delta
        )
    else:
        # Fallback to extreme position
        ...
```

### Phase 2: Integrate scipy Peak Detection

**Goal**: Replace greedy clustering with proper peak detection

**Changes**:

1. Add scipy dependency check (already in requirements)
2. Create `src/peak_detection.py` module:
   - Extract `find_peaks`-based logic from `polya_position_model.py`
   - Make it work standalone without XGBoost
   - Parameters: `min_prominence`, `min_width`, `min_height`

3. Update `cluster_terminal_peaks()`:
```python
def cluster_terminal_peaks(self, position_dict, ...):
    # Convert position_dict to numpy array
    positions = sorted(position_dict.keys())
    counts = [position_dict[p] for p in positions]

    # Use scipy peak detection
    from .peak_detection import find_position_peaks
    peaks = find_position_peaks(positions, counts,
                                min_prominence=self.min_peak_prominence)

    # Convert peaks back to (position, count, read_ids) format
    ...
```

### Phase 3: Optional XGBoost Classification

**Goal**: Add ML-based peak filtering (optional, for high-quality data)

**Changes**:

1. Keep XGBoost model from `polya_position` branch
2. Add CLI flag: `--use_polya_model PATH`
3. When enabled:
   - Load pre-trained model
   - Score each detected peak
   - Filter by confidence threshold

**Note**: This is optional enhancement. Phase 1-2 provide the core improvement.

---

## Testing Strategy

1. **Unit tests**: Add tests for `cluster_terminal_peaks()` with various distributions
2. **Integration**: Run on SIRV dataset (known transcript ends)
3. **Comparison**: Compare novel transcript 5' ends before/after change
4. **FL read validation**: Check that FL read terminal positions are preserved

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/intron_graph.py` | Unify clustering, add peak detection |
| `src/peak_detection.py` | New module (extract from polya_position branch) |
| `src/polya_position_model.py` | Keep for optional XGBoost (from branch) |
| `tests/test_intron_graph.py` | Uncomment and update tests |

---

## Parameters to Expose

Current hardcoded values to make configurable:

| Parameter | Current | Suggested CLI |
|-----------|---------|---------------|
| `apa_delta` | 50 | `--apa_delta` (exists) |
| min_prominence | N/A | `--min_peak_prominence` |
| min_peak_height | N/A | `--min_peak_height` |

---

## Priority

1. **Phase 1** (High): Symmetric clustering for 5' and 3' - fixes immediate issue
2. **Phase 2** (Medium): scipy peak detection - better algorithm
3. **Phase 3** (Low): XGBoost filtering - optional enhancement

---

## Related Issues

From `.claude/GRAPH_MODEL_ISSUES.md`:
- Issue #2: Suboptimal PolyA/TSS Clustering (this plan addresses it)
- Issue #3: Underutilization of FL read information (Phase 1 addresses it)
