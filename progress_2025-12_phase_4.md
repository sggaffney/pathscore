# Phase 4: Bokeh 3.x Migration - Completed

**Date:** 2025-12-29

## Summary

Successfully migrated from Bokeh 0.12/2.x to Bokeh 3.x. The main changes involve API renames and JavaScript callback updates for the new selection model.

## Changes Made

### 1. app/plot_fns.py - Plot Configuration and Python 3 Fixes

**Figure parameters renamed:**
```python
# Before (Bokeh 2.x)
plot_config = dict(plot_height=400, plot_width=600, ...)

# After (Bokeh 3.x)
plot_config = dict(height=400, width=600, ...)
```

**Python 3.x file mode:**
```python
# Before
with open(names_path, 'rU') as f:

# After
with open(names_path, 'r') as f:
```

**Python 3 integer division:**
```python
# Before
for group in range(len(vals)/n + 1):

# After
for group in range(len(vals) // n + 1):
```

**Python 3 range comparison:**
```python
# Before
assert js_inds == range(len(all_ids))
plot_inds = range(len(all_ids))

# After
assert js_inds == list(range(len(all_ids)))
plot_inds = list(range(len(all_ids)))
```

### 2. app/pway/routes.py - Figure and Callback Updates

**Figure constructor (lines 289-300):**
```python
# Before
p = figure(plot_width=DIM_COMP_W, plot_height=DIM_COMP_H, ...)

# After
p = figure(width=DIM_COMP_W, height=DIM_COMP_H, ...)
```

**JavaScript callback - Selection API (Q-filtering callback):**
```javascript
// Before (Bokeh 2.x)
var prv_selected = source.selected['1d'].indices;
source.selected['1d'].indices = new_selected;
source.trigger('change');

// After (Bokeh 3.x)
var prv_selected = source.selected.indices;
source.selected.indices = new_selected;
source.change.emit();
```

**pandas 2.x fix:**
```python
# Before (pd.np deprecated)
return pd.np.nan

# After
return np.nan
```

### 3. app/demo/routes.py - Same Changes as pway/routes.py

- Figure `plot_width`/`plot_height` → `width`/`height`
- JavaScript callback selection API updated
- `pd.np.nan` → `np.nan`

## JavaScript Callback Analysis

The Q-filtering callback performs the following operations:

1. **Get previous selection**: Stores indices of previously selected glyphs
2. **Filter by q-value**: Iterates through all pathways, keeping only those where q1 or q2 ≤ cutoff
3. **Update scatter_array**: Global JS array mapping visible glyph indices to original data indices
4. **Preserve selection**: If a previously selected item is still visible after filtering, keep it selected
5. **Emit change**: Trigger Bokeh's reactive update system
6. **Call app function**: `updateIfSelectionChange_afterWait()` - app-specific function that loads pathway images

The `selectPathwaysByGenes()` function (called by cb_inclusion and cb_genes callbacks) is defined in app JavaScript, not in the Python-generated CustomJS.

## Files Modified

| File | Changes |
|------|---------|
| [app/plot_fns.py](app/plot_fns.py) | `height`/`width`, file mode, integer division, range() |
| [app/pway/routes.py](app/pway/routes.py) | `height`/`width`, selection API, `source.change.emit()`, `np.nan` |
| [app/demo/routes.py](app/demo/routes.py) | Same as pway/routes.py |

## Bokeh 3.x API Changes Summary

| Old API (Bokeh 2.x) | New API (Bokeh 3.x) |
|---------------------|---------------------|
| `plot_width` | `width` |
| `plot_height` | `height` |
| `source.selected['1d'].indices` | `source.selected.indices` |
| `source.trigger('change')` | `source.change.emit()` |

## Verification

- Docker container builds successfully
- Flask starts without import errors
- Demo index page loads at http://localhost:5001/demo/
- No Bokeh deprecation warnings in logs

## Note on Full Visualization Testing

Full interactive testing of Bokeh scatter/MDS plots requires:
1. Populated refs database with pathway data
2. Demo project with completed analysis

The code changes are syntactically correct and follow Bokeh 3.x patterns. Full runtime testing can be done once demo data is available.

## Next Steps

- Phase 5: Python Modernization (`applymap` → `map`, remaining file mode updates)
- Phase 6: Testing Infrastructure
