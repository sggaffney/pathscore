# Phase 5: Python Modernization - Completed

**Date:** 2025-12-29

## Summary

Successfully modernized Python code for Python 3.12 compatibility. Updated deprecated pandas patterns and Python 2-era file handling.

## Changes Made

### 1. app/plot/plot_matrix.py - pandas 2.x Updates

**`applymap` → `map`:**
```python
# Before (pandas 1.x)
m_status = m.applymap(lambda x: x is not False)
self.m_status = matrix_df.applymap(lambda x: x is not False)

# After (pandas 2.x)
m_status = m.map(lambda x: x is not False)
self.m_status = matrix_df.map(lambda x: x is not False)
```

### 2. app/uploads.py - Python 3 File Modes

**Deprecated `'rU'` mode:**
```python
# Before (Python 2)
with open(self.temp_path, 'rU') as tempfile:

# After (Python 3)
with open(self.temp_path, 'r') as tempfile:
```

Changed at lines 38 and 86.

### 3. app/pathways_js.py - Python 3 Print

**Print statement → function:**
```python
# Before (Python 2)
print cmd

# After (Python 3)
print(cmd)
```

Changed at lines 50 and 63.

### 4. app/templates/pway/bmr.html - pandas 2.x Template Fix

**`iteritems()` → `items()`:**
```jinja2
{# Before (pandas 1.x) #}
{% for header, val in row.iteritems() %}

{# After (pandas 2.x) #}
{% for header, val in row.items() %}
```

## Files Modified

| File | Changes |
|------|---------|
| [app/plot/plot_matrix.py](app/plot/plot_matrix.py) | `applymap` → `map` (2 occurrences) |
| [app/uploads.py](app/uploads.py) | `'rU'` → `'r'` (2 occurrences) |
| [app/pathways_js.py](app/pathways_js.py) | `print` → `print()` (2 occurrences) |
| [app/templates/pway/bmr.html](app/templates/pway/bmr.html) | `iteritems()` → `items()` |
| [CLAUDE.md](CLAUDE.md) | Updated anti-patterns, removed completed debt items |

## Python 3.12 Compatibility Summary

| Old Pattern | New Pattern | Files Affected |
|-------------|-------------|----------------|
| `print x` | `print(x)` | pathways_js.py |
| `'rU'` file mode | `'r'` | uploads.py, plot_fns.py |
| `applymap()` | `map()` | plot_matrix.py |
| `iteritems()` | `items()` | bmr.html template |
| `pd.np.nan` | `np.nan` | routes.py (Phase 4) |
| Integer division `/` | `//` | plot_fns.py (Phase 4) |

## Verification

- No Python 2 `print` statements remain (except comments)
- No `'rU'` file mode usage remains
- No `applymap()` calls remain
- No `iteritems()` calls remain in templates
- No `xrange`, `basestring`, or `unicode()` usage found

## CLAUDE.md Updates

Added new anti-patterns to document:
- `pd.np` → `np` (pandas 2.x)
- `applymap` → `map` (pandas 2.x)
- `iteritems()` → `items()` (pandas 2.x)
- `'rU'` → `'r'` (Python 3)

Added section documenting Bokeh JavaScript callbacks:
- Selection API patterns
- Q-filtering callback logic
- App-defined JavaScript function interactions

Removed completed technical debt items:
- Bokeh 3.x migration (completed in Phase 4)
- Celery 5.x migration (completed in Phase 3)

## Next Steps

- Phase 6: Testing Infrastructure (pytest fixtures, test suites)
