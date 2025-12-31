# Phase 7b: LCalculator Refactoring - Completed

**Date:** 2025-12-31

## Summary

Deleted the `LCalculator` class and integrated the new `analyze_pathway()` pure function directly into the main analysis loop. This completes the connection between the Phase 7 architecture refactoring and the actual analysis pipeline.

## Problem Statement

The `LCalculator` class (175 lines) was a transitional artifact that:
- Took a `PathwaySummary` object in `__init__`
- Extracted patient data to build numpy arrays
- Stored results in mutable instance attributes (`self.pvalue`, `self.ne`, etc.)
- Mutated the input object (`pway.n_effective = self.ne`)

This pattern conflated input, computation, and output, making the code harder to test and reason about.

## Solution: Clean Break

Rather than wrap `LCalculator` with an adapter pattern, we chose a **clean break**:

1. Only 1 caller instantiated `LCalculator` - trivial to update
2. An adapter would add ~20 lines of mapping code for no benefit
3. Single pattern is clearer than "new way + wrapper for old way"
4. No risk of adapter becoming permanent tech debt

## Changes Made

### 1. Added `p_array` Property to PathwayAnalysisResult

**File:** [app/computation/results.py](app/computation/results.py)

Added convenience accessor to avoid awkward `.patient_probabilities.probabilities`:

```python
@property
def p_array(self) -> Optional[np.ndarray]:
    """Per-patient probabilities as numpy array."""
    if self.patient_probabilities is not None:
        return self.patient_probabilities.probabilities
    return None
```

### 2. Updated `write_pvalue_file` Method

**File:** [app/get_effective_pathways.py](app/get_effective_pathways.py) (lines 481-505)

Changed signature from `write_pvalue_file(self, lcalc, runtime)` to `write_pvalue_file(self, result, runtime)`:

| Old (lcalc attributes) | New (PathwayAnalysisResult) |
|------------------------|----------------------------|
| `lcalc.pway.path_id` | `result.pathway_id` |
| `lcalc.pvalue` | `result.p_value` |
| `lcalc.pway.n_actual` | `result.n_actual` |
| `lcalc.ne` | `result.n_effective` |
| `lcalc.likelihood` | `result.log_likelihood_actual` |
| `lcalc.ne_ll` | `result.log_likelihood_effective` |
| `lcalc.D` | `result.d_statistic` |
| `lcalc.ne_low` | `result.ci_low` |
| `lcalc.ne_high` | `result.ci_high` |
| `sum(lcalc.is_mutated_array)` | `result.patients_covered` |

### 3. Updated Main Analysis Loop

**File:** [app/get_effective_pathways.py](app/get_effective_pathways.py) (lines 835-859)

Before:
```python
for pathway_number in all_path_ids:
    start = timeit.default_timer()
    pway = PathwaySummary(pathway_number, table_list,
                          expressed_table=None,
                          ignore_genes=ignore_genes)
    pway.set_pathway_size(path_size_dict)
    pway.patients = get_patient_list(pathway_number, patient_size_dict,
                                     path_patient_dict)
    current_patient_names = [pa.patient_id for pa in pway.patients]
    lcalc = LCalculator(pway, genome_size)
    lcalc.run()
    runtime = timeit.default_timer() - start
    df_p.loc[pathway_number, current_patient_names] = lcalc.p_array
    basic_writer.write_pvalue_file(lcalc, runtime)
```

After:
```python
for pathway_number in all_path_ids:
    start = timeit.default_timer()

    # Get patient list and build arrays for computation
    patients = get_patient_list(pathway_number, patient_size_dict,
                                path_patient_dict)
    current_patient_names = [p.patient_id for p in patients]
    n_mutated_array = np.array([p.n_mutated for p in patients], dtype=np.int_)
    is_mutated_array = np.array([p.is_mutated for p in patients], dtype=np.int_)

    # Run analysis using pure computation function
    result = analyze_pathway(
        pathway_id=pathway_number,
        pathway_size=path_size_dict[pathway_number],
        genome_size=genome_size,
        n_mutated_array=n_mutated_array,
        is_mutated_array=is_mutated_array,
        include_patient_probs=True
    )

    runtime = timeit.default_timer() - start
    df_p.loc[pathway_number, current_patient_names] = result.p_array
    basic_writer.write_pvalue_file(result, runtime)
```

### 4. Deleted LCalculator Class

Removed entire class (175 lines) from [app/get_effective_pathways.py](app/get_effective_pathways.py).

### 5. Removed Unused Imports

```python
# Removed:
from scipy import stats
from scipy.optimize import minimize_scalar, brentq
```

These were only used by `LCalculator`. The equivalent functionality now lives in [app/computation/likelihood.py](app/computation/likelihood.py).

## Why PathwaySummary Mutations Aren't Needed

The old `LCalculator.run()` set `pway.n_effective` and `pway.p_value`, but investigation showed these values were never used after the loop. The flow is:

1. `write_pvalue_file()` writes results to `.txt` file
2. `PathwayListAssembler.get_ordered_pway_list()` reads file, creates new `PathwaySummaryParsed` objects
3. Downstream code uses those new objects (not the ones from the loop)

So the mutations were effectively dead code.

## Test Results

```
======================== 62 passed, 5 skipped in 5.26s =========================
```

All existing tests continue to pass.

## Code Metrics

| Metric | Value |
|--------|-------|
| Lines deleted | ~175 (LCalculator class) |
| Lines added | ~8 (new loop code) |
| Imports removed | 2 (scipy.stats, scipy.optimize functions) |
| Net reduction | ~165 lines |

## Architecture After Changes

```
┌─────────────────────────────────────────────────────────────────┐
│                    MAIN ANALYSIS LOOP                           │
│  get_patient_list() → Patient objects                           │
│  Build n_mutated_array, is_mutated_array                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              PURE COMPUTATION (app.computation)                 │
│  analyze_pathway(pathway_id, pathway_size, genome_size,         │
│                  n_mutated_array, is_mutated_array)             │
│  → PathwayAnalysisResult                                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    FILE OUTPUT                                  │
│  write_pvalue_file(result, runtime)                             │
│  → Tab-separated values to .txt file                            │
└─────────────────────────────────────────────────────────────────┘
```

## Files Modified

| File | Changes |
|------|---------|
| [app/computation/results.py](app/computation/results.py) | Added `p_array` property |
| [app/get_effective_pathways.py](app/get_effective_pathways.py) | Deleted LCalculator, updated loop and write_pvalue_file, removed imports |

## Next Steps

1. Consider removing `PathwaySummary` class if no longer needed elsewhere
2. Update `PathwayListAssembler` to return domain objects instead of `PathwaySummaryParsed`
3. Migrate remaining callers to use repository pattern from Phase 7
