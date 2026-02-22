# Code Review: Post-Phase 7 Audit

**Date:** February 22, 2026
**Status:** Complete
**Reviewer:** Claude Opus 4.6

## Background

A review of all code produced by the earlier modernization phases (1–7b, Dec 2025), prompted by CLAUDE.md being out of date with Phase 7/7b changes. After updating CLAUDE.md, a thorough review of key files identified bugs, dead code, architectural concerns, and test quality issues.

## Findings

### Bugs

#### B1. `is_significant` property has illegal signature
**File:** `app/computation/results.py` ~line 80
**Severity:** High (runtime TypeError if called)

```python
@property
def is_significant(self, alpha: float = 0.05) -> bool:
    return self.p_value < alpha
```

Python properties cannot accept parameters beyond `self`. This will raise `TypeError` on any access. Currently unused, but a landmine.

**Fix:** Either remove `@property` (making it a regular method) or remove the `alpha` parameter.

#### B2. Unsafe `eval()` on file content
**File:** `app/get_effective_pathways.py` ~lines 973-978
**Severity:** High (code injection risk)

```python
pway.exclusive = set(eval(vals[11]))
pway.cooccurring = set(eval(vals[12]))
pc_list = list(eval(pc_str.lstrip('struct')))
```

Result files contain MATLAB-style struct strings. `eval()` executes arbitrary Python code. Should use `ast.literal_eval()`.

#### B3. Mutable default arguments
**File:** `app/get_effective_pathways.py`
**Severity:** High (shared state between instances)

`PathwaySummary.__init__` and `PathwayListAssembler.__init__` use `patient_ids=list()` and `ignore_genes=list()` as defaults. All instances share the same list object when defaults are used.

**Fix:** Use `None` defaults with initialization inside the method body.

#### B4. String concatenation with int in `as_string_js()`
**File:** `app/get_effective_pathways.py` ~line 684
**Severity:** High (runtime TypeError)

```python
outstr = "{ id:" + self.path_id + ", "
```

`path_id` is an int, not a string. Would raise `TypeError` at runtime.

**Fix:** `str(self.path_id)`

#### B5. SQLAlchemy 1.x pattern in test fixtures
**File:** `tests/conftest.py` ~line 87
**Severity:** Medium (contradicts project conventions)

```python
general_role = db_session.query(Role).filter_by(name='general').first()
```

Should use `select().where()` per anti-pattern #1 in CLAUDE.md.

### Dead Code

#### D1. Cython import in get_effective_pathways.py
**Lines:** 13–18
**Context:** After Phase 7b deleted `LCalculator`, nothing in this file calls `get_pway_likelihood_cython`. The computation module has its own import.

#### D2. Three unused filter methods on PathwaySummary
**Methods:** `_build_patient_filter_str()`, `_build_ignore_gene_filter_str()`, `_build_expressed_filter_str()`
**Context:** The calling code (in `set_up_from_file`) is commented out. These methods appear to be remnants of an unused filtering feature.

#### D3. Duplicate initialization in PathwaySummaryParsed
Two consecutive `self.gene_coverage = OrderedDict()` assignments in `__init__`.

#### D4. Stale Python 2 comment
Line ~10: `# project_id = raw_input("Project id? ")  # <TODO:sgg> change to input in python3`

### Architecture Observations

#### A1. Domain layer and repositories are speculative infrastructure
`Gene`, `GeneSet`, `Pathway`, `MutationDataset`, `PatientAnalysisData`, `GeneRepository`, `PathwayRepository` — ~900 lines of code with no production callers. The only Phase 7 code wired into the production path is `analyze_pathway()` and `PathwayAnalysisResult`. Everything else exists only in tests and REPL setup, duplicating logic already in `db_lookups.py`.

Not necessarily wrong (deliberate incremental migration), but worth noting that adoption stalled after Phase 7b. The migration plan documented in `progress_2025-12_phase_7.md` lists Phases B–D as "Future."

#### A2. `RefInfo` god object persists
Despite repositories being built as the replacement, `get_effective_pathways.py` still has 16 direct accesses to `ref_info.*` for pathway sizes, names, metadata, and background stats.

#### A3. Text-file serialization is the real migration bottleneck
The write → parse → reconstruct flow (`PathwaySummary` → tab-separated text → `PathwaySummaryParsed`) with magic column indices and `eval()` parsing is what makes legacy classes hard to remove. This pattern, not the computation layer, is what blocks further simplification.

### Test Quality

#### T1. `test_computation.py` doesn't verify mathematical correctness
Tests check shapes and signs (`assert ll < 0`, `assert result.ci_low <= result.n_effective`) but never compare against known expected values. The effect size test is tautological — it computes the expected value from the result itself:

```python
expected_effect = result.n_effective / result.n_actual
assert abs(result.effect_size - expected_effect) < 0.001
```

#### T2. `test_routes.py` has assertions that can't fail
`assert response.status_code in [200, 404]` passes for both success and failure. The `if response.status_code == 200:` pattern means auth failures silently pass.

#### T3. `authenticated_client` fixture in conftest.py is broken and unused
Sets session state outside request context using incorrect Flask-Security-Too patterns. Never called by any test.

#### T4. `test_domain.py` is the strongest test file
Proper edge cases, immutability checks, equality semantics. No issues found.

---

## To-Do List

### Bugs (fix now)

- [x] B1: Fix `is_significant` property signature in `app/computation/results.py` — removed `@property`, kept as regular method
- [x] B2: Replace `eval()` with `ast.literal_eval()` in `app/get_effective_pathways.py`
- [x] B3: Fix mutable default arguments in `PathwaySummary` and `PathwayListAssembler`
- [x] B4: ~~Fix string/int concatenation in `as_string_js()`~~ — false positive: `path_id` is always a string from file parsing in the only construction path
- [x] B5: Fix SQLAlchemy 1.x pattern in `tests/conftest.py`

### Dead code (clean up)

- [x] D1: Remove dead Cython import from `get_effective_pathways.py`
- [x] D2: Remove unused filter methods from `PathwaySummary`
- [x] D3: Remove duplicate `gene_coverage` initialization
- [x] D4: Remove stale Python 2 comment (was in `pathways_js.py`, not `get_effective_pathways.py`)

### Test improvements

- [x] T1: Rewrite `test_computation.py` with hand-calculated expected values using the hypergeometric-like model
- [x] T2: Fix always-passing assertions in `test_routes.py`; moved API auth test out of integration-only class
- [x] T3: Remove broken `authenticated_client` fixture from `conftest.py`

### Documentation

- [x] Update CLAUDE.md to reflect Phase 7/7b changes (done earlier in this session)

## Commits

1. `5553b04` — Fix bugs found in post-Phase 7 code review
2. `e8876d4` — Remove dead code found in post-Phase 7 review
3. `5741748` — Improve test quality from post-Phase 7 review
