# Phase 7: Architecture Refactoring - Completed

**Date:** 2025-12-31

## Summary

Implemented new domain-driven architecture with clean separation between database persistence, computation, and visualization. Created immutable domain objects, pure computation functions, and repository pattern for data access.

## Problem Statement

See [docs/pathway_architecture_problem_statement.md](docs/pathway_architecture_problem_statement.md) for detailed analysis of the architectural issues addressed.

Key issues solved:
- Pathway classes were tightly coupled to SQLAlchemy models
- No coherent Gene abstraction existed
- "Size" semantics varied by algorithm without explicit typing
- `PathwaySummary` conflated input, output, and metadata
- `RefInfo` was a god object holding all cached data

## New Packages Created

### 1. `app/domain/` - Immutable Domain Objects

| File | Classes | Purpose |
|------|---------|---------|
| [gene.py](app/domain/gene.py) | `Gene` | Immutable gene with entrez_id, symbol, length_bp, mutation_rate |
| [geneset.py](app/domain/geneset.py) | `GeneSet` | Immutable gene collection with algorithm-aware size calculations |
| [pathway.py](app/domain/pathway.py) | `Pathway`, `PathwayMetadata` | Pathway as named gene collection with metadata |
| [mutations.py](app/domain/mutations.py) | `MutationDataset`, `PatientAnalysisData`, `MutationStats`, `RejectedGene` | Patient mutation data structures |
| [enums.py](app/domain/enums.py) | `SizeAlgorithm` | Enum for GENE_COUNT, GENE_LENGTH, BMR_LENGTH |

### 2. `app/computation/` - Pure Computation Functions

| File | Functions/Classes | Purpose |
|------|-------------------|---------|
| [likelihood.py](app/computation/likelihood.py) | `compute_pathway_likelihood()`, `compute_effective_size()`, `compute_confidence_interval()`, `analyze_pathway()` | Pure functions taking primitives, returning results |
| [results.py](app/computation/results.py) | `PathwayAnalysisResult`, `PatientProbabilities` | Immutable result dataclasses |

### 3. `app/repositories/` - Data Access Layer

| File | Classes | Purpose |
|------|---------|---------|
| [gene_repository.py](app/repositories/gene_repository.py) | `GeneRepository` | Load genes from database, caching, background size calculation |
| [pathway_repository.py](app/repositories/pathway_repository.py) | `PathwayRepository` | Load pathways with genes, metadata, size dictionaries |

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    DOMAIN LAYER (Immutable)                     │
│  Gene, GeneSet, Pathway, PathwayMetadata                        │
│  - Loaded once from database at startup                         │
│  - Cached, read-only during analysis                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                 DATA PREPARATION LAYER                          │
│  MutationDataset, PatientAnalysisData                           │
│  - Validates and filters user mutations                         │
│  - Builds patient arrays for computation                        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              PURE COMPUTATION LAYER                             │
│  analyze_pathway(pathway_id, pathway_size, genome_size,         │
│                  n_mutated_array, is_mutated_array)             │
│  - Pure functions, no side effects, no database access          │
│  - Takes primitives, returns PathwayAnalysisResult              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  RESULT OBJECTS LAYER                           │
│  PathwayAnalysisResult, PatientProbabilities                    │
│  - Immutable dataclasses holding computation output             │
│  - Can be serialized (to_dict() method)                         │
└─────────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

### 1. Immutable Domain Objects

All domain classes use `@dataclass(frozen=True)`:
```python
@dataclass(frozen=True)
class Gene:
    entrez_id: int
    symbol: str
    length_bp: int = 0
    mutation_rate: float = 1.0

    @property
    def effective_length(self) -> float:
        return self.length_bp * self.mutation_rate
```

### 2. Algorithm-Aware Size Calculation

`SizeAlgorithm` enum makes semantics explicit:
```python
class SizeAlgorithm(Enum):
    GENE_COUNT = 'gene_count'
    GENE_LENGTH = 'gene_length'
    BMR_LENGTH = 'bmr_length'

# Usage
pathway.get_mutational_target_size(SizeAlgorithm.GENE_COUNT)
gene_set.get_size(SizeAlgorithm.BMR_LENGTH)
```

### 3. Pure Computation Functions

Core computation takes only 5 primitive values:
```python
def analyze_pathway(
    pathway_id: int,
    pathway_size: int,
    genome_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray,
    ...
) -> PathwayAnalysisResult:
```

### 4. Repository Pattern

Database access abstracted from domain objects:
```python
gene_repo = GeneRepository()
genes = gene_repo.get_genes({673, 3845, 7157})
background = gene_repo.get_background_size(SizeAlgorithm.GENE_COUNT)

pathway_repo = PathwayRepository(gene_repo)
pathway = pathway_repo.get_pathway(path_id=123)
```

### 5. Pre-computed Patient Arrays

`PatientAnalysisData` separates shared data from per-pathway data:
```python
# Computed once (shared across all pathways)
patient_data.n_mutated_array  # total mutations per patient

# Computed per pathway
patient_data.get_is_mutated_array(pathway_entrez_ids)
```

## Test Results

```
======================== 62 passed, 5 skipped in 4.80s =========================
```

### New Tests Added

| File | Tests | Coverage |
|------|-------|----------|
| [tests/test_domain.py](tests/test_domain.py) | 32 | Gene, GeneSet, Pathway, SizeAlgorithm, MutationDataset, PatientAnalysisData |
| [tests/test_computation.py](tests/test_computation.py) | 11 | compute_pathway_likelihood, compute_effective_size, analyze_pathway |

## Usage Examples

### Creating Domain Objects

```python
from app.domain import Gene, GeneSet, Pathway, SizeAlgorithm

# Create genes
braf = Gene(entrez_id=673, symbol='BRAF', length_bp=2301, mutation_rate=5.0)
kras = Gene(entrez_id=3845, symbol='KRAS', length_bp=567, mutation_rate=3.0)

# Create gene set
genes = GeneSet([braf, kras])
print(genes.get_size(SizeAlgorithm.GENE_COUNT))  # 2.0
print(genes.get_size(SizeAlgorithm.BMR_LENGTH))  # 13206.0

# Create pathway
pathway = Pathway(path_id=1, name='MAPK_SIGNALING', genes=genes)
print(pathway.gene_count)  # 2
```

### Running Analysis

```python
from app.computation import analyze_pathway
import numpy as np

result = analyze_pathway(
    pathway_id=123,
    pathway_size=100,
    genome_size=18000,
    n_mutated_array=np.array([50, 30, 20, 40]),
    is_mutated_array=np.array([1, 0, 1, 0])
)

print(result.p_value)           # 0.0234
print(result.effect_size)       # 1.45
print(result.patients_covered)  # 2
print(result.to_dict())         # serializable dict
```

### Using Repositories

```python
from app.repositories import GeneRepository, PathwayRepository
from app.domain import SizeAlgorithm

gene_repo = GeneRepository()
pathway_repo = PathwayRepository(gene_repo)

# Load pathway with genes
pathway = pathway_repo.get_pathway(path_id=123)
print(pathway.genes.gene_count)

# Get all pathway sizes efficiently
sizes = pathway_repo.get_pathway_sizes(SizeAlgorithm.GENE_LENGTH)
```

## Migration Path

The new architecture is **additive** - it coexists with existing code:

1. **Phase A** (Complete): Create domain, computation, repository packages
2. **Phase B** (Future): Have `LCalculator` use `analyze_pathway()` internally
3. **Phase C** (Future): Update callers to use new objects
4. **Phase D** (Future): Deprecate `PathwaySummary` classes

## Files Created

| File | Lines | Purpose |
|------|-------|---------|
| `app/domain/__init__.py` | 34 | Package exports |
| `app/domain/gene.py` | 53 | Gene dataclass |
| `app/domain/geneset.py` | 108 | GeneSet dataclass |
| `app/domain/pathway.py` | 130 | Pathway, PathwayMetadata dataclasses |
| `app/domain/mutations.py` | 175 | MutationDataset, PatientAnalysisData |
| `app/domain/enums.py` | 29 | SizeAlgorithm enum |
| `app/computation/__init__.py` | 22 | Package exports |
| `app/computation/likelihood.py` | 310 | Pure computation functions |
| `app/computation/results.py` | 98 | Result dataclasses |
| `app/repositories/__init__.py` | 17 | Package exports |
| `app/repositories/gene_repository.py` | 130 | GeneRepository |
| `app/repositories/pathway_repository.py` | 175 | PathwayRepository |
| `tests/test_domain.py` | 230 | Domain object tests |
| `tests/test_computation.py` | 115 | Computation tests |

## Running Tests

```bash
# Run all tests
docker compose run --rm flask pytest -v

# Run only new architecture tests
docker compose run --rm flask pytest tests/test_domain.py tests/test_computation.py -v

# Run with coverage
docker compose run --rm flask pytest --cov=app.domain --cov=app.computation -v
```

## Next Steps

1. Integrate `analyze_pathway()` into existing `LCalculator`
2. Update `MutationTable` to produce `MutationDataset`
3. Replace `RefInfo` with repository instances
4. Update routes to use `PathwayDisplayData` (to be created)
