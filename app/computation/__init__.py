"""
Pure computation functions for pathway analysis.

This module contains stateless, pure functions that perform the core
statistical computations. They take primitive values (ints, floats,
numpy arrays) and return results.

Key design principles:
- Pure functions: no side effects, no database access
- Takes primitives: not pathway objects or SQLAlchemy models
- Returns dataclasses: structured, immutable results
"""

from .likelihood import (
    compute_pathway_likelihood,
    compute_effective_size,
    compute_confidence_interval,
    analyze_pathway,
)
from .results import (
    PathwayAnalysisResult,
    PatientProbabilities,
)

__all__ = [
    'compute_pathway_likelihood',
    'compute_effective_size',
    'compute_confidence_interval',
    'analyze_pathway',
    'PathwayAnalysisResult',
    'PatientProbabilities',
]
