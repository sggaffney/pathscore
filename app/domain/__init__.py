"""
Domain objects for PathScore.

This module contains immutable domain objects that represent the core
concepts in pathway analysis: genes, gene sets, pathways, and mutations.

These objects are separate from SQLAlchemy models and are designed for:
- Clear semantics (a pathway IS a collection of genes)
- Immutability (safe to cache and share)
- Algorithm-aware size calculations
"""

from .gene import Gene
from .geneset import GeneSet
from .pathway import Pathway, PathwayMetadata
from .enums import SizeAlgorithm
from .mutations import (
    MutationDataset,
    MutationStats,
    PatientAnalysisData,
    RejectedGene,
)

__all__ = [
    'Gene',
    'GeneSet',
    'Pathway',
    'PathwayMetadata',
    'SizeAlgorithm',
    'MutationDataset',
    'MutationStats',
    'PatientAnalysisData',
    'RejectedGene',
]
