"""
Repository pattern for data access.

This module provides repository classes that abstract database access
for domain objects. Repositories:
- Load data from the database
- Convert to domain objects
- Provide caching where appropriate
- Hide SQLAlchemy details from calling code

Key principles:
- Repositories return domain objects, not SQLAlchemy models
- Queries are encapsulated within repository methods
- Caching is transparent to callers
"""

from .gene_repository import GeneRepository
from .pathway_repository import PathwayRepository

__all__ = [
    'GeneRepository',
    'PathwayRepository',
]
