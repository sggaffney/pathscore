"""GeneSet domain object."""

from dataclasses import dataclass, field
from typing import FrozenSet, Iterable, Iterator, Set

from .gene import Gene
from .enums import SizeAlgorithm


@dataclass(frozen=True)
class GeneSet:
    """
    An immutable collection of genes with aggregate properties.

    GeneSet provides a clean abstraction for working with collections of genes.
    It supports different size calculations depending on the algorithm:

    - GENE_COUNT: Number of genes
    - GENE_LENGTH: Sum of gene lengths in base pairs
    - BMR_LENGTH: Sum of effective lengths (mutation-rate weighted)

    This class is immutable and hashable, making it safe to cache and use
    as a dictionary key.

    Attributes:
        genes: Frozen set of Gene objects
    """
    genes: FrozenSet[Gene] = field(default_factory=frozenset)

    def __init__(self, genes: Iterable[Gene] = ()):
        """
        Create a GeneSet from an iterable of Gene objects.

        Args:
            genes: Iterable of Gene objects (list, set, generator, etc.)
        """
        # Use object.__setattr__ because dataclass is frozen
        object.__setattr__(self, 'genes', frozenset(genes))

    @property
    def gene_count(self) -> int:
        """Number of genes in the set."""
        return len(self.genes)

    @property
    def total_length(self) -> int:
        """Sum of gene lengths in base pairs."""
        return sum(g.length_bp for g in self.genes)

    @property
    def total_effective_length(self) -> float:
        """Sum of effective lengths (mutation-rate weighted)."""
        return sum(g.effective_length for g in self.genes)

    def get_size(self, algorithm: SizeAlgorithm) -> float:
        """
        Return size according to specified algorithm.

        Args:
            algorithm: SizeAlgorithm enum value

        Returns:
            Size as int (gene_count) or float (lengths)

        Raises:
            ValueError: If algorithm is not recognized
        """
        if algorithm == SizeAlgorithm.GENE_COUNT:
            return float(self.gene_count)
        elif algorithm == SizeAlgorithm.GENE_LENGTH:
            return float(self.total_length)
        elif algorithm == SizeAlgorithm.BMR_LENGTH:
            return self.total_effective_length
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

    def get_entrez_ids(self) -> Set[int]:
        """Return set of entrez IDs in this gene set."""
        return {g.entrez_id for g in self.genes}

    def get_symbols(self) -> Set[str]:
        """Return set of gene symbols in this gene set."""
        return {g.symbol for g in self.genes}

    def intersection(self, other: 'GeneSet') -> 'GeneSet':
        """Return genes present in both sets."""
        return GeneSet(self.genes & other.genes)

    def union(self, other: 'GeneSet') -> 'GeneSet':
        """Return genes present in either set."""
        return GeneSet(self.genes | other.genes)

    def difference(self, other: 'GeneSet') -> 'GeneSet':
        """Return genes in this set but not in other."""
        return GeneSet(self.genes - other.genes)

    def filter_by_entrez_ids(self, entrez_ids: Set[int]) -> 'GeneSet':
        """Return genes whose entrez_id is in the given set."""
        return GeneSet(g for g in self.genes if g.entrez_id in entrez_ids)

    def exclude_by_symbols(self, symbols: Set[str]) -> 'GeneSet':
        """Return genes whose symbol is NOT in the given set."""
        return GeneSet(g for g in self.genes if g.symbol not in symbols)

    def __len__(self) -> int:
        return len(self.genes)

    def __iter__(self) -> Iterator[Gene]:
        return iter(self.genes)

    def __contains__(self, gene: Gene) -> bool:
        return gene in self.genes

    def __bool__(self) -> bool:
        return len(self.genes) > 0

    def __repr__(self) -> str:
        if len(self.genes) <= 3:
            symbols = ', '.join(sorted(g.symbol for g in self.genes))
            return f"GeneSet({{{symbols}}})"
        else:
            return f"GeneSet({len(self.genes)} genes)"
