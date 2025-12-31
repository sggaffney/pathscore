"""Gene domain object."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Gene:
    """
    An immutable gene with its mutational properties.

    This is the core domain object representing a gene. It encapsulates:
    - Identity: entrez_id and symbol
    - Size: length in base pairs
    - Mutation rate: background mutation rate per megabase

    The effective_length property calculates the mutation-rate-weighted
    length, which is used in the BMR_LENGTH algorithm.

    Attributes:
        entrez_id: NCBI Entrez gene ID (unique identifier)
        symbol: HUGO gene symbol (e.g., 'BRAF', 'KRAS')
        length_bp: Gene length in base pairs (CDS length)
        mutation_rate: Background mutation rate per megabase (per_Mb)
                      Default is 1.0 for gene_count/gene_length algorithms
    """
    entrez_id: int
    symbol: str
    length_bp: int = 0
    mutation_rate: float = 1.0

    @property
    def effective_length(self) -> float:
        """
        Mutation-rate weighted length.

        This is length_bp × mutation_rate, used in the BMR_LENGTH algorithm.
        Genes with higher mutation rates have larger effective lengths,
        reflecting their greater chance of acquiring passenger mutations.
        """
        return self.length_bp * self.mutation_rate

    def __hash__(self) -> int:
        """Hash by entrez_id for set operations."""
        return hash(self.entrez_id)

    def __eq__(self, other: object) -> bool:
        """Equality by entrez_id."""
        if not isinstance(other, Gene):
            return NotImplemented
        return self.entrez_id == other.entrez_id

    def __repr__(self) -> str:
        return f"Gene({self.entrez_id}, '{self.symbol}')"
