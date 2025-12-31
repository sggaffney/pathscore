"""Pathway domain object."""

from dataclasses import dataclass
from typing import Optional

from .geneset import GeneSet
from .enums import SizeAlgorithm


@dataclass(frozen=True)
class PathwayMetadata:
    """
    Metadata for a pathway (display information).

    This is separated from the core Pathway class to allow pathways
    to be used in computation without loading display metadata.

    Attributes:
        url: External URL for pathway information
        description: Brief description of the pathway
        contributor: Source/contributor of the pathway definition
    """
    url: str = ''
    description: str = ''
    contributor: str = ''


@dataclass(frozen=True)
class Pathway:
    """
    A named collection of genes with metadata.

    This is the core domain object representing a biological pathway.
    A pathway:
    - Has a unique identifier (path_id)
    - Has a human-readable name
    - Contains a set of genes (GeneSet)
    - Has optional metadata (URL, description, contributor)

    The pathway can calculate its mutational target size using different
    algorithms, delegating to the underlying GeneSet.

    Attributes:
        path_id: Unique pathway identifier
        name: Pathway name (e.g., "KEGG_MAPK_SIGNALING_PATHWAY")
        genes: GeneSet containing the genes in this pathway
        metadata: Optional PathwayMetadata for display information
    """
    path_id: int
    name: str
    genes: GeneSet = None
    metadata: PathwayMetadata = None

    def __post_init__(self):
        """Set defaults for None values (since frozen=True)."""
        if self.genes is None:
            object.__setattr__(self, 'genes', GeneSet())
        if self.metadata is None:
            object.__setattr__(self, 'metadata', PathwayMetadata())

    def get_mutational_target_size(self, algorithm: SizeAlgorithm) -> float:
        """
        Get the mutational target size for this pathway.

        Delegates to the underlying GeneSet.

        Args:
            algorithm: Size calculation algorithm

        Returns:
            Pathway size according to the algorithm
        """
        return self.genes.get_size(algorithm)

    @property
    def gene_count(self) -> int:
        """Number of genes in this pathway."""
        return self.genes.gene_count

    @property
    def nice_name(self) -> str:
        """
        Clean pathway name for display.

        Strips common suffixes like contributor names.
        Example: "KEGG_MAPK_SIGNALING_PATHWAY" -> "MAPK SIGNALING PATHWAY"
        """
        name = self.name

        # Strip common prefixes
        prefixes = ['KEGG_', 'REACTOME_', 'BIOCARTA_', 'PID_', 'ST_']
        for prefix in prefixes:
            if name.startswith(prefix):
                name = name[len(prefix):]
                break

        # Replace underscores with spaces
        name = name.replace('_', ' ')

        return name

    @property
    def url(self) -> str:
        """Shortcut to metadata.url."""
        return self.metadata.url if self.metadata else ''

    @property
    def description(self) -> str:
        """Shortcut to metadata.description."""
        return self.metadata.description if self.metadata else ''

    @property
    def contributor(self) -> str:
        """Shortcut to metadata.contributor."""
        return self.metadata.contributor if self.metadata else ''

    def with_genes(self, genes: GeneSet) -> 'Pathway':
        """Return a new Pathway with different genes (for filtering)."""
        return Pathway(
            path_id=self.path_id,
            name=self.name,
            genes=genes,
            metadata=self.metadata,
        )

    def exclude_genes(self, symbols_to_exclude: set) -> 'Pathway':
        """Return a new Pathway with specified genes removed."""
        filtered = self.genes.exclude_by_symbols(symbols_to_exclude)
        return self.with_genes(filtered)

    def __hash__(self) -> int:
        """Hash by path_id for set operations."""
        return hash(self.path_id)

    def __eq__(self, other: object) -> bool:
        """Equality by path_id."""
        if not isinstance(other, Pathway):
            return NotImplemented
        return self.path_id == other.path_id

    def __repr__(self) -> str:
        return f"Pathway({self.path_id}, '{self.name}', {self.gene_count} genes)"
