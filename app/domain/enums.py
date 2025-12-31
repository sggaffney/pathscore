"""Enumerations for domain objects."""

from enum import Enum


class SizeAlgorithm(Enum):
    """
    Algorithm for calculating pathway mutational target size.

    The "size" of a pathway varies by algorithm:
    - GENE_COUNT: Number of genes in the pathway
    - GENE_LENGTH: Sum of gene lengths in base pairs (CDS length)
    - BMR_LENGTH: Sum of effective lengths (length × mutation rate)
    """
    GENE_COUNT = 'gene_count'
    GENE_LENGTH = 'gene_length'
    BMR_LENGTH = 'bmr_length'

    @classmethod
    def from_string(cls, value: str) -> 'SizeAlgorithm':
        """Convert string to enum, with backwards compatibility."""
        mapping = {
            'gene_count': cls.GENE_COUNT,
            'gene_length': cls.GENE_LENGTH,
            'bmr_length': cls.BMR_LENGTH,
        }
        if value not in mapping:
            raise ValueError(
                f"Unknown algorithm: {value}. "
                f"Must be one of: {list(mapping.keys())}"
            )
        return mapping[value]
