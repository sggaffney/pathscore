"""Mutation data structures for patient analysis."""

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np


@dataclass(frozen=True)
class RejectedGene:
    """
    A gene that was rejected during mutation loading.

    Attributes:
        entrez_id: The entrez ID from the input file
        symbol: The gene symbol from the input file
        reason: Why the gene was rejected
    """
    entrez_id: int
    symbol: str
    reason: str


@dataclass(frozen=True)
class MutationStats:
    """
    Statistics about a loaded mutation dataset.

    Attributes:
        n_patients: Number of unique patients
        n_genes: Number of unique genes (after filtering)
        n_mutations: Total number of mutation records
        n_rejected: Number of rejected gene entries
        n_ignored: Number of entries for genes not in any pathway
    """
    n_patients: int
    n_genes: int
    n_mutations: int
    n_rejected: int = 0
    n_ignored: int = 0

    @property
    def n_loaded(self) -> int:
        """Number of mutations successfully loaded."""
        return self.n_mutations - self.n_rejected - self.n_ignored


@dataclass
class MutationDataset:
    """
    Validated mutations ready for analysis.

    This class holds the validated and filtered mutation data,
    ready to be used for pathway analysis. It provides methods
    to compute the patient arrays needed for computation.

    Attributes:
        patient_gene_pairs: Set of (patient_id, entrez_id) pairs
        patient_ids: Tuple of unique patient IDs (ordered)
        rejected_genes: List of genes that were rejected
        stats: Statistics about the loaded data
    """
    patient_gene_pairs: FrozenSet[Tuple[str, int]]
    patient_ids: Tuple[str, ...]
    rejected_genes: Tuple[RejectedGene, ...] = ()
    stats: MutationStats = None

    def __post_init__(self):
        """Compute stats if not provided."""
        if self.stats is None:
            n_genes = len({pair[1] for pair in self.patient_gene_pairs})
            object.__setattr__(self, 'stats', MutationStats(
                n_patients=len(self.patient_ids),
                n_genes=n_genes,
                n_mutations=len(self.patient_gene_pairs),
                n_rejected=len(self.rejected_genes)
            ))

    def get_patient_mutation_counts(self) -> Dict[str, int]:
        """
        Get total mutations per patient.

        Returns:
            Dict mapping patient_id to number of unique genes mutated
        """
        counts: Dict[str, int] = {pid: 0 for pid in self.patient_ids}
        for patient_id, _ in self.patient_gene_pairs:
            counts[patient_id] = counts.get(patient_id, 0) + 1
        return counts

    def get_genes_by_patient(self) -> Dict[str, Set[int]]:
        """
        Get set of mutated genes for each patient.

        Returns:
            Dict mapping patient_id to set of entrez_ids
        """
        result: Dict[str, Set[int]] = {pid: set() for pid in self.patient_ids}
        for patient_id, entrez_id in self.patient_gene_pairs:
            result[patient_id].add(entrez_id)
        return result

    def get_patients_by_gene(self) -> Dict[int, Set[str]]:
        """
        Get set of patients for each gene.

        Returns:
            Dict mapping entrez_id to set of patient_ids
        """
        result: Dict[int, Set[str]] = {}
        for patient_id, entrez_id in self.patient_gene_pairs:
            if entrez_id not in result:
                result[entrez_id] = set()
            result[entrez_id].add(patient_id)
        return result

    def get_all_entrez_ids(self) -> Set[int]:
        """Get set of all entrez IDs in the dataset."""
        return {pair[1] for pair in self.patient_gene_pairs}


@dataclass(frozen=True)
class PatientAnalysisData:
    """
    Pre-computed arrays for pathway analysis.

    This class holds the patient data in numpy arrays, ready for
    use in the Cython computation. The n_mutated_array is shared
    across all pathways, while is_mutated_array is computed per pathway.

    Attributes:
        patient_ids: Tuple of patient identifiers (ordered)
        n_patients: Number of patients
        n_mutated_array: Total mutations per patient (numpy array)
        _patient_genes: Dict mapping patient to their mutated genes
    """
    patient_ids: Tuple[str, ...]
    n_patients: int
    n_mutated_array: np.ndarray
    _patient_genes: Dict[str, FrozenSet[int]] = field(default_factory=dict)

    @classmethod
    def from_mutation_dataset(cls, dataset: MutationDataset) -> 'PatientAnalysisData':
        """
        Create PatientAnalysisData from a MutationDataset.

        Args:
            dataset: The validated mutation dataset

        Returns:
            PatientAnalysisData ready for analysis
        """
        patient_ids = dataset.patient_ids
        n_patients = len(patient_ids)

        # Build n_mutated_array
        counts = dataset.get_patient_mutation_counts()
        n_mutated_array = np.array(
            [counts.get(pid, 0) for pid in patient_ids],
            dtype=np.int_
        )

        # Build patient_genes map
        genes_by_patient = dataset.get_genes_by_patient()
        patient_genes = {
            pid: frozenset(genes_by_patient.get(pid, set()))
            for pid in patient_ids
        }

        return cls(
            patient_ids=patient_ids,
            n_patients=n_patients,
            n_mutated_array=n_mutated_array,
            _patient_genes=patient_genes
        )

    def get_is_mutated_array(self, pathway_entrez_ids: Set[int]) -> np.ndarray:
        """
        Compute the is_mutated array for a specific pathway.

        Args:
            pathway_entrez_ids: Set of entrez IDs in the pathway

        Returns:
            numpy array where 1 means patient has mutation in pathway
        """
        result = np.zeros(self.n_patients, dtype=np.int_)
        for i, pid in enumerate(self.patient_ids):
            patient_genes = self._patient_genes.get(pid, frozenset())
            if patient_genes & pathway_entrez_ids:  # intersection
                result[i] = 1
        return result

    def get_patients_in_pathway(self, pathway_entrez_ids: Set[int]) -> Set[str]:
        """
        Get set of patients with mutations in the pathway.

        Args:
            pathway_entrez_ids: Set of entrez IDs in the pathway

        Returns:
            Set of patient IDs with at least one mutation in pathway
        """
        return {
            pid for pid, genes in self._patient_genes.items()
            if genes & pathway_entrez_ids
        }
