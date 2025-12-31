"""Tests for domain objects."""

import pytest
import numpy as np

from app.domain import (
    Gene, GeneSet, Pathway, PathwayMetadata, SizeAlgorithm,
    MutationDataset, MutationStats, PatientAnalysisData, RejectedGene,
)


class TestGene:
    """Tests for Gene domain object."""

    def test_create_gene(self):
        """Test basic gene creation."""
        gene = Gene(
            entrez_id=673,
            symbol='BRAF',
            length_bp=2301,
            mutation_rate=5.0
        )
        assert gene.entrez_id == 673
        assert gene.symbol == 'BRAF'
        assert gene.length_bp == 2301
        assert gene.mutation_rate == 5.0

    def test_gene_effective_length(self):
        """Test effective length calculation."""
        gene = Gene(
            entrez_id=673,
            symbol='BRAF',
            length_bp=1000,
            mutation_rate=2.5
        )
        assert gene.effective_length == 2500.0

    def test_gene_equality(self):
        """Test gene equality is by entrez_id."""
        gene1 = Gene(entrez_id=673, symbol='BRAF')
        gene2 = Gene(entrez_id=673, symbol='BRAF_V2')  # different symbol
        gene3 = Gene(entrez_id=999, symbol='BRAF')     # different id
        assert gene1 == gene2
        assert gene1 != gene3

    def test_gene_hashable(self):
        """Test genes can be used in sets."""
        gene1 = Gene(entrez_id=673, symbol='BRAF')
        gene2 = Gene(entrez_id=3845, symbol='KRAS')
        gene_set = {gene1, gene2}
        assert len(gene_set) == 2
        assert gene1 in gene_set

    def test_gene_immutable(self):
        """Test gene is immutable (frozen dataclass)."""
        gene = Gene(entrez_id=673, symbol='BRAF')
        with pytest.raises(AttributeError):
            gene.symbol = 'CHANGED'


class TestGeneSet:
    """Tests for GeneSet domain object."""

    @pytest.fixture
    def sample_genes(self):
        """Create sample genes for testing."""
        return [
            Gene(entrez_id=673, symbol='BRAF', length_bp=2301, mutation_rate=5.0),
            Gene(entrez_id=3845, symbol='KRAS', length_bp=567, mutation_rate=3.0),
            Gene(entrez_id=7157, symbol='TP53', length_bp=1182, mutation_rate=2.0),
        ]

    def test_create_geneset(self, sample_genes):
        """Test basic GeneSet creation."""
        gs = GeneSet(sample_genes)
        assert gs.gene_count == 3
        assert len(gs) == 3

    def test_geneset_empty(self):
        """Test empty GeneSet."""
        gs = GeneSet()
        assert gs.gene_count == 0
        assert not gs  # bool is False

    def test_geneset_total_length(self, sample_genes):
        """Test total length calculation."""
        gs = GeneSet(sample_genes)
        expected = 2301 + 567 + 1182
        assert gs.total_length == expected

    def test_geneset_total_effective_length(self, sample_genes):
        """Test total effective length calculation."""
        gs = GeneSet(sample_genes)
        expected = (2301 * 5.0) + (567 * 3.0) + (1182 * 2.0)
        assert gs.total_effective_length == expected

    def test_geneset_get_size_by_algorithm(self, sample_genes):
        """Test size calculation by algorithm."""
        gs = GeneSet(sample_genes)

        assert gs.get_size(SizeAlgorithm.GENE_COUNT) == 3.0
        assert gs.get_size(SizeAlgorithm.GENE_LENGTH) == 4050.0
        # BMR_LENGTH = (2301*5.0) + (567*3.0) + (1182*2.0) = 11505 + 1701 + 2364 = 15570
        assert gs.get_size(SizeAlgorithm.BMR_LENGTH) == 15570.0

    def test_geneset_entrez_ids(self, sample_genes):
        """Test getting entrez IDs."""
        gs = GeneSet(sample_genes)
        ids = gs.get_entrez_ids()
        assert ids == {673, 3845, 7157}

    def test_geneset_symbols(self, sample_genes):
        """Test getting gene symbols."""
        gs = GeneSet(sample_genes)
        symbols = gs.get_symbols()
        assert symbols == {'BRAF', 'KRAS', 'TP53'}

    def test_geneset_intersection(self, sample_genes):
        """Test GeneSet intersection."""
        gs1 = GeneSet(sample_genes[:2])  # BRAF, KRAS
        gs2 = GeneSet(sample_genes[1:])  # KRAS, TP53
        result = gs1.intersection(gs2)
        assert result.gene_count == 1
        assert 3845 in result.get_entrez_ids()  # KRAS

    def test_geneset_filter_by_entrez(self, sample_genes):
        """Test filtering by entrez IDs."""
        gs = GeneSet(sample_genes)
        filtered = gs.filter_by_entrez_ids({673, 7157})
        assert filtered.gene_count == 2
        assert 3845 not in filtered.get_entrez_ids()

    def test_geneset_exclude_by_symbols(self, sample_genes):
        """Test excluding genes by symbol."""
        gs = GeneSet(sample_genes)
        filtered = gs.exclude_by_symbols({'BRAF'})
        assert filtered.gene_count == 2
        assert 'BRAF' not in filtered.get_symbols()


class TestPathway:
    """Tests for Pathway domain object."""

    @pytest.fixture
    def sample_pathway(self):
        """Create a sample pathway for testing."""
        genes = GeneSet([
            Gene(entrez_id=673, symbol='BRAF', length_bp=2301),
            Gene(entrez_id=3845, symbol='KRAS', length_bp=567),
        ])
        return Pathway(
            path_id=1,
            name='KEGG_MAPK_SIGNALING_PATHWAY',
            genes=genes,
            metadata=PathwayMetadata(
                url='http://example.com/mapk',
                description='MAPK signaling',
                contributor='KEGG'
            )
        )

    def test_pathway_creation(self, sample_pathway):
        """Test basic pathway creation."""
        assert sample_pathway.path_id == 1
        assert sample_pathway.name == 'KEGG_MAPK_SIGNALING_PATHWAY'
        assert sample_pathway.gene_count == 2

    def test_pathway_metadata_shortcuts(self, sample_pathway):
        """Test metadata property shortcuts."""
        assert sample_pathway.url == 'http://example.com/mapk'
        assert sample_pathway.description == 'MAPK signaling'
        assert sample_pathway.contributor == 'KEGG'

    def test_pathway_nice_name(self, sample_pathway):
        """Test nice_name strips prefix and replaces underscores."""
        assert sample_pathway.nice_name == 'MAPK SIGNALING PATHWAY'

    def test_pathway_get_size(self, sample_pathway):
        """Test getting pathway size by algorithm."""
        assert sample_pathway.get_mutational_target_size(SizeAlgorithm.GENE_COUNT) == 2.0
        assert sample_pathway.get_mutational_target_size(SizeAlgorithm.GENE_LENGTH) == 2868.0

    def test_pathway_exclude_genes(self, sample_pathway):
        """Test creating pathway with genes excluded."""
        filtered = sample_pathway.exclude_genes({'BRAF'})
        assert filtered.gene_count == 1
        assert filtered.path_id == sample_pathway.path_id  # same pathway id

    def test_pathway_equality(self):
        """Test pathway equality is by path_id."""
        p1 = Pathway(path_id=1, name='A')
        p2 = Pathway(path_id=1, name='B')
        p3 = Pathway(path_id=2, name='A')
        assert p1 == p2
        assert p1 != p3


class TestSizeAlgorithm:
    """Tests for SizeAlgorithm enum."""

    def test_from_string(self):
        """Test converting string to enum."""
        assert SizeAlgorithm.from_string('gene_count') == SizeAlgorithm.GENE_COUNT
        assert SizeAlgorithm.from_string('gene_length') == SizeAlgorithm.GENE_LENGTH
        assert SizeAlgorithm.from_string('bmr_length') == SizeAlgorithm.BMR_LENGTH

    def test_from_string_invalid(self):
        """Test invalid string raises error."""
        with pytest.raises(ValueError):
            SizeAlgorithm.from_string('invalid')


class TestMutationDataset:
    """Tests for MutationDataset."""

    @pytest.fixture
    def sample_dataset(self):
        """Create a sample mutation dataset."""
        pairs = frozenset([
            ('patient_A', 673),   # BRAF
            ('patient_A', 3845),  # KRAS
            ('patient_B', 673),   # BRAF
            ('patient_C', 7157),  # TP53
        ])
        return MutationDataset(
            patient_gene_pairs=pairs,
            patient_ids=('patient_A', 'patient_B', 'patient_C')
        )

    def test_patient_mutation_counts(self, sample_dataset):
        """Test getting mutation counts per patient."""
        counts = sample_dataset.get_patient_mutation_counts()
        assert counts['patient_A'] == 2
        assert counts['patient_B'] == 1
        assert counts['patient_C'] == 1

    def test_genes_by_patient(self, sample_dataset):
        """Test getting genes for each patient."""
        genes = sample_dataset.get_genes_by_patient()
        assert genes['patient_A'] == {673, 3845}
        assert genes['patient_B'] == {673}
        assert genes['patient_C'] == {7157}

    def test_patients_by_gene(self, sample_dataset):
        """Test getting patients for each gene."""
        patients = sample_dataset.get_patients_by_gene()
        assert patients[673] == {'patient_A', 'patient_B'}
        assert patients[3845] == {'patient_A'}
        assert patients[7157] == {'patient_C'}

    def test_all_entrez_ids(self, sample_dataset):
        """Test getting all unique gene IDs."""
        ids = sample_dataset.get_all_entrez_ids()
        assert ids == {673, 3845, 7157}


class TestPatientAnalysisData:
    """Tests for PatientAnalysisData."""

    @pytest.fixture
    def sample_analysis_data(self):
        """Create sample analysis data from mutation dataset."""
        pairs = frozenset([
            ('patient_A', 673),   # BRAF
            ('patient_A', 3845),  # KRAS
            ('patient_B', 673),   # BRAF
            ('patient_C', 7157),  # TP53
        ])
        dataset = MutationDataset(
            patient_gene_pairs=pairs,
            patient_ids=('patient_A', 'patient_B', 'patient_C')
        )
        return PatientAnalysisData.from_mutation_dataset(dataset)

    def test_n_mutated_array(self, sample_analysis_data):
        """Test n_mutated_array contains correct counts."""
        expected = np.array([2, 1, 1])  # patient_A, B, C
        np.testing.assert_array_equal(
            sample_analysis_data.n_mutated_array,
            expected
        )

    def test_is_mutated_array_braf(self, sample_analysis_data):
        """Test is_mutated_array for BRAF pathway."""
        pathway_genes = {673}  # BRAF only
        result = sample_analysis_data.get_is_mutated_array(pathway_genes)
        expected = np.array([1, 1, 0])  # A and B have BRAF
        np.testing.assert_array_equal(result, expected)

    def test_is_mutated_array_tp53(self, sample_analysis_data):
        """Test is_mutated_array for TP53 pathway."""
        pathway_genes = {7157}  # TP53 only
        result = sample_analysis_data.get_is_mutated_array(pathway_genes)
        expected = np.array([0, 0, 1])  # only C has TP53
        np.testing.assert_array_equal(result, expected)

    def test_is_mutated_array_multi_gene(self, sample_analysis_data):
        """Test is_mutated_array for pathway with multiple genes."""
        pathway_genes = {673, 3845}  # BRAF and KRAS
        result = sample_analysis_data.get_is_mutated_array(pathway_genes)
        expected = np.array([1, 1, 0])  # A and B have at least one
        np.testing.assert_array_equal(result, expected)

    def test_patients_in_pathway(self, sample_analysis_data):
        """Test getting patients with mutations in pathway."""
        pathway_genes = {673}  # BRAF
        patients = sample_analysis_data.get_patients_in_pathway(pathway_genes)
        assert patients == {'patient_A', 'patient_B'}
