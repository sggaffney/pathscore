"""Tests for computation functions."""

import pytest
import numpy as np

from app.computation import (
    compute_pathway_likelihood,
    compute_effective_size,
    analyze_pathway,
    PathwayAnalysisResult,
)


class TestComputePathwayLikelihood:
    """Tests for the likelihood computation function."""

    def test_no_mutations_zero_likelihood(self):
        """Test likelihood when no patients have mutations."""
        genome_size = 18000
        pathway_size = 100
        n_mutated_array = np.array([50, 30, 20], dtype=np.int_)
        is_mutated_array = np.array([0, 0, 0], dtype=np.int_)

        ll, probs = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=True
        )

        # All patients unmutated, so likelihood should be product of
        # (1 - P(mutation)) for each patient
        assert ll < 0  # log likelihood is negative
        assert probs is not None
        assert len(probs) == 3

    def test_all_mutations_likelihood(self):
        """Test likelihood when all patients have mutations."""
        genome_size = 18000
        pathway_size = 1000  # large pathway
        n_mutated_array = np.array([100, 100, 100], dtype=np.int_)
        is_mutated_array = np.array([1, 1, 1], dtype=np.int_)

        ll, _ = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=False
        )

        assert ll < 0  # log likelihood is negative

    def test_mixed_mutations(self):
        """Test likelihood with mixed mutation status."""
        genome_size = 18000
        pathway_size = 100
        n_mutated_array = np.array([50, 30, 10], dtype=np.int_)
        is_mutated_array = np.array([1, 0, 0], dtype=np.int_)

        ll, probs = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=True
        )

        assert ll < 0
        assert len(probs) == 3
        # All probabilities should be between 0 and 1
        assert all(0 <= p <= 1 for p in probs)


class TestComputeEffectiveSize:
    """Tests for effective size calculation."""

    def test_zero_pathway_size(self):
        """Test effective size when actual size is zero."""
        genome_size = 18000
        actual_size = 0
        n_mutated_array = np.array([50, 30, 20], dtype=np.int_)
        is_mutated_array = np.array([0, 0, 0], dtype=np.int_)

        ne, ll = compute_effective_size(
            genome_size, actual_size,
            n_mutated_array, is_mutated_array
        )

        assert ne == 0
        assert ll == 0.0

    def test_all_mutated_large_effective(self):
        """Test effective size when all patients have mutations."""
        genome_size = 18000
        actual_size = 100
        n_mutated_array = np.array([50, 30, 20], dtype=np.int_)
        is_mutated_array = np.array([1, 1, 1], dtype=np.int_)

        ne, ll = compute_effective_size(
            genome_size, actual_size,
            n_mutated_array, is_mutated_array
        )

        # When all patients are mutated, effective size should be large
        assert ne == genome_size

    def test_no_mutations_small_effective(self):
        """Test effective size when no patients have mutations."""
        genome_size = 18000
        actual_size = 100
        n_mutated_array = np.array([50, 30, 20], dtype=np.int_)
        is_mutated_array = np.array([0, 0, 0], dtype=np.int_)

        ne, ll = compute_effective_size(
            genome_size, actual_size,
            n_mutated_array, is_mutated_array
        )

        # When no patients are mutated, effective size should be 0
        assert ne == 0


class TestAnalyzePathway:
    """Tests for the complete pathway analysis function."""

    def test_basic_analysis(self):
        """Test basic pathway analysis returns correct structure."""
        result = analyze_pathway(
            pathway_id=123,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20, 40]),
            is_mutated_array=np.array([1, 0, 1, 0])
        )

        assert isinstance(result, PathwayAnalysisResult)
        assert result.pathway_id == 123
        assert result.n_actual == 100.0
        assert 0 <= result.p_value <= 1
        assert result.patients_covered == 2
        assert result.coverage_fraction == 0.5

    def test_analysis_with_patient_probs(self):
        """Test analysis with patient probabilities."""
        patient_ids = ('A', 'B', 'C', 'D')
        result = analyze_pathway(
            pathway_id=123,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20, 40]),
            is_mutated_array=np.array([1, 0, 1, 0]),
            patient_ids=patient_ids,
            include_patient_probs=True
        )

        assert result.patient_probabilities is not None
        assert len(result.patient_probabilities.probabilities) == 4
        assert result.patient_probabilities.patient_ids == patient_ids

    def test_effect_size(self):
        """Test effect size calculation."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20]),
            is_mutated_array=np.array([1, 0, 0])
        )

        # Effect size should be n_effective / n_actual
        expected_effect = result.n_effective / result.n_actual
        assert abs(result.effect_size - expected_effect) < 0.001

    def test_to_dict(self):
        """Test result serialization to dict."""
        result = analyze_pathway(
            pathway_id=123,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50]),
            is_mutated_array=np.array([1])
        )

        d = result.to_dict()
        assert 'pathway_id' in d
        assert 'p_value' in d
        assert 'effect_size' in d
        assert d['pathway_id'] == 123

    def test_confidence_interval_bounds(self):
        """Test confidence interval has sensible bounds."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20, 40, 60]),
            is_mutated_array=np.array([1, 0, 1, 0, 1])
        )

        # CI should contain the effective size
        assert result.ci_low <= result.n_effective
        assert result.ci_high >= result.n_effective
        # CI should be within genome bounds
        assert result.ci_low >= 0
        assert result.ci_high <= 18000
