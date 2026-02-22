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

    def test_no_mutations_known_value(self):
        """Verify likelihood against hand-calculated value when no patients are mutated.

        For each patient: p = product(1 - pway_size/(G-j) for j in range(n_mut))
        With G=18000, pway_size=100, n_mut=50:
            p = product(1 - 100/(18000-j) for j in range(50))
        Log-likelihood = sum of log(p) for all patients.
        """
        genome_size = 18000
        pathway_size = 100
        n_mutated_array = np.array([50, 30, 20], dtype=np.int_)
        is_mutated_array = np.array([0, 0, 0], dtype=np.int_)

        ll, probs = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=True
        )

        # Compute expected probabilities for each patient
        # Patient 0: 50 mutations, unmutated in pathway
        p0 = np.prod([(1 - 100 / (18000 - j)) for j in range(50)])
        # Patient 1: 30 mutations, unmutated in pathway
        p1 = np.prod([(1 - 100 / (18000 - j)) for j in range(30)])
        # Patient 2: 20 mutations, unmutated in pathway
        p2 = np.prod([(1 - 100 / (18000 - j)) for j in range(20)])

        expected_ll = np.log(p0) + np.log(p1) + np.log(p2)

        np.testing.assert_almost_equal(ll, expected_ll, decimal=6)
        np.testing.assert_almost_equal(probs[0], p0, decimal=10)
        np.testing.assert_almost_equal(probs[1], p1, decimal=10)
        np.testing.assert_almost_equal(probs[2], p2, decimal=10)

    def test_all_mutations_known_value(self):
        """Verify likelihood when all patients have mutations.

        For mutated patients: p = 1 - product(1 - pway_size/(G-j) for j in range(n_mut))
        """
        genome_size = 18000
        pathway_size = 1000
        n_mutated_array = np.array([100, 100, 100], dtype=np.int_)
        is_mutated_array = np.array([1, 1, 1], dtype=np.int_)

        ll, probs = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=True
        )

        # Each patient: p = 1 - product(1 - 1000/(18000-j) for j in range(100))
        p_no_mut = np.prod([(1 - 1000 / (18000 - j)) for j in range(100)])
        expected_p = 1 - p_no_mut
        expected_ll = 3 * np.log(expected_p)

        np.testing.assert_almost_equal(ll, expected_ll, decimal=6)
        for p in probs:
            np.testing.assert_almost_equal(p, expected_p, decimal=10)

    def test_mixed_mutations_known_value(self):
        """Verify likelihood with a mix of mutated and unmutated patients."""
        genome_size = 18000
        pathway_size = 100
        n_mutated_array = np.array([50, 30, 10], dtype=np.int_)
        is_mutated_array = np.array([1, 0, 0], dtype=np.int_)

        ll, probs = compute_pathway_likelihood(
            genome_size, pathway_size,
            n_mutated_array, is_mutated_array,
            return_patient_probs=True
        )

        # Patient 0: 50 mutations, mutated -> p = 1 - product(...)
        p0_no = np.prod([(1 - 100 / (18000 - j)) for j in range(50)])
        p0 = 1 - p0_no
        # Patient 1: 30 mutations, unmutated -> p = product(...)
        p1 = np.prod([(1 - 100 / (18000 - j)) for j in range(30)])
        # Patient 2: 10 mutations, unmutated -> p = product(...)
        p2 = np.prod([(1 - 100 / (18000 - j)) for j in range(10)])

        expected_ll = np.log(p0) + np.log(p1) + np.log(p2)

        np.testing.assert_almost_equal(ll, expected_ll, decimal=6)
        np.testing.assert_almost_equal(probs[0], p0, decimal=10)
        np.testing.assert_almost_equal(probs[1], p1, decimal=10)
        np.testing.assert_almost_equal(probs[2], p2, decimal=10)

    def test_probs_not_returned_when_not_requested(self):
        """Verify probs are None when return_patient_probs=False."""
        ll, probs = compute_pathway_likelihood(
            18000, 100,
            np.array([50], dtype=np.int_),
            np.array([0], dtype=np.int_),
            return_patient_probs=False
        )
        assert probs is None
        assert ll < 0

    def test_single_patient(self):
        """Test with a single patient."""
        ll, probs = compute_pathway_likelihood(
            18000, 100,
            np.array([1], dtype=np.int_),
            np.array([1], dtype=np.int_),
            return_patient_probs=True
        )
        # 1 mutation, pathway size 100, genome 18000
        # p = 1 - (1 - 100/18000) = 100/18000
        expected_p = 100 / 18000
        np.testing.assert_almost_equal(probs[0], expected_p, decimal=10)
        np.testing.assert_almost_equal(ll, np.log(expected_p), decimal=6)


class TestComputeEffectiveSize:
    """Tests for effective size calculation."""

    def test_zero_pathway_size(self):
        """Effective size of a zero-size pathway should be zero."""
        ne, ll = compute_effective_size(
            18000, 0,
            np.array([50, 30, 20], dtype=np.int_),
            np.array([0, 0, 0], dtype=np.int_)
        )
        assert ne == 0
        assert ll == 0.0

    def test_all_mutated_returns_genome_size(self):
        """When all patients are mutated, effective size should be genome_size."""
        ne, ll = compute_effective_size(
            18000, 100,
            np.array([50, 30, 20], dtype=np.int_),
            np.array([1, 1, 1], dtype=np.int_)
        )
        assert ne == 18000

    def test_no_mutations_returns_zero(self):
        """When no patients are mutated, effective size should be 0."""
        ne, ll = compute_effective_size(
            18000, 100,
            np.array([50, 30, 20], dtype=np.int_),
            np.array([0, 0, 0], dtype=np.int_)
        )
        assert ne == 0

    def test_partial_mutations_between_bounds(self):
        """With some mutations, effective size should be between 0 and genome_size."""
        ne, ll = compute_effective_size(
            18000, 100,
            np.array([50, 30, 20, 40, 60], dtype=np.int_),
            np.array([1, 0, 1, 0, 1], dtype=np.int_)
        )
        assert 0 < ne < 18000


class TestAnalyzePathway:
    """Tests for the complete pathway analysis function."""

    def test_basic_analysis_structure(self):
        """Test analysis returns correct structure with correct values."""
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

    def test_effect_size_is_ratio(self):
        """Effect size should be n_effective / n_actual."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20]),
            is_mutated_array=np.array([1, 0, 0])
        )

        assert result.n_actual > 0
        np.testing.assert_almost_equal(
            result.effect_size,
            result.n_effective / result.n_actual,
            decimal=10
        )

    def test_d_statistic_nonnegative(self):
        """D statistic should be >= 0 (likelihood ratio test)."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20, 40, 60]),
            is_mutated_array=np.array([1, 0, 1, 0, 1])
        )
        assert result.d_statistic >= 0

    def test_patient_probs_returned_when_requested(self):
        """Patient probabilities should be returned when requested."""
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

    def test_confidence_interval_contains_estimate(self):
        """CI should contain the effective size estimate."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20, 40, 60]),
            is_mutated_array=np.array([1, 0, 1, 0, 1])
        )

        assert result.ci_low <= result.n_effective
        assert result.ci_high >= result.n_effective
        assert result.ci_low >= 0
        assert result.ci_high <= 18000

    def test_to_dict_completeness(self):
        """to_dict should contain all key fields with correct types."""
        result = analyze_pathway(
            pathway_id=123,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50]),
            is_mutated_array=np.array([1])
        )

        d = result.to_dict()
        assert d['pathway_id'] == 123
        assert isinstance(d['p_value'], float)
        assert isinstance(d['effect_size'], float)
        assert isinstance(d['n_effective'], (int, float))
        # Verify no numpy types leak into dict
        for key, val in d.items():
            if val is not None:
                assert not isinstance(val, np.generic), \
                    f"{key} is numpy type {type(val)}"

    def test_is_significant_method(self):
        """is_significant() should work as a regular method (not property)."""
        result = analyze_pathway(
            pathway_id=1,
            pathway_size=100,
            genome_size=18000,
            n_mutated_array=np.array([50, 30, 20]),
            is_mutated_array=np.array([1, 1, 1])
        )
        # Should be callable with custom alpha
        assert isinstance(result.is_significant(), bool)
        assert isinstance(result.is_significant(alpha=0.01), bool)
