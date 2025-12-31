"""Result dataclasses for computation output."""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class PatientProbabilities:
    """
    Per-patient probability results.

    Attributes:
        patient_ids: List of patient identifiers (optional)
        probabilities: Probability of observing mutation pattern for each patient
    """
    probabilities: np.ndarray
    patient_ids: Optional[Tuple[str, ...]] = None

    def __post_init__(self):
        """Ensure probabilities is a numpy array."""
        if not isinstance(self.probabilities, np.ndarray):
            object.__setattr__(
                self, 'probabilities',
                np.array(self.probabilities, dtype=np.float64)
            )


@dataclass(frozen=True)
class PathwayAnalysisResult:
    """
    Immutable result from pathway analysis computation.

    This dataclass contains all outputs from analyzing a single pathway,
    including the statistical test results and derived metrics.

    Attributes:
        pathway_id: Identifier for the pathway analyzed
        n_actual: Actual pathway size (gene count, bp, or effective bp)
        n_effective: Maximum likelihood estimate of effective size
        log_likelihood_actual: Log-likelihood at actual size
        log_likelihood_effective: Log-likelihood at effective size
        d_statistic: Likelihood ratio test statistic (D = -2*ll_actual + 2*ll_effective)
        p_value: P-value from chi-squared test on D statistic
        ci_low: Lower bound of 95% confidence interval for n_effective
        ci_high: Upper bound of 95% confidence interval for n_effective
        patients_covered: Number of patients with mutations in this pathway
        coverage_fraction: Fraction of patients with mutations (0.0 to 1.0)
        computation_time_seconds: Time taken for computation (optional)
        patient_probabilities: Per-patient probabilities (optional)
    """
    pathway_id: int
    n_actual: float
    n_effective: float
    log_likelihood_actual: float
    log_likelihood_effective: float
    d_statistic: float
    p_value: float
    ci_low: float
    ci_high: float
    patients_covered: int
    coverage_fraction: float
    computation_time_seconds: Optional[float] = None
    patient_probabilities: Optional[PatientProbabilities] = None

    @property
    def effect_size(self) -> float:
        """
        Effect size: ratio of effective to actual size.

        Values > 1 indicate the pathway is mutated more than expected.
        Values < 1 indicate less mutation than expected.
        """
        if self.n_actual == 0:
            return 0.0
        return self.n_effective / self.n_actual

    @property
    def is_significant(self, alpha: float = 0.05) -> bool:
        """Check if result is statistically significant at given alpha."""
        return self.p_value < alpha

    @property
    def coverage_percentage(self) -> float:
        """Coverage as percentage (0-100)."""
        return self.coverage_fraction * 100

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'pathway_id': self.pathway_id,
            'n_actual': self.n_actual,
            'n_effective': self.n_effective,
            'log_likelihood_actual': self.log_likelihood_actual,
            'log_likelihood_effective': self.log_likelihood_effective,
            'd_statistic': self.d_statistic,
            'p_value': self.p_value,
            'ci_low': self.ci_low,
            'ci_high': self.ci_high,
            'patients_covered': self.patients_covered,
            'coverage_fraction': self.coverage_fraction,
            'effect_size': self.effect_size,
            'computation_time_seconds': self.computation_time_seconds,
        }

    def __repr__(self) -> str:
        return (
            f"PathwayAnalysisResult(pathway_id={self.pathway_id}, "
            f"p_value={self.p_value:.2e}, "
            f"effect_size={self.effect_size:.2f})"
        )
