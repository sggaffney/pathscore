"""
Pure computation functions for pathway likelihood analysis.

These functions are stateless and take only primitive values.
They wrap the Cython implementation for performance.
"""

from typing import Optional, Tuple
import time

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar, brentq

from .results import PathwayAnalysisResult, PatientProbabilities

# Import the Cython implementation
try:
    from ..comb_functions import get_pway_likelihood_cython
except ImportError:
    # Fallback to pure Python if Cython not compiled
    def get_pway_likelihood_cython(
        G: int,
        pway_size: int,
        n_patients: int,
        n_mutated_array: np.ndarray,
        is_mutated_array: np.ndarray,
        get_pvals: int = 0
    ):
        """Pure Python fallback for likelihood calculation."""
        prob_array = np.zeros(n_patients, dtype=np.float64)
        log_sum = 0.0

        for i in range(n_patients):
            n_patient = n_mutated_array[i]
            if G - pway_size >= n_patient:
                # Calculate probability of no mutations in pathway
                p_no_mut = 1.0
                for j in range(n_patient):
                    p_no_mut *= (1 - pway_size / (G - j))
                p = 1 - p_no_mut if is_mutated_array[i] > 0 else p_no_mut
            else:
                p = 1.0

            prob_array[i] = p
            log_sum += -np.inf if p == 0 else np.log(p)

        if get_pvals > 0:
            return log_sum, prob_array
        return log_sum


def compute_pathway_likelihood(
    genome_size: int,
    pathway_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray,
    return_patient_probs: bool = False
) -> Tuple[float, Optional[np.ndarray]]:
    """
    Compute log-likelihood of observing mutations at given pathway size.

    This is the core computation: given the background genome size and
    pathway size, calculate how likely it is to observe the mutation
    pattern in the patients.

    Args:
        genome_size: Total number of genes/bp in background genome (G)
        pathway_size: Size of the pathway being tested
        n_mutated_array: Total mutations per patient (across all genes)
        is_mutated_array: Binary indicator of pathway mutation (1=mutated, 0=not)
        return_patient_probs: If True, also return per-patient probabilities

    Returns:
        Tuple of (log_likelihood, patient_probabilities or None)
    """
    n_patients = len(n_mutated_array)

    if return_patient_probs:
        ll, probs = get_pway_likelihood_cython(
            genome_size, pathway_size, n_patients,
            n_mutated_array, is_mutated_array, get_pvals=1
        )
        return ll, probs
    else:
        ll = get_pway_likelihood_cython(
            genome_size, pathway_size, n_patients,
            n_mutated_array, is_mutated_array, get_pvals=0
        )
        return ll, None


def _likelihood_neg(
    pathway_size: float,
    genome_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray
) -> float:
    """Negative log-likelihood for minimization."""
    ll, _ = compute_pathway_likelihood(
        genome_size, int(pathway_size),
        n_mutated_array, is_mutated_array, False
    )
    if ll == -np.inf:
        return np.nan
    return -ll


def compute_effective_size(
    genome_size: int,
    actual_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray
) -> Tuple[int, float]:
    """
    Find the maximum likelihood estimate of pathway size.

    Uses Brent's method to find the pathway size that maximizes
    the likelihood of observing the mutation pattern.

    Args:
        genome_size: Total background genome size
        actual_size: Actual pathway size
        n_mutated_array: Total mutations per patient
        is_mutated_array: Binary pathway mutation indicator

    Returns:
        Tuple of (effective_size, log_likelihood_at_effective)
    """
    # Edge case: pathway has no genes
    if actual_size == 0:
        return 0, 0.0

    # Edge case: all patients have mutations in pathway
    if np.all(is_mutated_array > 0):
        return genome_size, 0.0

    # Edge case: no patients have mutations
    if np.all(is_mutated_array == 0):
        return 0, 0.0

    # Use bounded method to find maximum
    result = minimize_scalar(
        _likelihood_neg,
        args=(genome_size, n_mutated_array, is_mutated_array),
        method='bounded',
        bounds=(1, genome_size)
    )

    # Try both floor and ceiling of result
    candidates = [int(np.floor(result.x)), int(np.ceil(result.x))]
    candidates = [c for c in candidates if 1 <= c <= genome_size]

    best_ne = candidates[0]
    best_ll = -np.inf
    for ne in candidates:
        ll, _ = compute_pathway_likelihood(
            genome_size, ne, n_mutated_array, is_mutated_array, False
        )
        if ll > best_ll:
            best_ll = ll
            best_ne = ne

    return best_ne, best_ll


def _likelihood_for_ci(
    pathway_size: float,
    genome_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray,
    target_ll: float
) -> float:
    """Log-likelihood minus target, for CI boundary finding."""
    ll, _ = compute_pathway_likelihood(
        genome_size, int(pathway_size),
        n_mutated_array, is_mutated_array, False
    )
    # 1.92 is chi2(1, 0.95) / 2 for 95% CI
    return ll - target_ll + 1.92


def compute_confidence_interval(
    genome_size: int,
    actual_size: int,
    effective_size: int,
    ll_effective: float,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray
) -> Tuple[int, int]:
    """
    Compute 95% confidence interval for effective pathway size.

    Uses likelihood ratio to find the range of pathway sizes
    where the log-likelihood is within 1.92 of the maximum
    (corresponding to 95% CI for chi-squared with 1 df).

    Args:
        genome_size: Total background genome size
        actual_size: Actual pathway size
        effective_size: MLE of effective size
        ll_effective: Log-likelihood at effective size
        n_mutated_array: Total mutations per patient
        is_mutated_array: Binary pathway mutation indicator

    Returns:
        Tuple of (ci_low, ci_high)
    """
    ci_low = None
    ci_high = None

    # Edge cases
    if np.all(is_mutated_array > 0):  # everyone mutated
        ci_high = genome_size
    if np.all(is_mutated_array == 0):  # no one mutated
        ci_low = 1
    if effective_size >= genome_size:
        ci_high = genome_size

    # Find lower bound
    if ci_low is None:
        try:
            ci_low = int(np.floor(brentq(
                _likelihood_for_ci,
                1, effective_size,
                args=(genome_size, n_mutated_array, is_mutated_array, ll_effective)
            )))
        except ValueError:
            ci_low = 1

    # Find upper bound
    if ci_high is None:
        try:
            ci_high = int(np.ceil(brentq(
                _likelihood_for_ci,
                effective_size, genome_size - actual_size,
                args=(genome_size, n_mutated_array, is_mutated_array, ll_effective)
            )))
        except ValueError:
            ci_high = genome_size

    return ci_low, ci_high


def analyze_pathway(
    pathway_id: int,
    pathway_size: int,
    genome_size: int,
    n_mutated_array: np.ndarray,
    is_mutated_array: np.ndarray,
    patient_ids: Optional[Tuple[str, ...]] = None,
    include_patient_probs: bool = False
) -> PathwayAnalysisResult:
    """
    Perform complete pathway analysis and return results.

    This is the main entry point for pathway analysis. It:
    1. Computes log-likelihood at actual pathway size
    2. Finds the maximum likelihood effective size
    3. Computes confidence interval
    4. Calculates test statistic and p-value
    5. Calculates coverage statistics

    Args:
        pathway_id: Identifier for the pathway
        pathway_size: Actual pathway size (gene count or bp)
        genome_size: Background genome size
        n_mutated_array: Total mutations per patient
        is_mutated_array: Binary pathway mutation indicator (0/1)
        patient_ids: Optional tuple of patient identifiers
        include_patient_probs: Whether to include per-patient probabilities

    Returns:
        PathwayAnalysisResult with all computed statistics
    """
    start_time = time.perf_counter()

    # Ensure arrays are numpy arrays with correct dtype
    n_mutated_array = np.asarray(n_mutated_array, dtype=np.int_)
    is_mutated_array = np.asarray(is_mutated_array, dtype=np.int_)

    n_patients = len(n_mutated_array)

    # Calculate log-likelihood at actual size
    ll_actual, probs = compute_pathway_likelihood(
        genome_size, pathway_size,
        n_mutated_array, is_mutated_array,
        return_patient_probs=include_patient_probs
    )

    # Find effective size (MLE)
    n_effective, ll_effective = compute_effective_size(
        genome_size, pathway_size,
        n_mutated_array, is_mutated_array
    )

    # Compute confidence interval
    ci_low, ci_high = compute_confidence_interval(
        genome_size, pathway_size, n_effective, ll_effective,
        n_mutated_array, is_mutated_array
    )

    # Calculate test statistic and p-value
    d_statistic = -2 * ll_actual + 2 * ll_effective
    p_value = 1 - stats.chi2.cdf(d_statistic, 1)

    # Calculate coverage
    patients_covered = int(np.sum(is_mutated_array > 0))
    coverage_fraction = patients_covered / n_patients if n_patients > 0 else 0.0

    computation_time = time.perf_counter() - start_time

    # Build patient probabilities if requested
    patient_probs = None
    if include_patient_probs and probs is not None:
        patient_probs = PatientProbabilities(
            probabilities=probs,
            patient_ids=patient_ids
        )

    return PathwayAnalysisResult(
        pathway_id=pathway_id,
        n_actual=float(pathway_size),
        n_effective=float(n_effective),
        log_likelihood_actual=ll_actual,
        log_likelihood_effective=ll_effective,
        d_statistic=d_statistic,
        p_value=p_value,
        ci_low=float(ci_low),
        ci_high=float(ci_high),
        patients_covered=patients_covered,
        coverage_fraction=coverage_fraction,
        computation_time_seconds=computation_time,
        patient_probabilities=patient_probs
    )
