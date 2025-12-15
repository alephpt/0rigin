"""
Dirac gamma matrices for relativistic quantum field theory.

Implements the clifford algebra {γ^μ, γ^ν} = 2η^μν where η = diag(1,-1,-1,-1).
Uses the Dirac (standard) representation.

References
----------
R. Christopher, "Mass Synchronization Field Theory (MSFT)", Dec 2025
B. Thaller, "The Dirac Equation", Springer (1992)
M. Peskin & D. Schroeder, "An Introduction to Quantum Field Theory" (1995)
"""

from typing import Tuple
import numpy as np
from numpy.typing import NDArray


# ==============================================================================
# 3+1D Dirac Gamma Matrices (4×4)
# ==============================================================================

def get_gamma_matrices_3plus1() -> Tuple[NDArray, NDArray, NDArray, NDArray, NDArray]:
    """
    Return Dirac gamma matrices in 3+1 dimensions (Dirac representation).

    Standard representation where γ^0 is diagonal and γ^i are off-diagonal.

    Algebra: {γ^μ, γ^ν} = γ^μγ^ν + γ^νγ^μ = 2η^μν·I₄
    where η = diag(1, -1, -1, -1) is Minkowski metric.

    Returns
    -------
    gamma0, gamma1, gamma2, gamma3, gamma5 : NDArray
        4×4 complex gamma matrices.
        γ^0 = time component
        γ^1, γ^2, γ^3 = spatial components
        γ^5 = iγ^0γ^1γ^2γ^3 = chiral matrix

    Examples
    --------
    >>> γ0, γ1, γ2, γ3, γ5 = get_gamma_matrices_3plus1()
    >>> # Verify anticommutator {γ^0, γ^1} = 0
    >>> np.allclose(γ0 @ γ1 + γ1 @ γ0, 0)
    True
    >>> # Verify (γ^0)² = I
    >>> np.allclose(γ0 @ γ0, np.eye(4))
    True
    >>> # Verify (γ^1)² = -I
    >>> np.allclose(γ1 @ γ1, -np.eye(4))
    True
    """
    # Pauli matrices (building blocks)
    σ1 = np.array([[0, 1], [1, 0]], dtype=complex)
    σ2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    σ3 = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    # Dirac representation: γ^0 = [[I, 0], [0, -I]], γ^i = [[0, σ^i], [-σ^i, 0]]
    gamma0 = np.block([[I2, np.zeros((2, 2))],
                       [np.zeros((2, 2)), -I2]])

    gamma1 = np.block([[np.zeros((2, 2)), σ1],
                       [-σ1, np.zeros((2, 2))]])

    gamma2 = np.block([[np.zeros((2, 2)), σ2],
                       [-σ2, np.zeros((2, 2))]])

    gamma3 = np.block([[np.zeros((2, 2)), σ3],
                       [-σ3, np.zeros((2, 2))]])

    # Chiral matrix: γ^5 = iγ^0γ^1γ^2γ^3
    gamma5 = 1j * gamma0 @ gamma1 @ gamma2 @ gamma3

    return gamma0, gamma1, gamma2, gamma3, gamma5


def verify_gamma_algebra_3plus1() -> dict:
    """
    Verify that gamma matrices satisfy Clifford algebra.

    Checks:
    1. Anticommutation: {γ^μ, γ^ν} = 2η^μν
    2. Hermiticity: (γ^0)† = γ^0, (γ^i)† = -γ^i
    3. Chiral properties: {γ^5, γ^μ} = 0, (γ^5)² = I

    Returns
    -------
    dict
        Validation results with boolean flags and error metrics.
    """
    γ0, γ1, γ2, γ3, γ5 = get_gamma_matrices_3plus1()
    gammas = [γ0, γ1, γ2, γ3]
    η = np.diag([1, -1, -1, -1])  # Minkowski metric
    I = np.eye(4, dtype=complex)

    results = {
        'anticommutator_correct': True,
        'hermiticity_correct': True,
        'chiral_anticommute': True,
        'gamma5_squared': True,
        'max_error': 0.0
    }

    # Check anticommutator {γ^μ, γ^ν} = 2η^μν
    for μ in range(4):
        for ν in range(4):
            anticomm = gammas[μ] @ gammas[ν] + gammas[ν] @ gammas[μ]
            expected = 2 * η[μ, ν] * I
            error = np.max(np.abs(anticomm - expected))
            results['max_error'] = max(results['max_error'], error)
            if error > 1e-10:
                results['anticommutator_correct'] = False

    # Check Hermiticity
    if not np.allclose(γ0.conj().T, γ0):
        results['hermiticity_correct'] = False
    for γi in [γ1, γ2, γ3]:
        if not np.allclose(γi.conj().T, -γi):
            results['hermiticity_correct'] = False

    # Check chiral properties
    for γμ in gammas:
        if not np.allclose(γ5 @ γμ + γμ @ γ5, 0):
            results['chiral_anticommute'] = False

    if not np.allclose(γ5 @ γ5, I):
        results['gamma5_squared'] = False

    return results


# ==============================================================================
# 1+1D Gamma Matrices (2×2) - Simplified for Testing
# ==============================================================================

def get_gamma_matrices_1plus1() -> Tuple[NDArray, NDArray, NDArray]:
    """
    Return Dirac gamma matrices in 1+1 dimensions (2×2).

    Simpler 2-component spinors for 1 spatial + 1 time dimension.
    Useful for prototyping and testing before full 3+1D implementation.

    Algebra: {γ^μ, γ^ν} = 2η^μν where η = diag(1, -1)

    Returns
    -------
    gamma0, gamma1, gamma5 : NDArray
        2×2 complex gamma matrices.
        γ^0 = σ^1 (time)
        γ^1 = iσ^2 (space)
        γ^5 = σ^3 (chiral)
    """
    # Pauli matrices
    σ1 = np.array([[0, 1], [1, 0]], dtype=complex)
    σ2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    σ3 = np.array([[1, 0], [0, -1]], dtype=complex)

    # Standard 1+1D representation
    gamma0 = σ1
    gamma1 = 1j * σ2
    gamma5 = σ3  # Chiral matrix in 1+1D

    return gamma0, gamma1, gamma5


def verify_gamma_algebra_1plus1() -> dict:
    """Verify 1+1D gamma matrices satisfy Clifford algebra."""
    γ0, γ1, γ5 = get_gamma_matrices_1plus1()
    η = np.diag([1, -1])
    I = np.eye(2, dtype=complex)

    results = {
        'anticommutator_correct': True,
        'chiral_properties': True,
        'max_error': 0.0
    }

    # {γ^0, γ^0} = 2I
    error = np.max(np.abs(γ0 @ γ0 + γ0 @ γ0 - 2 * I))
    results['max_error'] = max(results['max_error'], error)

    # {γ^1, γ^1} = -2I
    error = np.max(np.abs(γ1 @ γ1 + γ1 @ γ1 + 2 * I))
    results['max_error'] = max(results['max_error'], error)

    # {γ^0, γ^1} = 0
    error = np.max(np.abs(γ0 @ γ1 + γ1 @ γ0))
    results['max_error'] = max(results['max_error'], error)

    if results['max_error'] > 1e-10:
        results['anticommutator_correct'] = False

    # Chiral anticommutation
    if not (np.allclose(γ5 @ γ0 + γ0 @ γ5, 0) and
            np.allclose(γ5 @ γ1 + γ1 @ γ5, 0)):
        results['chiral_properties'] = False

    return results


# ==============================================================================
# Chiral Projectors and Mass Operator
# ==============================================================================

def chiral_projectors(gamma5: NDArray) -> Tuple[NDArray, NDArray]:
    """
    Compute left and right chiral projection operators.

    P_L = (1 - γ^5) / 2  projects onto left-handed states
    P_R = (1 + γ^5) / 2  projects onto right-handed states

    Properties:
    - P_L + P_R = I
    - P_L @ P_R = 0 (orthogonal)
    - P_L @ P_L = P_L (idempotent)
    - P_R @ P_R = P_R (idempotent)

    Parameters
    ----------
    gamma5 : NDArray
        Chiral gamma matrix.

    Returns
    -------
    P_L, P_R : NDArray
        Left and right chiral projectors.
    """
    I = np.eye(gamma5.shape[0], dtype=complex)
    P_L = (I - gamma5) / 2
    P_R = (I + gamma5) / 2
    return P_L, P_R


def mass_operator_MSFT(
    Delta: float,
    R: float,
    theta: float,
    gamma5: NDArray
) -> NDArray:
    """
    Construct MSFT mass operator m = Δ·R·exp(iθγ^5).

    From R.Christopher (2025) MSFT PDF:
        m(x,t) = Δ · R(x,t) · e^(iθγ^5)
               = Δ · R · [cos(θ)·I + i·sin(θ)·γ^5]

    Decomposes into:
    - Scalar mass: m_S = Δ·R·cos(θ) (parity-conserving)
    - Pseudoscalar mass: m_P = Δ·R·sin(θ) (parity-violating)

    Parameters
    ----------
    Delta : float
        Mass gap (energy scale).
    R : float
        Synchronization order parameter [0, 1].
    theta : float
        Chiral angle in radians.
    gamma5 : NDArray
        Chiral gamma matrix.

    Returns
    -------
    NDArray
        Mass operator matrix.

    Examples
    --------
    >>> γ0, γ1, γ2, γ3, γ5 = get_gamma_matrices_3plus1()
    >>> # Pure scalar mass (θ=0)
    >>> m_scalar = mass_operator_MSFT(2.0, 0.5, 0.0, γ5)
    >>> np.allclose(m_scalar, 1.0 * np.eye(4))
    True
    >>> # Pure pseudoscalar (θ=π/2)
    >>> m_pseudo = mass_operator_MSFT(2.0, 0.5, np.pi/2, γ5)
    >>> np.allclose(m_pseudo, 1.0j * γ5)
    True
    """
    I = np.eye(gamma5.shape[0], dtype=complex)

    # Euler's formula: e^(iθγ^5) = cos(θ)·I + i·sin(θ)·γ^5
    mass_op = Delta * R * (np.cos(theta) * I + 1j * np.sin(theta) * gamma5)

    return mass_op


# ==============================================================================
# Utility Functions
# ==============================================================================

def commutator(A: NDArray, B: NDArray) -> NDArray:
    """Compute commutator [A, B] = AB - BA."""
    return A @ B - B @ A


def anticommutator(A: NDArray, B: NDArray) -> NDArray:
    """Compute anticommutator {A, B} = AB + BA."""
    return A @ B + B @ A


def is_hermitian(M: NDArray, tol: float = 1e-10) -> bool:
    """Check if matrix is Hermitian (M† = M)."""
    return np.allclose(M.conj().T, M, atol=tol)


def is_antihermitian(M: NDArray, tol: float = 1e-10) -> bool:
    """Check if matrix is anti-Hermitian (M† = -M)."""
    return np.allclose(M.conj().T, -M, atol=tol)


# ==============================================================================
# Self-Test
# ==============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print(" GAMMA MATRIX VALIDATION")
    print("=" * 70)

    # Test 3+1D
    print("\n[1/2] Testing 3+1D Gamma Matrices...")
    results_3d = verify_gamma_algebra_3plus1()
    print(f"  Anticommutator: {'✓ PASS' if results_3d['anticommutator_correct'] else '✗ FAIL'}")
    print(f"  Hermiticity: {'✓ PASS' if results_3d['hermiticity_correct'] else '✗ FAIL'}")
    print(f"  Chiral anticommute: {'✓ PASS' if results_3d['chiral_anticommute'] else '✗ FAIL'}")
    print(f"  (γ^5)² = I: {'✓ PASS' if results_3d['gamma5_squared'] else '✗ FAIL'}")
    print(f"  Max error: {results_3d['max_error']:.2e}")

    # Test 1+1D
    print("\n[2/2] Testing 1+1D Gamma Matrices...")
    results_1d = verify_gamma_algebra_1plus1()
    print(f"  Anticommutator: {'✓ PASS' if results_1d['anticommutator_correct'] else '✗ FAIL'}")
    print(f"  Chiral properties: {'✓ PASS' if results_1d['chiral_properties'] else '✗ FAIL'}")
    print(f"  Max error: {results_1d['max_error']:.2e}")

    # Test MSFT mass operator
    print("\n[3/3] Testing MSFT Mass Operator...")
    γ0, γ1, γ2, γ3, γ5 = get_gamma_matrices_3plus1()

    # Pure scalar (θ=0)
    m_s = mass_operator_MSFT(2.0, 0.5, 0.0, γ5)
    scalar_correct = np.allclose(m_s, 1.0 * np.eye(4))
    print(f"  Scalar mass (θ=0): {'✓ PASS' if scalar_correct else '✗ FAIL'}")

    # Pure pseudoscalar (θ=π/2)
    m_p = mass_operator_MSFT(2.0, 0.5, np.pi/2, γ5)
    pseudo_correct = np.allclose(m_p, 1.0j * γ5)
    print(f"  Pseudoscalar mass (θ=π/2): {'✓ PASS' if pseudo_correct else '✗ FAIL'}")

    # Mixed (θ=π/4)
    m_mix = mass_operator_MSFT(2.0, 0.5, np.pi/4, γ5)
    expected_mix = 1.0 * (np.cos(np.pi/4) * np.eye(4) + 1j * np.sin(np.pi/4) * γ5)
    mixed_correct = np.allclose(m_mix, expected_mix)
    print(f"  Mixed mass (θ=π/4): {'✓ PASS' if mixed_correct else '✗ FAIL'}")

    print("\n" + "=" * 70)
    all_pass = (results_3d['anticommutator_correct'] and
                results_3d['hermiticity_correct'] and
                results_3d['chiral_anticommute'] and
                results_3d['gamma5_squared'] and
                results_1d['anticommutator_correct'] and
                results_1d['chiral_properties'] and
                scalar_correct and pseudo_correct and mixed_correct)

    if all_pass:
        print(" ✓✓✓ ALL TESTS PASSED ✓✓✓")
    else:
        print(" ✗✗✗ SOME TESTS FAILED ✗✗✗")
    print("=" * 70)
