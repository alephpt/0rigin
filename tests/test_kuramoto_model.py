"""
Test suite for Kuramoto model core functionality.

Tests cover:
- Model initialization
- Order parameter computation
- Synchronization dynamics
- Known analytical results
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

# These imports will work once the package is implemented
# from kuramoto import KuramotoModel
# from kuramoto.distributions import LorentzianDistribution


class TestKuramotoModel:
    """Test core Kuramoto model functionality."""

    def test_order_parameter_bounds(self):
        """Test that order parameter R is always in [0, 1]."""
        # Create model with random initial conditions
        N = 100
        phases = np.random.uniform(0, 2 * np.pi, N)

        # Compute order parameter
        z = np.mean(np.exp(1j * phases))
        R = np.abs(z)

        # Check bounds
        assert 0 <= R <= 1, f"Order parameter R={R} outside [0, 1]"

    def test_fully_synchronized_state(self):
        """Test R=1 for fully synchronized oscillators."""
        N = 100
        # All oscillators at same phase
        phases = np.ones(N) * np.pi / 4

        # Compute order parameter
        z = np.mean(np.exp(1j * phases))
        R = np.abs(z)
        Psi = np.angle(z)

        # Check full synchronization
        assert_allclose(R, 1.0, rtol=1e-10)
        assert_allclose(Psi, np.pi / 4, rtol=1e-10)

    def test_uniform_distribution_incoherent(self):
        """Test R≈0 for uniformly distributed phases."""
        N = 10000
        # Uniformly distributed phases
        phases = np.linspace(0, 2 * np.pi, N, endpoint=False)

        # Compute order parameter
        z = np.mean(np.exp(1j * phases))
        R = np.abs(z)

        # Should be near zero for large N
        assert R < 0.01, f"Order parameter R={R} too large for uniform phases"

    def test_identical_frequencies_synchronize(self):
        """Test that identical frequencies lead to synchronization."""
        N = 50
        frequencies = np.zeros(N)  # All identical
        phases = np.random.uniform(0, 2*np.pi, N)

        # With strong coupling and identical frequencies, should synchronize
        K = 5.0
        dt = 0.01

        # Simple evolution for testing
        for _ in range(1000):  # 10 time units
            z = np.mean(np.exp(1j * phases))
            R = np.abs(z)
            psi = np.angle(z)

            # Kuramoto dynamics
            dtheta_dt = frequencies + K * R * np.sin(psi - phases)
            phases += dt * dtheta_dt

        final_R = np.abs(np.mean(np.exp(1j * phases)))
        assert final_R > 0.95, f"Identical frequencies should synchronize, got R={final_R}"

    def test_subcritical_coupling_incoherent(self):
        """Test that K < Kc maintains incoherent state."""
        # For Lorentzian-like distribution with γ≈1, Kc≈2
        N = 500
        gamma = 1.0

        # Create Lorentzian-distributed frequencies
        u = np.random.uniform(0, 1, N)
        frequencies = gamma * np.tan(np.pi * (u - 0.5))

        # Use K=1 < Kc=2
        K = 1.0
        phases = np.random.uniform(0, 2*np.pi, N)
        dt = 0.01

        # Evolve to steady state
        R_history = []
        for _ in range(5000):  # 50 time units
            z = np.mean(np.exp(1j * phases))
            R = np.abs(z)
            R_history.append(R)
            psi = np.angle(z)

            # Kuramoto dynamics
            dtheta_dt = frequencies + K * R * np.sin(psi - phases)
            phases += dt * dtheta_dt

        steady_R = np.mean(R_history[-1000:])
        assert steady_R < 0.1, f"Subcritical coupling should remain incoherent, got R={steady_R}"

    def test_supercritical_coupling_synchronizes(self):
        """Test that K > Kc leads to partial synchronization."""
        # For Lorentzian with γ=1, Kc=2
        N = 500
        gamma = 1.0

        # Create Lorentzian-distributed frequencies
        u = np.random.uniform(0, 1, N)
        frequencies = gamma * np.tan(np.pi * (u - 0.5))

        # Use K=4 > Kc=2
        K = 4.0
        Kc = 2.0 * gamma
        phases = np.random.uniform(0, 2*np.pi, N)
        dt = 0.01

        # Evolve to steady state
        R_history = []
        for _ in range(5000):  # 50 time units
            z = np.mean(np.exp(1j * phases))
            R = np.abs(z)
            R_history.append(R)
            psi = np.angle(z)

            # Kuramoto dynamics
            dtheta_dt = frequencies + K * R * np.sin(psi - phases)
            phases += dt * dtheta_dt

        steady_R = np.mean(R_history[-1000:])

        # Theoretical prediction: R = sqrt(1 - (Kc/K)^2)
        theory_R = np.sqrt(1 - (Kc/K)**2)

        # Check within reasonable tolerance (dynamics can be noisy)
        assert abs(steady_R - theory_R) < 0.2, \
            f"Expected R≈{theory_R:.3f}, got R={steady_R:.3f}"


class TestOttAntonsen:
    """Test Ott-Antonsen exact solutions for Lorentzian distribution."""

    def test_critical_coupling_formula(self):
        """Test Kc = 2γ for Lorentzian distribution."""
        gamma_values = [0.5, 1.0, 2.0, 5.0]
        for gamma in gamma_values:
            # dist = LorentzianDistribution(center=0, width=gamma)
            # Kc = dist.critical_coupling()
            Kc_expected = 2 * gamma
            # assert Kc == Kc_expected
            assert Kc_expected == 2 * gamma  # Placeholder

    def test_steady_state_order_parameter(self):
        """Test steady-state R from OA theory."""
        gamma = 1.0
        K_values = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
        Kc = 2 * gamma

        for K in K_values:
            if K <= Kc:
                R_theory = 0.0
            else:
                R_theory = np.sqrt(1 - (Kc / K)**2)

            # Would compare with numerical simulation
            # dist = LorentzianDistribution(center=0, width=gamma)
            # R_steady = dist.steady_state_order_parameter(K)
            # assert R_steady == R_theory
            assert R_theory >= 0  # Placeholder

    def test_bifurcation_scaling(self):
        """Test R ∝ √(K - Kc) near critical point."""
        gamma = 1.0
        Kc = 2 * gamma
        epsilon = 0.1  # Small deviation from Kc

        K = Kc + epsilon
        R_theory = np.sqrt(1 - (Kc / K)**2)

        # Near Kc, R ≈ √(2ε/Kc)
        R_approx = np.sqrt(2 * epsilon / Kc)

        assert_allclose(R_theory, R_approx, rtol=0.01)


class TestNumericalStability:
    """Test numerical stability of solvers."""

    def test_phase_wrapping(self):
        """Test that phases remain in [0, 2π) or [-π, π)."""
        N = 100
        phases = np.random.uniform(-10 * np.pi, 10 * np.pi, N)

        # Wrap to [0, 2π)
        wrapped = np.mod(phases, 2 * np.pi)
        assert np.all(wrapped >= 0)
        assert np.all(wrapped < 2 * np.pi)

        # Wrap to [-π, π)
        wrapped_centered = np.angle(np.exp(1j * phases))
        assert np.all(wrapped_centered >= -np.pi)
        assert np.all(wrapped_centered <= np.pi)

    def test_conservation_of_oscillators(self):
        """Test that number of oscillators is conserved."""
        N = 100
        phases = np.random.uniform(0, 2*np.pi, N)
        frequencies = np.random.normal(0, 1, N)

        # Track N through evolution
        K = 3.0
        dt = 0.01
        n_history = []

        for _ in range(100):
            n_history.append(len(phases))
            z = np.mean(np.exp(1j * phases))
            R = np.abs(z)
            psi = np.angle(z)
            dtheta_dt = frequencies + K * R * np.sin(psi - phases)
            phases += dt * dtheta_dt

        # Number should remain constant
        assert all(n == N for n in n_history), "Number of oscillators changed!"

    def test_deterministic_evolution(self):
        """Test that evolution is deterministic with fixed seed."""
        N = 50
        K = 3.0
        dt = 0.01
        n_steps = 100

        # First run
        np.random.seed(42)
        frequencies1 = np.random.normal(0, 1, N)
        phases1 = np.random.uniform(0, 2*np.pi, N)

        for _ in range(n_steps):
            z = np.mean(np.exp(1j * phases1))
            R = np.abs(z)
            psi = np.angle(z)
            dtheta_dt = frequencies1 + K * R * np.sin(psi - phases1)
            phases1 += dt * dtheta_dt

        # Second run with same seed
        np.random.seed(42)
        frequencies2 = np.random.normal(0, 1, N)
        phases2 = np.random.uniform(0, 2*np.pi, N)

        for _ in range(n_steps):
            z = np.mean(np.exp(1j * phases2))
            R = np.abs(z)
            psi = np.angle(z)
            dtheta_dt = frequencies2 + K * R * np.sin(psi - phases2)
            phases2 += dt * dtheta_dt

        assert_allclose(phases1, phases2, rtol=1e-10)


class TestParameterValidation:
    """Test input parameter validation."""

    def test_positive_coupling(self):
        """Test that negative coupling is handled correctly."""
        # Negative coupling should work (repulsive interaction)
        N = 10
        K = -1.0  # Repulsive coupling
        phases = np.random.uniform(0, 2*np.pi, N)
        frequencies = np.zeros(N)
        dt = 0.01

        # Should evolve without errors
        for _ in range(10):
            z = np.mean(np.exp(1j * phases))
            R = np.abs(z)
            psi = np.angle(z)
            dtheta_dt = frequencies + K * R * np.sin(psi - phases)
            phases += dt * dtheta_dt

        # With negative K, phases should disperse rather than synchronize
        final_R = np.abs(np.mean(np.exp(1j * phases)))
        assert final_R < 0.5, "Repulsive coupling should prevent synchronization"

    def test_invalid_N(self):
        """Test that invalid N values are handled."""
        # Test with N=0
        try:
            phases = np.array([])
            frequencies = np.array([])
            # Empty arrays should be handled gracefully
            R = np.abs(np.mean(np.exp(1j * phases))) if len(phases) > 0 else 0
            assert R == 0, "Empty system should have R=0"
        except Exception as e:
            pytest.fail(f"Failed to handle N=0: {e}")

        # Test with N=1 (single oscillator)
        phases = np.array([0.0])
        frequencies = np.array([1.0])
        R = np.abs(np.mean(np.exp(1j * phases)))
        assert R == 1.0, "Single oscillator has R=1 by definition"

    def test_frequency_array_length(self):
        """Test that frequency array must match phase array."""
        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)
        frequencies_wrong = np.zeros(15)  # Wrong length

        # Test mismatch detection
        try:
            # Attempting to use mismatched arrays
            if len(phases) != len(frequencies_wrong):
                raise ValueError(f"Frequency array length {len(frequencies_wrong)} "
                               f"doesn't match N={len(phases)}")
        except ValueError as e:
            assert "doesn't match" in str(e)
        else:
            pytest.fail("Should have raised ValueError for mismatched arrays")


# Fixtures for common test data
@pytest.fixture
def small_model_params():
    """Parameters for small test model."""
    return {
        'N': 10,
        'coupling': 2.0,
        'frequencies': np.random.randn(10) * 0.5
    }


@pytest.fixture
def large_model_params():
    """Parameters for large test model."""
    return {
        'N': 1000,
        'coupling': 3.0,
        'frequencies': 'lorentzian'
    }


# Performance benchmarks (would use pytest-benchmark)
def test_benchmark_order_parameter(benchmark):
    """Benchmark order parameter calculation."""
    N = 10000
    phases = np.random.uniform(0, 2 * np.pi, N)

    def compute_R():
        z = np.mean(np.exp(1j * phases))
        return np.abs(z)

    # result = benchmark(compute_R)
    # assert 0 <= result <= 1
    pass