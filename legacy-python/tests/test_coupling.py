"""
Comprehensive test suite for Kuramoto coupling functions.

Tests all coupling classes in kuramoto.core.coupling module:
- SinusoidalCoupling (standard all-to-all)
- NetworkCoupling (graph-based)
- PulseCoupling (event-driven)
- CustomCoupling (user-defined)
- HigherHarmonicCoupling (generalized)
"""

import sys
import os
import numpy as np
from numpy.testing import assert_allclose
import unittest

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from kuramoto.core.coupling import (
    Coupling, SinusoidalCoupling, NetworkCoupling,
    PulseCoupling, CustomCoupling, HigherHarmonicCoupling
)


class TestSinusoidalCoupling(unittest.TestCase):
    """Test standard all-to-all sinusoidal coupling."""

    def test_initialization(self):
        """Test coupling initialization with different strengths."""
        coupling = SinusoidalCoupling(strength=2.5)
        self.assertEqual(coupling.strength, 2.5)
        self.assertEqual(coupling.K, 2.5)  # K is alias for strength

    def test_zero_coupling(self):
        """Test that zero coupling strength gives zero field."""
        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)
        coupling = SinusoidalCoupling(strength=0.0)
        field = coupling.compute_field(phases)
        assert_allclose(field, np.zeros(N), atol=1e-10)

    def test_symmetric_configuration(self):
        """Test coupling field for symmetric phase configurations."""
        # Two oscillators in anti-phase should have opposite fields
        phases = np.array([0, np.pi])
        coupling = SinusoidalCoupling(strength=2.0)
        field = coupling.compute_field(phases)

        # Fields should be equal and opposite
        assert_allclose(field[0], -field[1], rtol=1e-10)

    def test_fully_synchronized(self):
        """Test that synchronized oscillators have zero coupling field."""
        N = 20
        phases = np.ones(N) * np.pi/3  # All at same phase
        coupling = SinusoidalCoupling(strength=3.0)
        field = coupling.compute_field(phases)

        # No phase differences -> zero field
        assert_allclose(field, np.zeros(N), atol=1e-10)

    def test_mean_field_property(self):
        """Test mean-field nature: sum of fields is zero."""
        N = 50
        phases = np.random.uniform(0, 2*np.pi, N)
        coupling = SinusoidalCoupling(strength=1.5)
        field = coupling.compute_field(phases)

        # Conservation: total field sums to zero
        assert_allclose(np.sum(field), 0.0, atol=1e-10)

    def test_coupling_strength_scaling(self):
        """Test that field scales linearly with coupling strength."""
        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)

        K1, K2 = 1.0, 3.0
        coupling1 = SinusoidalCoupling(strength=K1)
        coupling2 = SinusoidalCoupling(strength=K2)

        field1 = coupling1.compute_field(phases)
        field2 = coupling2.compute_field(phases)

        # Field should scale with K
        assert_allclose(field2, field1 * (K2/K1), rtol=1e-10)

    def test_large_system(self):
        """Test coupling computation for large N."""
        N = 1000
        phases = np.random.uniform(0, 2*np.pi, N)
        coupling = SinusoidalCoupling(strength=2.0)

        field = coupling.compute_field(phases)

        # Should return field for all oscillators
        self.assertEqual(len(field), N)
        # Field should be bounded by K
        self.assertTrue(np.all(np.abs(field) <= coupling.strength))

    def test_specific_two_oscillator_case(self):
        """Test analytical result for two oscillators."""
        # For two oscillators: H_1 = (K/2)sin(θ_2 - θ_1)
        theta1, theta2 = np.pi/4, np.pi/3
        phases = np.array([theta1, theta2])
        K = 4.0

        coupling = SinusoidalCoupling(strength=K)
        field = coupling.compute_field(phases)

        expected_field1 = (K/2) * np.sin(theta2 - theta1)
        expected_field2 = (K/2) * np.sin(theta1 - theta2)

        assert_allclose(field[0], expected_field1, rtol=1e-10)
        assert_allclose(field[1], expected_field2, rtol=1e-10)


class TestNetworkCoupling(unittest.TestCase):
    """Test network-based coupling on graphs."""

    def test_initialization(self):
        """Test network coupling initialization."""
        adj = np.array([[0, 1, 1],
                        [1, 0, 1],
                        [1, 1, 0]])
        coupling = NetworkCoupling(adj, strength=2.0)

        self.assertEqual(coupling.strength, 2.0)
        self.assertEqual(coupling.N, 3)
        assert_allclose(coupling.adjacency, adj)

    def test_isolated_node(self):
        """Test that isolated nodes have zero field."""
        # Node 2 is isolated
        adj = np.array([[0, 1, 0],
                        [1, 0, 0],
                        [0, 0, 0]])
        phases = np.random.uniform(0, 2*np.pi, 3)

        coupling = NetworkCoupling(adj, strength=1.0, normalize=True)
        field = coupling.compute_field(phases)

        # Isolated node should have zero field
        assert_allclose(field[2], 0.0, atol=1e-10)

    def test_chain_network(self):
        """Test coupling on a chain network."""
        # Chain: 0-1-2-3
        N = 4
        adj = np.zeros((N, N))
        for i in range(N-1):
            adj[i, i+1] = 1
            adj[i+1, i] = 1

        # Use non-symmetric phase configuration to ensure non-zero field
        phases = np.array([0, np.pi/3, 2*np.pi/3, np.pi])
        coupling = NetworkCoupling(adj, strength=2.0, normalize=True)
        field = coupling.compute_field(phases)

        # End nodes have one neighbor, middle nodes have two
        # Field magnitude should reflect connectivity
        self.assertNotAlmostEqual(field[0], 0.0)  # Has neighbor
        self.assertNotAlmostEqual(field[1], 0.0)  # Has 2 neighbors

    def test_weighted_network(self):
        """Test coupling with weighted adjacency matrix."""
        # Weighted triangle
        adj = np.array([[0, 2, 1],
                        [2, 0, 3],
                        [1, 3, 0]])
        phases = np.array([0, np.pi/2, np.pi])

        coupling = NetworkCoupling(adj, strength=1.0, normalize=False)
        field = coupling.compute_field(phases)

        # Field should incorporate weights
        # Stronger connections should have larger influence
        self.assertEqual(len(field), 3)

    def test_normalization_effect(self):
        """Test effect of degree normalization."""
        # Star network: node 0 connected to all others
        N = 5
        adj = np.zeros((N, N))
        adj[0, 1:] = 1
        adj[1:, 0] = 1

        phases = np.random.uniform(0, 2*np.pi, N)

        # Without normalization
        coupling_unnorm = NetworkCoupling(adj, strength=1.0, normalize=False)
        field_unnorm = coupling_unnorm.compute_field(phases)

        # With normalization
        coupling_norm = NetworkCoupling(adj, strength=1.0, normalize=True)
        field_norm = coupling_norm.compute_field(phases)

        # Central node field should be reduced by normalization
        # (divided by degree 4)
        self.assertLess(np.abs(field_norm[0]), np.abs(field_unnorm[0]))

    def test_complete_graph_matches_all_to_all(self):
        """Test that complete graph gives same result as all-to-all."""
        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)
        K = 2.5

        # Complete graph (all connections)
        adj = np.ones((N, N)) - np.eye(N)
        # Network coupling without normalization sums over (N-1) neighbors
        # Sinusoidal coupling divides by N
        # To match, we need network strength K*(N-1)/N
        network_coupling = NetworkCoupling(adj, strength=K/N, normalize=False)
        network_field = network_coupling.compute_field(phases)

        # All-to-all sinusoidal
        sinusoidal_coupling = SinusoidalCoupling(strength=K)
        sinusoidal_field = sinusoidal_coupling.compute_field(phases)

        # Now they should match
        assert_allclose(network_field, sinusoidal_field, rtol=1e-10)


class TestPulseCoupling(unittest.TestCase):
    """Test pulse-coupled oscillators."""

    def test_initialization(self):
        """Test pulse coupling initialization."""
        coupling = PulseCoupling(strength=0.1)
        self.assertEqual(coupling.strength, 0.1)
        self.assertIsNotNone(coupling.pulse_shape)

    def test_no_firing(self):
        """Test that no pulses occur when no oscillators are firing."""
        N = 10
        phases = np.random.uniform(0, np.pi, N)  # All far from 2π
        coupling = PulseCoupling(strength=0.2)
        field = coupling.compute_field(phases)

        # No firing -> zero field
        assert_allclose(field, np.zeros(N), atol=1e-10)

    def test_single_firing(self):
        """Test pulse when one oscillator fires."""
        N = 5
        phases = np.array([np.pi/2, np.pi, 3*np.pi/2, 1.95*np.pi, 0.1])
        phases[3] = 2*np.pi - 0.05  # Near firing

        coupling = PulseCoupling(strength=0.3)
        field = coupling.compute_field(phases)

        # Firing oscillator should have zero field
        self.assertAlmostEqual(field[3], 0.0)
        # Others should receive pulse influence
        for i in [0, 1, 2, 4]:
            self.assertGreater(field[i], 0)

    def test_multiple_firing(self):
        """Test pulses when multiple oscillators fire."""
        N = 6
        phases = np.ones(N) * np.pi
        phases[1] = 2*np.pi - 0.05  # Firing
        phases[3] = 2*np.pi - 0.02  # Firing

        coupling = PulseCoupling(strength=0.25)
        field = coupling.compute_field(phases)

        # Non-firing oscillators receive pulses from two sources
        n_firing = 2
        expected_pulse = coupling.strength * n_firing / N

        for i in [0, 2, 4, 5]:
            assert_allclose(field[i], expected_pulse, rtol=0.1)

    def test_custom_pulse_shape(self):
        """Test custom pulse response function."""
        def custom_pulse(phase):
            return min(phase + 0.5, 2*np.pi)

        coupling = PulseCoupling(strength=0.2, pulse_shape=custom_pulse)
        self.assertEqual(coupling.pulse_shape, custom_pulse)

        # Test it's used
        test_phase = np.pi
        result = coupling.pulse_shape(test_phase)
        self.assertEqual(result, test_phase + 0.5)


class TestCustomCoupling(unittest.TestCase):
    """Test user-defined custom coupling."""

    def test_constant_coupling(self):
        """Test custom coupling with constant field."""
        def constant_field(phases):
            return np.ones_like(phases) * 0.5

        coupling = CustomCoupling(constant_field)
        phases = np.random.uniform(0, 2*np.pi, 10)
        field = coupling.compute_field(phases)

        assert_allclose(field, 0.5 * np.ones(10), rtol=1e-10)

    def test_linear_coupling(self):
        """Test custom linear coupling."""
        def linear_field(phases):
            # Field proportional to phase
            return phases / (2*np.pi)

        coupling = CustomCoupling(linear_field)
        phases = np.linspace(0, 2*np.pi, 5)
        field = coupling.compute_field(phases)

        expected = phases / (2*np.pi)
        assert_allclose(field, expected, rtol=1e-10)

    def test_complex_custom_function(self):
        """Test complex custom coupling function."""
        def complex_coupling(phases):
            N = len(phases)
            # Combine multiple effects
            mean_phase = np.mean(phases)
            dispersion = np.std(phases)
            field = np.sin(phases - mean_phase) * dispersion
            return field

        coupling = CustomCoupling(complex_coupling)
        phases = np.array([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        field = coupling.compute_field(phases)

        # Should return array of same length
        self.assertEqual(len(field), len(phases))
        # Check it computed something non-trivial
        self.assertFalse(np.allclose(field, 0))


class TestHigherHarmonicCoupling(unittest.TestCase):
    """Test coupling with higher harmonics."""

    def test_single_harmonic(self):
        """Test single harmonic (standard Kuramoto)."""
        harmonics = {1: (2.0, 0)}  # n=1, K=2, α=0
        coupling = HigherHarmonicCoupling(harmonics)

        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)
        field = coupling.compute_field(phases)

        # Compare with standard sinusoidal coupling
        standard = SinusoidalCoupling(strength=2.0)
        standard_field = standard.compute_field(phases)

        assert_allclose(field, standard_field, rtol=1e-10)

    def test_second_harmonic(self):
        """Test second harmonic coupling."""
        harmonics = {2: (1.5, 0)}  # n=2, K=1.5, α=0
        coupling = HigherHarmonicCoupling(harmonics)

        # Test with simple two-oscillator system
        phases = np.array([0, np.pi/2])
        field = coupling.compute_field(phases)

        # Manually compute expected field
        # H_1 = (K/N) sin(2*(θ_2 - θ_1))
        expected_1 = (1.5/2) * np.sin(2*(phases[1] - phases[0]))
        expected_2 = (1.5/2) * np.sin(2*(phases[0] - phases[1]))

        assert_allclose(field[0], expected_1, rtol=1e-10)
        assert_allclose(field[1], expected_2, rtol=1e-10)

    def test_multiple_harmonics(self):
        """Test coupling with multiple harmonics."""
        harmonics = {
            1: (2.0, 0),        # First harmonic
            2: (0.5, np.pi/4),  # Second with phase shift
            3: (0.2, 0)         # Third harmonic
        }
        coupling = HigherHarmonicCoupling(harmonics)

        N = 20
        phases = np.random.uniform(0, 2*np.pi, N)
        field = coupling.compute_field(phases)

        # Field should be superposition of all harmonics
        self.assertEqual(len(field), N)
        # Check that field is bounded (rough check)
        max_possible = sum(abs(K) for n, (K, alpha) in harmonics.items())
        self.assertTrue(np.all(np.abs(field) <= max_possible * 1.5))

    def test_phase_shift_effect(self):
        """Test effect of phase shift α in harmonics."""
        # Same coupling strength, different phase shifts
        harmonics1 = {1: (1.0, 0)}
        harmonics2 = {1: (1.0, np.pi/2)}

        coupling1 = HigherHarmonicCoupling(harmonics1)
        coupling2 = HigherHarmonicCoupling(harmonics2)

        phases = np.array([0, np.pi/3, 2*np.pi/3])
        field1 = coupling1.compute_field(phases)
        field2 = coupling2.compute_field(phases)

        # Fields should be different due to phase shift
        self.assertFalse(np.allclose(field1, field2))

    def test_zero_harmonic_coupling(self):
        """Test that zero coupling gives zero field."""
        harmonics = {1: (0.0, 0), 2: (0.0, np.pi/4)}
        coupling = HigherHarmonicCoupling(harmonics)

        phases = np.random.uniform(0, 2*np.pi, 15)
        field = coupling.compute_field(phases)

        assert_allclose(field, np.zeros(15), atol=1e-10)


class TestCouplingIntegration(unittest.TestCase):
    """Integration tests for coupling classes."""

    def test_all_couplings_implement_interface(self):
        """Test that all coupling classes implement the base interface."""
        # Create instances of each coupling type
        couplings = [
            SinusoidalCoupling(strength=1.0),
            NetworkCoupling(np.ones((3,3)), strength=1.0),
            PulseCoupling(strength=0.1),
            CustomCoupling(lambda p: np.zeros_like(p)),
            HigherHarmonicCoupling({1: (1.0, 0)})
        ]

        phases = np.array([0, np.pi/2, np.pi])

        for coupling in couplings:
            # Should be instance of base class
            self.assertIsInstance(coupling, Coupling)
            # Should have compute_field method
            self.assertTrue(hasattr(coupling, 'compute_field'))
            # Method should work
            field = coupling.compute_field(phases)
            self.assertEqual(len(field), len(phases))

    def test_coupling_composition(self):
        """Test combining multiple coupling types."""
        N = 10
        phases = np.random.uniform(0, 2*np.pi, N)

        # Create different couplings
        sinusoidal = SinusoidalCoupling(strength=1.0)
        harmonic = HigherHarmonicCoupling({2: (0.5, 0)})

        # Compute combined field
        field1 = sinusoidal.compute_field(phases)
        field2 = harmonic.compute_field(phases)
        combined_field = field1 + field2

        # Combined field should be well-defined
        self.assertEqual(len(combined_field), N)
        self.assertFalse(np.any(np.isnan(combined_field)))


if __name__ == '__main__':
    unittest.main(verbosity=2)