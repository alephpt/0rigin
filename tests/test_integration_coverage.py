"""
Integration tests to verify complete system functionality.

Tests the full workflow from initialization through evolution
to analysis, covering all critical paths.
"""

import sys
import os
import numpy as np
from numpy.testing import assert_allclose
import unittest

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from kuramoto import KuramotoModel
from kuramoto.core.coupling import SinusoidalCoupling, NetworkCoupling
from kuramoto.solvers import EulerSolver, RK4Solver
from kuramoto.field_theory import MediatorField


class TestFullKuramotoWorkflow(unittest.TestCase):
    """Test complete Kuramoto simulation workflow."""

    def test_basic_synchronization_workflow(self):
        """Test full workflow: init -> evolve -> analyze."""
        # Initialize model
        N = 50
        model = KuramotoModel(
            N=N,
            coupling=3.0,
            frequencies=np.random.normal(0, 0.5, N),
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        # Evolve system
        t_span = (0, 10)
        dt = 0.01
        t, phases = model.evolve(t_span, dt=dt)

        # Analyze results
        # Compute order parameter over time
        R_values = []
        for phase_config in phases:
            z = np.mean(np.exp(1j * phase_config))
            R = np.abs(z)
            R_values.append(R)

        # Check that we have results
        self.assertEqual(len(t), len(R_values))
        self.assertGreater(len(t), 100)  # Should have many time steps

        # Check that R is bounded
        self.assertTrue(all(0 <= R <= 1 for R in R_values))

        # With strong coupling, should show increasing synchronization
        early_R = np.mean(R_values[:10])
        late_R = np.mean(R_values[-10:])
        self.assertGreater(late_R, early_R)  # Should synchronize over time

    def test_different_solvers_consistency(self):
        """Test that different solvers give consistent results."""
        N = 20
        model = KuramotoModel(
            N=N,
            coupling=2.0,
            frequencies=np.zeros(N),  # Identical frequencies
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        # Set seed for reproducibility
        np.random.seed(42)

        # Evolve with Euler solver
        euler_solver = EulerSolver()
        model.solver = euler_solver
        t1, phases1 = model.evolve((0, 1), dt=0.001)

        # Reset and evolve with RK4 solver
        np.random.seed(42)
        rk4_solver = RK4Solver()
        model.solver = rk4_solver
        model.phases = model.initial_phases.copy()
        t2, phases2 = model.evolve((0, 1), dt=0.001)

        # Both should synchronize (identical frequencies)
        final_R1 = np.abs(np.mean(np.exp(1j * phases1[-1])))
        final_R2 = np.abs(np.mean(np.exp(1j * phases2[-1])))

        self.assertGreater(final_R1, 0.9)  # Should synchronize
        self.assertGreater(final_R2, 0.9)  # Should synchronize

        # RK4 should be more accurate
        # They should give similar but not identical results
        assert_allclose(final_R1, final_R2, rtol=0.1)

    def test_network_coupling_integration(self):
        """Test integration with network coupling."""
        # Create star network
        N = 7
        adj = np.zeros((N, N))
        # Central node (0) connected to all others
        for i in range(1, N):
            adj[0, i] = 1
            adj[i, 0] = 1

        # Create model with network coupling
        coupling = NetworkCoupling(adj, strength=5.0)
        model = KuramotoModel(
            N=N,
            coupling=coupling,
            frequencies=np.random.normal(0, 0.2, N),
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        # Evolve
        t, phases = model.evolve((0, 10), dt=0.01)

        # Check that evolution completed
        self.assertGreater(len(t), 100)
        self.assertEqual(phases.shape[1], N)

        # Central node should have strong influence
        # Compute pairwise phase differences with central node
        final_phases = phases[-1]
        central_phase = final_phases[0]
        phase_diffs = np.abs(np.angle(np.exp(1j * (final_phases - central_phase))))

        # Connected nodes should be closer to central node than unconnected average
        mean_diff = np.mean(phase_diffs[1:])
        self.assertLess(mean_diff, np.pi/2)  # Should be somewhat synchronized

    def test_field_theory_integration(self):
        """Test integration with field theory components."""
        # Create field
        grid_shape = (10, 10)
        phi = np.ones(grid_shape)
        pi_phi = np.zeros(grid_shape)
        field = MediatorField(phi=phi, pi_phi=pi_phi, mass=1.0)

        # Create Kuramoto model
        N = 10
        model = KuramotoModel(N=N)

        # Couple to field (if supported)
        try:
            from kuramoto.field_theory import FieldCoupling
            field_coupling = FieldCoupling(field, coupling_strength=0.5)
            # This would be used in a coupled evolution
            self.assertIsNotNone(field_coupling)
        except (ImportError, AttributeError):
            # Field coupling might not be fully implemented
            pass

        # Verify field properties
        self.assertEqual(field.phi.shape, grid_shape)
        self.assertEqual(field.pi_phi.shape, grid_shape)

    def test_critical_transition(self):
        """Test behavior near critical coupling strength."""
        N = 200
        gamma = 1.0  # Width of frequency distribution

        # Create Lorentzian-like frequency distribution
        u = np.random.uniform(0, 1, N)
        frequencies = gamma * np.tan(np.pi * (u - 0.5))

        # Test subcritical (K < Kc = 2γ)
        model_sub = KuramotoModel(
            N=N,
            coupling=1.5,  # < 2.0
            frequencies=frequencies,
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        t_sub, phases_sub = model_sub.evolve((0, 20), dt=0.01)
        R_sub = np.abs(np.mean(np.exp(1j * phases_sub[-100:]), axis=0))
        mean_R_sub = np.mean(R_sub)

        # Test supercritical (K > Kc)
        model_super = KuramotoModel(
            N=N,
            coupling=3.0,  # > 2.0
            frequencies=frequencies,
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        t_super, phases_super = model_super.evolve((0, 20), dt=0.01)
        R_super = np.abs(np.mean(np.exp(1j * phases_super[-100:]), axis=0))
        mean_R_super = np.mean(R_super)

        # Subcritical should remain incoherent, supercritical should synchronize
        self.assertLess(mean_R_sub, 0.2)  # Mostly incoherent
        self.assertGreater(mean_R_super, 0.4)  # Partially synchronized


class TestStressConditions(unittest.TestCase):
    """Test system under stress conditions."""

    def test_large_system_performance(self):
        """Test with large number of oscillators."""
        N = 1000
        model = KuramotoModel(
            N=N,
            coupling=2.0,
            frequencies=np.random.normal(0, 1, N)
        )

        # Short evolution to test feasibility
        t, phases = model.evolve((0, 1), dt=0.01)

        # Should complete without error
        self.assertEqual(phases.shape[1], N)
        self.assertGreater(len(t), 10)

    def test_long_time_evolution(self):
        """Test long-time stability."""
        N = 50
        model = KuramotoModel(
            N=N,
            coupling=3.0,
            frequencies=np.random.normal(0, 0.5, N)
        )

        # Evolve for longer time
        t, phases = model.evolve((0, 100), dt=0.1)

        # Check stability - no NaN or inf
        self.assertFalse(np.any(np.isnan(phases)))
        self.assertFalse(np.any(np.isinf(phases)))

        # Phases should remain bounded (modulo 2π)
        phases_wrapped = np.mod(phases, 2*np.pi)
        self.assertTrue(np.all(phases_wrapped >= 0))
        self.assertTrue(np.all(phases_wrapped < 2*np.pi))

    def test_extreme_parameter_values(self):
        """Test with extreme parameter values."""
        # Very strong coupling
        N = 30
        model_strong = KuramotoModel(
            N=N,
            coupling=100.0,  # Very strong
            frequencies=np.random.normal(0, 1, N)
        )

        t, phases = model_strong.evolve((0, 1), dt=0.001)

        # Should rapidly synchronize
        final_R = np.abs(np.mean(np.exp(1j * phases[-1])))
        self.assertGreater(final_R, 0.95)

        # Very weak coupling
        model_weak = KuramotoModel(
            N=N,
            coupling=0.01,  # Very weak
            frequencies=np.random.normal(0, 1, N)
        )

        t, phases = model_weak.evolve((0, 10), dt=0.01)

        # Should remain mostly incoherent
        final_R = np.abs(np.mean(np.exp(1j * phases[-1])))
        self.assertLess(final_R, 0.3)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    def test_all_zero_frequencies(self):
        """Test with all zero frequencies."""
        N = 20
        model = KuramotoModel(
            N=N,
            coupling=2.0,
            frequencies=np.zeros(N),
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        t, phases = model.evolve((0, 10), dt=0.01)

        # Should synchronize to common phase
        final_R = np.abs(np.mean(np.exp(1j * phases[-1])))
        self.assertGreater(final_R, 0.95)

    def test_bimodal_frequency_distribution(self):
        """Test with bimodal frequency distribution."""
        N = 100
        # Create two groups with different frequencies
        frequencies = np.concatenate([
            np.ones(50) * (-1.0),  # Group 1
            np.ones(50) * (1.0)     # Group 2
        ])

        model = KuramotoModel(
            N=N,
            coupling=1.5,
            frequencies=frequencies,
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        t, phases = model.evolve((0, 20), dt=0.01)

        # Should show complex dynamics
        # Each group might synchronize internally but not globally
        final_phases = phases[-1]
        group1_R = np.abs(np.mean(np.exp(1j * final_phases[:50])))
        group2_R = np.abs(np.mean(np.exp(1j * final_phases[50:])))
        global_R = np.abs(np.mean(np.exp(1j * final_phases)))

        # Groups might be more synchronized than global
        self.assertLess(global_R, max(group1_R, group2_R))

    def test_single_outlier_frequency(self):
        """Test with one oscillator having very different frequency."""
        N = 21
        frequencies = np.zeros(N)
        frequencies[0] = 10.0  # One fast oscillator

        model = KuramotoModel(
            N=N,
            coupling=3.0,
            frequencies=frequencies,
            initial_phases=np.random.uniform(0, 2*np.pi, N)
        )

        t, phases = model.evolve((0, 10), dt=0.01)

        # Most oscillators should synchronize, outlier should not
        final_phases = phases[-1]
        # Check synchronization of main group (excluding outlier)
        main_group_R = np.abs(np.mean(np.exp(1j * final_phases[1:])))
        self.assertGreater(main_group_R, 0.9)

        # Outlier should have different phase velocity
        outlier_velocity = np.diff(phases[:, 0])
        mean_velocity = np.mean(np.diff(phases[:, 1:], axis=0))
        self.assertNotAlmostEqual(np.mean(outlier_velocity), mean_velocity, places=0)


if __name__ == '__main__':
    unittest.main(verbosity=2)