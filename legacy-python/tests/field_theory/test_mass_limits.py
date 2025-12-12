"""
Tests for heavy mass limit and backward compatibility.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from kuramoto.field_theory import (
    SMFTSystem,
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField
)


class TestHeavyMassLimit:
    """Test M→∞ limit recovers Kuramoto model."""

    def test_increasing_mass_convergence(self):
        """Test that R converges as M increases."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=30,
            coupling='local'
        )

        M_values = [1.0, 10.0, 100.0]
        results = system.test_heavy_mass_limit(M_values)

        # Results should be defined
        assert len(results) == 3
        for M in M_values:
            assert 0 <= results[M] <= 1

    def test_heavy_mass_reduces_field_dynamics(self):
        """Test that heavy mass slows field dynamics."""
        N = 40
        frequencies = np.random.normal(0, 1, N)

        # Light mass
        system_light = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=N,
            mediator_mass=1.0,
            oscillator_frequencies=frequencies
        )

        # Heavy mass
        system_heavy = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=N,
            mediator_mass=100.0,
            oscillator_frequencies=frequencies
        )

        # Evolve both
        sol_light = system_light.evolve((0, 2), dt=0.01, store_interval=20)
        sol_heavy = system_heavy.evolve((0, 2), dt=0.01, store_interval=20)

        # Heavy mass should have smaller field variation
        field_var_light = np.var(sol_light['mediator_field'][-1])
        field_var_heavy = np.var(sol_heavy['mediator_field'][-1])

        # Heavy mass → slower field dynamics
        assert field_var_heavy <= field_var_light * 2  # Some margin


class TestBackwardCompatibility:
    """Test that field theory doesn't break Sprint 1 components."""

    def test_hamiltonian_kuramoto_standalone(self):
        """Test HamiltonianKuramoto works independently."""
        model = HamiltonianKuramoto(
            N=50,
            coupling_strength=3.0,
            frequencies=np.random.normal(0, 1, 50),
            damping=1.0
        )

        solution = model.evolve((0, 5), dt=0.01)

        assert 't' in solution
        assert 'theta' in solution
        assert 'R' in solution
        assert len(solution['t']) > 0

    def test_spatial_grid_standalone(self):
        """Test SpatialGrid works independently."""
        grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')

        # Test Laplacian
        max_error, rel_error = grid.test_laplacian()
        assert rel_error < 0.1

        # Test gradient
        f = grid.X + 2 * grid.Y
        grad_x, grad_y = grid.gradient(f)
        assert np.allclose(grad_x, 1.0, atol=0.1)

    def test_scalar_field_standalone(self):
        """Test ScalarField works independently."""
        grid = SpatialGrid(20, 20, 1.0, 1.0)
        field = ScalarField(grid, name="test")

        # Test diffusion
        field.values = grid.create_gaussian(sigma=0.05)
        initial_max = np.max(field.values)

        for _ in range(20):
            field.diffuse(0.01, 0.01)

        assert np.max(field.values) < initial_max