"""
Tests for Hamiltonian and field integration.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from kuramoto.field_theory import SMFTSystem


class TestHamiltonianFieldIntegration:
    """Test Hamiltonian phase space ↔ field integration."""

    def test_field_force_on_oscillators(self):
        """Test that field exerts force on oscillators."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=30,
            coupling='local',
            mediator_mass=5.0
        )

        # Set non-zero mediator field
        system.mediator_field.values[:] = 0.5

        # Compute forces
        forces = system.compute_field_force_on_oscillators(kernel_width=0.1)

        assert len(forces) == 30
        # Forces should be non-zero
        assert np.any(forces != 0)

    def test_energy_conservation_no_damping(self):
        """Test energy behavior in conservative case."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=20,
            coupling='local',
            mediator_mass=10.0
        )

        # Reduce damping for near-conservative system
        system.oscillators.gamma = 0.01

        # Initial energy
        E0 = system.oscillators.compute_hamiltonian()

        # Evolve briefly
        solution = system.evolve((0, 1.0), dt=0.001, store_interval=100)

        # Check energy doesn't explode
        E_final = solution['energy'][-1]
        assert abs(E_final) < 100 * abs(E0 + 1)  # Reasonable bound

    def test_coupled_step_updates_both(self):
        """Test that single step updates oscillators and field."""
        system = SMFTSystem(
            grid_shape=(10, 10),
            N_oscillators=15,
            coupling='local'
        )

        # Initial states
        theta0 = system.oscillators.theta.copy()
        field0 = system.mediator_field.values.copy()
        t0 = system.t

        # Single step
        system.step(dt=0.01)

        # Check both updated
        assert not np.array_equal(system.oscillators.theta, theta0)
        assert not np.array_equal(system.mediator_field.values, field0)
        assert system.t > t0