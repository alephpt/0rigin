"""
Tests for SMFT system initialization and evolution.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from src.kuramoto.field_theory import (
    SMFTSystem,
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField,
    MediatorField,
    LocalFieldCoupling
)


class TestSMFTSystemInitialization:
    """Test SMFT system initialization and setup."""

    def test_system_creation(self):
        """Test basic system creation."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='local',
            mediator_mass=10.0
        )

        assert system.grid_shape == (20, 20)
        assert system.N == 50
        assert system.coupling_type == 'local'
        assert system.M == 10.0
        assert system.t == 0.0

    def test_component_creation(self):
        """Test that all components are created."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=30,
            coupling='global'
        )

        # Check grid
        assert system.grid is not None
        assert system.grid.Nx == 15
        assert system.grid.Ny == 15

        # Check oscillators
        assert system.oscillators is not None
        assert system.oscillators.N == 30

        # Check fields
        assert system.mediator_field is not None
        assert system.sync_field is not None
        assert system.phase_field is not None

        # Check positions
        assert system.oscillator_positions.shape == (30, 2)

    def test_custom_frequencies(self):
        """Test initialization with custom frequencies."""
        N = 25
        frequencies = np.random.normal(0, 2, N)

        system = SMFTSystem(
            grid_shape=(10, 10),
            N_oscillators=N,
            oscillator_frequencies=frequencies
        )

        np.testing.assert_array_equal(
            system.oscillators.frequencies,
            frequencies
        )


class TestFullSystemEvolution:
    """Test full integrated evolution with field theory."""

    def test_short_evolution(self):
        """Test short time evolution maintains consistency."""
        system = SMFTSystem(
            grid_shape=(10, 10),
            N_oscillators=20,
            coupling='local',
            g=1.0,
            mediator_mass=5.0
        )

        initial_phases = system.oscillators.theta.copy()
        initial_mediator = system.mediator_field.values.copy()

        # Evolve for short time
        t, states = system.evolve(duration=0.1, dt=0.01)

        # Check time array
        assert len(t) > 0
        assert t[-1] >= 0.1

        # Check phases changed
        assert not np.allclose(
            system.oscillators.theta,
            initial_phases
        )

        # Check mediator evolved
        assert not np.allclose(
            system.mediator_field.values,
            initial_mediator
        )

    def test_synchronization_emergence(self):
        """Test that synchronization can emerge under strong coupling."""
        N = 50
        # Use narrow frequency distribution
        frequencies = np.random.normal(0, 0.5, N)

        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=N,
            oscillator_frequencies=frequencies,
            coupling='global',
            g=10.0,  # Strong coupling
            mediator_mass=1.0  # Light mediator
        )

        # Evolve system
        t, states = system.evolve(duration=5.0, dt=0.01, save_every=50)

        # Extract final order parameter
        final_R, _ = system.compute_order_parameter()

        # Should have some synchronization
        assert final_R > 0.3  # Reasonable sync level

    @pytest.mark.slow
    def test_phase_transition(self):
        """Test phase transition as coupling strength increases."""
        grid_shape = (10, 10)
        N = 30

        coupling_strengths = np.linspace(0.1, 5.0, 10)
        final_R_values = []

        for g in coupling_strengths:
            system = SMFTSystem(
                grid_shape=grid_shape,
                N_oscillators=N,
                coupling='global',
                g=g,
                mediator_mass=2.0
            )

            # Evolve to steady state
            system.evolve(duration=2.0, dt=0.01, save_every=100)

            R, _ = system.compute_order_parameter()
            final_R_values.append(R)

        # Order parameter should increase with coupling
        correlation = np.corrcoef(coupling_strengths, final_R_values)[0, 1]
        assert correlation > 0.5  # Positive correlation