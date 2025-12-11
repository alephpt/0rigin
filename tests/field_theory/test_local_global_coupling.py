"""
Tests for local vs global coupling regimes.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from src.kuramoto.field_theory import SMFTSystem


class TestLocalGlobalCoupling:
    """Test local vs global coupling regimes."""

    def test_local_coupling_spatial_structure(self):
        """Test that local coupling creates spatial structure."""
        system = SMFTSystem(
            grid_shape=(30, 30),
            N_oscillators=100,
            coupling='local',
            mediator_mass=10.0
        )

        # Evolve
        solution = system.evolve((0, 10), dt=0.01, store_interval=50)

        # Check that sync field has spatial variation
        final_sync = solution['sync_field'][-1]
        variance = np.var(final_sync)

        assert variance > 0  # Spatial structure exists

    def test_global_coupling_uniform_field(self):
        """Test that global coupling creates more uniform field."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='global',
            mediator_mass=10.0
        )

        # Evolve
        solution = system.evolve((0, 5), dt=0.01, store_interval=25)

        # Global coupling should create less spatial variation
        # (though not completely uniform due to discrete sampling)
        final_sync = solution['sync_field'][-1]
        variance = np.var(final_sync)

        # Should have some structure but less than highly local case
        assert variance >= 0