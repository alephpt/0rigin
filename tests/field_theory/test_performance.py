"""
Tests for performance benchmarks and integration examples.
"""

import pytest
import numpy as np
import time
from numpy.testing import assert_allclose

from src.kuramoto.field_theory import SMFTSystem


class TestPerformance:
    """Test performance benchmarks."""

    def test_small_system_performance(self):
        """Test small system meets performance target."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='local',
            mediator_mass=10.0
        )

        start = time.time()
        system.step(dt=0.01)
        elapsed = time.time() - start

        # Single step should be fast (<100ms)
        assert elapsed < 0.1, f"Step too slow: {elapsed:.3f}s"

    def test_medium_system_performance(self):
        """Test medium system performance."""
        system = SMFTSystem(
            grid_shape=(50, 50),
            N_oscillators=200,
            coupling='local',
            mediator_mass=10.0
        )

        start = time.time()
        system.step(dt=0.01)
        elapsed = time.time() - start

        # Should complete in reasonable time (<1s)
        assert elapsed < 1.0, f"Step too slow: {elapsed:.3f}s"


def test_integration_example():
    """
    Integration test demonstrating full workflow.

    This is the example from the task specification.
    """
    # Create system
    system = SMFTSystem(
        grid_shape=(100, 100),
        N_oscillators=200,
        coupling='local',
        mediator_mass=10.0
    )

    # Evolve
    result = system.evolve(t_span=(0, 20), dt=0.01, store_interval=50)

    # Access fields
    R_field = result['sync_field']
    theta_field = result['theta']
    sigma_field = result['mediator_field']
    m_eff = system.compute_effective_mass()

    # Verify outputs
    assert R_field.shape[1:] == (100, 100)
    assert len(theta_field) > 0
    assert sigma_field.shape[1:] == (100, 100)
    assert m_eff.shape == (100, 100)