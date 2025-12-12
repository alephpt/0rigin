"""
Tests for discrete-continuum integration.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from kuramoto.field_theory import SMFTSystem


class TestDiscreteToContiuumIntegration:
    """Test discrete oscillators → continuous field mapping."""

    def test_local_order_parameter_computation(self):
        """Test that local R(x,y) is computed correctly."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=10,
            coupling='local'
        )

        # Set all oscillators to same phase
        system.oscillators.theta[:] = np.pi / 4

        # Compute local order parameter
        R_field, theta_field = system.compute_local_order_parameter(kernel_width=0.2)

        # Check shapes
        assert R_field.shape == (20, 20)
        assert theta_field.shape == (20, 20)

        # With all oscillators synchronized, R should be positive
        # Note: with sparse sampling (10 oscillators on 20x20 grid), R is lower
        assert np.mean(R_field) > 0.01

    def test_local_order_parameter_bounds(self):
        """Test that local R ∈ [0,1]."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=50,
            coupling='local'
        )

        R_field, theta_field = system.compute_local_order_parameter()

        # Check bounds
        assert np.all(R_field >= 0)
        assert np.all(R_field <= 1)

        # Check phase bounds
        assert np.all(theta_field >= -np.pi)
        assert np.all(theta_field <= np.pi)

    def test_continuum_limit_convergence(self):
        """Test convergence as N increases."""
        grid_shape = (30, 30)
        N_values = [50, 100, 200]

        variances = []
        for N in N_values:
            system = SMFTSystem(
                grid_shape=grid_shape,
                N_oscillators=N,
                coupling='local'
            )

            # Synchronized state
            system.oscillators.theta[:] = 0

            R_field, _ = system.compute_local_order_parameter(kernel_width=0.1)

            # Variance should decrease with more oscillators
            variances.append(np.var(R_field))

        # Check decreasing variance
        assert variances[1] < variances[0]
        assert variances[2] < variances[1]