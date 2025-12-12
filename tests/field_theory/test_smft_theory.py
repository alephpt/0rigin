"""
Tests validating SMFT theory alignment.

Tests the corrected mass formula implementation:
    m_eff = Δ · R  (NOT m ∝ 1/R)

Where:
    Δ = mass gap parameter (energy scale)
    R = synchronization order parameter [0,1]

Physical interpretation:
    High sync R → high mass (trapped energy E=mc²)
    Low sync R → low mass (approaches massless)
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from kuramoto.field_theory import SMFTSystem


class TestMassFormula:
    """Test correct SMFT mass formula: m = Δ·R."""

    def test_mass_proportional_to_sync(self):
        """Verify m_eff ∝ R (not m ∝ 1/R)."""
        system = SMFTSystem(N=10, mass_gap=2.0)

        # Set known sync values
        system.sync_field.values[:] = 0.5
        m_eff = system.compute_effective_mass()

        # Should be: m = Δ·R = 2.0 * 0.5 = 1.0
        expected = 2.0 * 0.5
        assert_allclose(m_eff, expected, err_msg=f"Expected {expected}, got {m_eff.mean()}")

    def test_mass_zero_at_zero_sync(self):
        """Test R=0 → m=0 (massless limit)."""
        system = SMFTSystem(N=10, mass_gap=2.0)

        # Zero sync
        system.sync_field.values[:] = 0.0
        m_eff = system.compute_effective_mass()

        assert_allclose(m_eff, 0.0, err_msg="R=0 should give m=0 (massless fermions)")

    def test_mass_maximum_at_full_sync(self):
        """Test R=1 → m=Δ (maximum mass)."""
        system = SMFTSystem(N=10, mass_gap=2.0)

        # Full sync
        system.sync_field.values[:] = 1.0
        m_eff = system.compute_effective_mass()

        assert_allclose(m_eff, 2.0, err_msg=f"R=1 should give m=Δ={2.0}")

    def test_mass_linear_scaling(self):
        """Test mass scales linearly with R."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        # Low sync
        system.sync_field.values[:] = 0.1
        m_low = system.compute_effective_mass().mean()

        # High sync
        system.sync_field.values[:] = 0.9
        m_high = system.compute_effective_mass().mean()

        # Should scale linearly
        ratio = m_high / m_low
        expected_ratio = 0.9 / 0.1

        assert_allclose(ratio, expected_ratio, rtol=1e-3,
                       err_msg=f"Mass should scale linearly with R: got ratio {ratio}, expected {expected_ratio}")


class TestMassGapParameter:
    """Test mass_gap (Δ) parameter."""

    def test_default_mass_gap(self):
        """Test default mass_gap value."""
        system = SMFTSystem(N=10)
        assert system.Delta == 1.0, "Default mass_gap should be 1.0"

    def test_custom_mass_gap(self):
        """Test custom mass_gap value."""
        system = SMFTSystem(N=10, mass_gap=5.0)
        assert system.Delta == 5.0, "Should use provided mass_gap"

    def test_mass_gap_validation_zero(self):
        """Test validation rejects zero mass_gap."""
        with pytest.raises(ValueError, match="mass_gap must be positive"):
            SMFTSystem(N=10, mass_gap=0.0)

    def test_mass_gap_validation_negative(self):
        """Test validation rejects negative mass_gap."""
        with pytest.raises(ValueError, match="mass_gap must be positive"):
            SMFTSystem(N=10, mass_gap=-1.0)

    def test_mass_gap_scales_mass(self):
        """Test that mass_gap scales the effective mass."""
        # System with Delta=1.0
        system1 = SMFTSystem(N=10, mass_gap=1.0)
        system1.sync_field.values[:] = 0.5
        m1 = system1.compute_effective_mass().mean()

        # System with Delta=3.0
        system2 = SMFTSystem(N=10, mass_gap=3.0)
        system2.sync_field.values[:] = 0.5
        m2 = system2.compute_effective_mass().mean()

        # Should scale by Delta ratio
        assert_allclose(m2 / m1, 3.0, rtol=1e-10,
                       err_msg="Mass should scale linearly with mass_gap")


class TestChiralDecomposition:
    """Test chiral mass decomposition: m = Δ·R·[cos(θ) + iγ⁵sin(θ)]."""

    def test_pure_scalar_mass(self):
        """Test pure scalar mass at theta=0."""
        system = SMFTSystem(N=10, mass_gap=2.0)
        system.sync_field.values[:] = 0.5

        m_s, m_p = system.compute_chiral_mass(theta=0.0)

        # Scalar mass: Δ·R·cos(0) = 2.0 * 0.5 * 1 = 1.0
        assert_allclose(m_s, 1.0, err_msg="Scalar mass should be Δ·R·cos(0) = 1.0")

        # Pseudoscalar should be zero
        assert_allclose(m_p, 0.0, atol=1e-15, err_msg="Pseudoscalar should be 0 at theta=0")

    def test_pure_pseudoscalar_mass(self):
        """Test pure pseudoscalar mass at theta=π/2."""
        system = SMFTSystem(N=10, mass_gap=2.0)
        system.sync_field.values[:] = 0.5

        m_s, m_p = system.compute_chiral_mass(theta=np.pi/2)

        # Scalar should be zero
        assert_allclose(m_s, 0.0, atol=1e-10, err_msg="Scalar should be 0 at theta=π/2")

        # Pseudoscalar: Δ·R·sin(π/2) = 2.0 * 0.5 * 1 = 1.0
        assert_allclose(m_p, 1.0, err_msg="Pseudoscalar should be Δ·R·sin(π/2) = 1.0")

    def test_equal_components(self):
        """Test equal scalar/pseudoscalar components at theta=π/4."""
        system = SMFTSystem(N=10, mass_gap=2.0)
        system.sync_field.values[:] = 0.5

        m_s, m_p = system.compute_chiral_mass(theta=np.pi/4)

        # Both should be Δ·R/√2 = 1.0/√2
        expected = 1.0 / np.sqrt(2)
        assert_allclose(m_s, expected, rtol=1e-10, err_msg=f"Expected scalar {expected}")
        assert_allclose(m_p, expected, rtol=1e-10, err_msg=f"Expected pseudoscalar {expected}")

        # Should be equal
        assert_allclose(m_s, m_p, rtol=1e-10, err_msg="Components should be equal at π/4")

    def test_orthogonality_condition(self):
        """Verify orthogonality: m_s² + m_p² = (Δ·R)²."""
        system = SMFTSystem(N=10, mass_gap=2.0)
        system.sync_field.values[:] = 0.5

        # Test at various angles
        angles = [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2]

        for theta in angles:
            m_s, m_p = system.compute_chiral_mass(theta=theta)

            # Check orthogonality
            total_sq = m_s**2 + m_p**2
            expected_sq = (2.0 * 0.5)**2  # (Δ·R)²

            assert_allclose(total_sq, expected_sq, rtol=1e-10,
                           err_msg=f"Orthogonality violated at theta={theta}: {total_sq} != {expected_sq}")

    def test_chiral_angle_sweep(self):
        """Test chiral decomposition across full angle range."""
        system = SMFTSystem(N=10, mass_gap=1.5)
        system.sync_field.values[:] = 0.8

        angles = np.linspace(0, 2*np.pi, 20)

        for theta in angles:
            m_s, m_p = system.compute_chiral_mass(theta=theta)

            # Verify bounds
            total_mass = 1.5 * 0.8
            assert np.all(np.abs(m_s) <= total_mass), "Scalar component out of bounds"
            assert np.all(np.abs(m_p) <= total_mass), "Pseudoscalar component out of bounds"

            # Verify orthogonality
            assert_allclose(m_s**2 + m_p**2, total_mass**2, rtol=1e-10,
                           err_msg=f"Orthogonality violated at theta={theta}")


class TestTheoryConsistency:
    """Test overall theory consistency."""

    def test_high_sync_gives_high_mass(self):
        """Verify high sync → high mass (not inverse)."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        # Low sync
        system.sync_field.values[:] = 0.1
        m_low = system.compute_effective_mass().mean()

        # High sync
        system.sync_field.values[:] = 0.9
        m_high = system.compute_effective_mass().mean()

        assert m_high > m_low, "High sync should give higher mass (not lower!)"

    def test_mass_monotonic_with_sync(self):
        """Test mass increases monotonically with sync."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        R_values = np.linspace(0, 1, 11)
        masses = []

        for R in R_values:
            system.sync_field.values[:] = R
            m = system.compute_effective_mass().mean()
            masses.append(m)

        masses = np.array(masses)

        # Check monotonicity
        diffs = np.diff(masses)
        assert np.all(diffs >= 0), "Mass should increase monotonically with sync"

        # Check linearity
        expected = R_values
        assert_allclose(masses, expected, rtol=1e-10, err_msg="Mass should be linear in R")

    def test_mass_field_stored(self):
        """Test that compute_effective_mass stores result."""
        system = SMFTSystem(N=10, mass_gap=1.0)
        system.sync_field.values[:] = 0.7

        # Initially None
        assert system.mass_field is None, "Mass field should be None initially"

        # Compute
        m_eff = system.compute_effective_mass()

        # Should be stored
        assert system.mass_field is not None, "Mass field should be stored after compute"
        assert_allclose(system.mass_field, m_eff, err_msg="Stored mass should match computed")

    def test_spatial_variation(self):
        """Test that spatially varying sync gives spatially varying mass."""
        system = SMFTSystem(N=10, grid_shape=(10, 10), mass_gap=1.0)

        # Create spatial variation in sync
        X, Y = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))
        system.sync_field.values = 0.5 + 0.3 * np.sin(2*np.pi*X)

        m_eff = system.compute_effective_mass()

        # Mass should vary spatially
        assert m_eff.std() > 0, "Mass should vary spatially when sync varies"

        # Check correlation with sync
        assert np.corrcoef(system.sync_field.values.flatten(), m_eff.flatten())[0,1] > 0.99, \
            "Mass should be highly correlated with sync"


class TestBugRegression:
    """Regression tests for the m=M/(R+ε) bug."""

    def test_no_division_by_sync(self):
        """Ensure mass formula doesn't divide by sync (old bug)."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        # Near-zero sync should give near-zero mass (not infinity!)
        system.sync_field.values[:] = 1e-10
        m_eff = system.compute_effective_mass()

        # Should be near zero, not large
        assert np.all(m_eff < 1e-8), "Near-zero sync should give near-zero mass, not infinity"

        # Should be proportional to R
        expected = 1.0 * 1e-10
        assert_allclose(m_eff, expected, rtol=1e-5, err_msg="Should be m=Δ·R, not m=M/(R+ε)")

    def test_no_epsilon_regularization(self):
        """Ensure no epsilon regularization is used (old bug)."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        # Exact zero sync
        system.sync_field.values[:] = 0.0
        m_eff = system.compute_effective_mass()

        # Should be exactly zero (no epsilon added)
        assert_allclose(m_eff, 0.0, atol=1e-15, err_msg="R=0 should give exactly m=0")

    def test_correct_direction_of_proportionality(self):
        """Test mass increases (not decreases) with sync."""
        system = SMFTSystem(N=10, mass_gap=1.0)

        # Test multiple sync levels
        sync_levels = [0.2, 0.4, 0.6, 0.8]
        masses = []

        for R in sync_levels:
            system.sync_field.values[:] = R
            m = system.compute_effective_mass().mean()
            masses.append(m)

        # Should be strictly increasing
        for i in range(len(masses) - 1):
            assert masses[i+1] > masses[i], \
                f"Mass should increase with sync: {masses[i]} -> {masses[i+1]} failed"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
