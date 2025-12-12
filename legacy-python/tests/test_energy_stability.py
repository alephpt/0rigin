"""
Energy Stability Validation Test

Tests that verify the energy explosion bug fix from Sprint 3.
"""

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')

import numpy as np
import pytest
from kuramoto.field_theory.fields import MediatorField, SpatialGrid


def test_energy_stability_with_damping():
    """Test that energy remains bounded with damping enabled."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=16.0, Ly=16.0)
    field = MediatorField(
        grid=grid,
        wave_speed=1.0,
        mass=1.0,
        coupling_constant=1.0,
        damping=0.1  # With damping
    )

    # Initialize with random field
    field.sigma = np.random.randn(32, 32) * 0.1
    field.sigma_dot = np.random.randn(32, 32) * 0.1

    E0 = field.compute_field_energy()
    print(f"\nInitial energy: {E0:.6f}")

    # Evolve for 100 steps (no external source)
    source = np.zeros((32, 32))
    for i in range(100):
        field.evolve_step(source_density=source, dt=0.01)

    E_final = field.compute_field_energy()
    print(f"Final energy: {E_final:.6f}")
    print(f"Energy ratio: {E_final/E0:.6f}")
    print(f"Energy drift: {abs(E_final - E0)/E0 * 100:.2f}%")

    # Check bounds
    assert not np.isnan(E_final), "Energy is NaN"
    assert not np.isinf(E_final), "Energy is Inf"
    assert E_final < 1000 * E0, f"Energy exploded: {E_final/E0}x initial"
    print("✓ Energy stability test PASSED")


def test_energy_stability_no_damping():
    """Test that energy remains bounded without damping (conservative system)."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=16.0, Ly=16.0)
    field = MediatorField(
        grid=grid,
        wave_speed=1.0,
        mass=1.0,
        coupling_constant=1.0,
        damping=0.0  # No damping - conservative
    )

    # Initialize with random field
    field.sigma = np.random.randn(32, 32) * 0.1
    field.sigma_dot = np.random.randn(32, 32) * 0.1

    E0 = field.compute_field_energy()
    print(f"\nInitial energy (no damping): {E0:.6f}")

    # Evolve for 100 steps (no external source)
    source = np.zeros((32, 32))
    for i in range(100):
        field.evolve_step(source_density=source, dt=0.01)

    E_final = field.compute_field_energy()
    print(f"Final energy: {E_final:.6f}")
    print(f"Energy ratio: {E_final/E0:.6f}")
    print(f"Energy conservation error: {abs(E_final - E0)/E0 * 100:.2f}%")

    # Check bounds (more strict for conservative system)
    assert not np.isnan(E_final), "Energy is NaN"
    assert not np.isinf(E_final), "Energy is Inf"
    assert E_final < 10 * E0, f"Energy grew too much: {E_final/E0}x initial"
    print("✓ Energy conservation test PASSED")


def test_field_values_remain_bounded():
    """Test that field values themselves don't explode."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=16.0, Ly=16.0)
    field = MediatorField(
        grid=grid,
        wave_speed=1.0,
        mass=1.0,
        coupling_constant=1.0,
        damping=0.1
    )

    # Initialize with random field
    field.sigma = np.random.randn(32, 32) * 0.1
    field.sigma_dot = np.random.randn(32, 32) * 0.1

    sigma_0_max = np.abs(field.sigma).max()
    print(f"\nInitial max |σ|: {sigma_0_max:.6f}")

    # Evolve for 100 steps (no external source)
    source = np.zeros((32, 32))
    for i in range(100):
        field.evolve_step(source_density=source, dt=0.01)

    sigma_final_max = np.abs(field.sigma).max()
    print(f"Final max |σ|: {sigma_final_max:.6f}")
    print(f"Field growth: {sigma_final_max/sigma_0_max:.6f}x")

    # Check field values
    assert not np.any(np.isnan(field.sigma)), "Field contains NaN"
    assert not np.any(np.isinf(field.sigma)), "Field contains Inf"
    assert sigma_final_max < 100 * sigma_0_max, f"Field exploded: {sigma_final_max/sigma_0_max}x"
    print("✓ Field boundedness test PASSED")


def test_long_evolution_stability():
    """Test stability over longer evolution time."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=8.0, Ly=8.0)  # Smaller grid for speed
    field = MediatorField(
        grid=grid,
        wave_speed=1.0,
        mass=1.0,
        coupling_constant=1.0,
        damping=0.1
    )

    # Initialize with random field
    field.sigma = np.random.randn(16, 16) * 0.05
    field.sigma_dot = np.random.randn(16, 16) * 0.05

    E0 = field.compute_field_energy()
    print(f"\nInitial energy (long run): {E0:.6f}")

    # Evolve for 1000 steps (no external source)
    source = np.zeros((16, 16))
    max_energy = E0
    for i in range(1000):
        field.evolve_step(source_density=source, dt=0.01)
        E = field.compute_field_energy()
        if E > max_energy:
            max_energy = E

    E_final = field.compute_field_energy()
    print(f"Final energy: {E_final:.6f}")
    print(f"Max energy reached: {max_energy:.6f}")
    print(f"Max energy ratio: {max_energy/E0:.6f}x")

    # Check that energy never exploded during evolution
    assert not np.isnan(E_final), "Energy is NaN"
    assert not np.isinf(E_final), "Energy is Inf"
    assert max_energy < 1000 * E0, f"Energy exploded during evolution: {max_energy/E0}x"
    print("✓ Long evolution stability test PASSED")


if __name__ == "__main__":
    test_energy_stability_with_damping()
    test_energy_stability_no_damping()
    test_field_values_remain_bounded()
    test_long_evolution_stability()
    print("\n" + "="*60)
    print("ALL ENERGY STABILITY TESTS PASSED")
    print("="*60)
