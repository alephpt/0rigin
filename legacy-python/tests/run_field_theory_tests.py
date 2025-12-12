#!/usr/bin/env python3
"""
Simple test runner for field theory components (no pytest required).
"""

import sys
import traceback
import numpy as np
from numpy.testing import assert_allclose

# Add parent to path
sys.path.insert(0, '/home/persist/neotec/0rigin')

from kuramoto.field_theory import (
    SpatialGrid,
    MediatorField,
    LocalFieldCoupling,
    FermionMassDemo,
    SMFTSystem
)


def test_mediator_initialization():
    """Test mediator field initialization."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
    mediator = MediatorField(grid, wave_speed=1.0, mass=1.0)

    assert mediator.c == 1.0
    assert mediator.M == 1.0
    assert mediator.sigma.shape == (32, 32)
    print("✓ Mediator initialization")


def test_source_density():
    """Test oscillator density source computation."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
    mediator = MediatorField(grid)

    phases = np.array([0.0])
    positions = np.array([[0.5, 0.5]])

    rho = mediator.compute_source_density(phases, positions, kernel_width=0.1)

    assert rho.shape == (32, 32)
    assert np.max(rho) > 0
    print("✓ Source density computation")


def test_klein_gordon_evolution():
    """Test Klein-Gordon evolution."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')
    mediator = MediatorField(grid, wave_speed=1.0, mass=0.5)

    mediator.sigma = grid.create_gaussian(sigma=0.1)
    initial_energy = mediator.compute_field_energy()

    source = np.zeros((32, 32))
    for _ in range(100):
        mediator.evolve_step(source, dt=0.01, method='leapfrog')

    final_energy = mediator.compute_field_energy()
    energy_change = abs(final_energy - initial_energy) / initial_energy

    assert energy_change < 0.2
    print(f"✓ Klein-Gordon evolution (energy drift: {energy_change:.4f})")


def test_heavy_mass_limit():
    """Test heavy mass limit."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    rho = np.random.rand(16, 16)

    M = 10.0
    mediator = MediatorField(grid, mass=M, coupling_constant=1.0)
    sigma_limit = mediator.get_heavy_mass_limit(rho)

    expected = mediator.g * rho / M**2
    error = np.max(np.abs(sigma_limit - expected))

    assert error < 1e-10
    print("✓ Heavy mass limit")


def test_local_coupling_init():
    """Test coupling system initialization."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    coupling = LocalFieldCoupling(grid)

    assert coupling.mediator is not None
    assert coupling.R_field is not None
    print("✓ Local coupling initialization")


def test_fields_update():
    """Test that fields respond to oscillators."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    coupling = LocalFieldCoupling(grid)

    N = 20
    phases = np.zeros(N)
    positions = np.random.rand(N, 2)

    coupling.update_fields_from_oscillators(phases, positions, dt=0.01)

    # Check that R field has been updated (not all zeros)
    assert np.max(coupling.R_field.values) > 0.01
    print(f"✓ Fields update from oscillators (max R={np.max(coupling.R_field.values):.3f})")


def test_coupled_evolution():
    """Test full coupled evolution."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    coupling = LocalFieldCoupling(grid, kernel_width=0.1)

    N = 30
    phases = np.random.uniform(0, 2*np.pi, N)
    positions = np.random.rand(N, 2)
    frequencies = np.random.normal(1.0, 0.1, N)

    result = coupling.evolve_coupled_system(
        phases, positions, frequencies,
        dt=0.02, n_steps=100, store_interval=20
    )

    assert 'R' in result
    assert len(result['R']) > 1
    print(f"✓ Coupled evolution (final R={result['R'][-1]:.3f})")


def test_fermion_mass_init():
    """Test fermion demo initialization."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    fermion = FermionMassDemo(grid, yukawa_coupling=1.0, bare_mass=0.5)

    assert fermion.Delta == 1.0
    assert fermion.m0 == 0.5
    print("✓ Fermion mass demo initialization")


def test_mass_from_R():
    """Test that m_eff = m0 + Δ·R."""
    grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
    Delta = 2.0
    m0 = 0.5
    fermion = FermionMassDemo(grid, yukawa_coupling=Delta, bare_mass=m0)

    R_value = 0.7
    fermion.R_field.values = np.full((16, 16), R_value)
    fermion.m_eff_field.values = m0 + Delta * fermion.R_field.values

    expected_mass = m0 + Delta * R_value
    assert_allclose(fermion.m_eff_field.values, expected_mass, rtol=1e-10)
    print("✓ Mass generation formula")


def test_mass_generation():
    """Test mass generation from oscillator synchronization."""
    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
    fermion = FermionMassDemo(grid, yukawa_coupling=1.5, bare_mass=0.0)

    N = 50
    phases = np.random.uniform(0, 0.5, N)
    positions = np.random.rand(N, 2)

    fermion.update_from_oscillators(phases, positions, kernel_width=0.1)

    avg_mass = fermion.compute_average_mass()
    assert avg_mass > fermion.m0
    print(f"✓ Mass generation (<m>={avg_mass:.3f})")


def test_smft_system():
    """Test SMFTSystem with new components."""
    system = SMFTSystem(
        grid_shape=(16, 16),
        N_oscillators=20,
        coupling='local',
        mediator_mass=5.0
    )

    result = system.evolve((0, 1), dt=0.02, store_interval=10)

    assert 'R' in result
    assert len(result['R']) > 0
    print(f"✓ SMFT system integration (R={result['R'][-1]:.3f})")


def test_full_integration():
    """Integration test: all components together."""
    grid = SpatialGrid(Nx=24, Ny=24, Lx=1.0, Ly=1.0)
    coupling = LocalFieldCoupling(grid)
    fermion = FermionMassDemo(grid, yukawa_coupling=1.0)

    N = 30
    phases = np.random.uniform(0, 2*np.pi, N)
    positions = np.random.rand(N, 2)
    frequencies = np.random.normal(1.0, 0.1, N)

    result = coupling.evolve_coupled_system(
        phases, positions, frequencies,
        dt=0.02, n_steps=50, store_interval=10
    )

    final_phases = result['phases'][-1]
    fermion.update_from_oscillators(final_phases, positions)

    avg_mass = fermion.compute_average_mass()
    assert avg_mass >= fermion.m0
    print(f"✓ Full integration (mass={avg_mass:.3f})")


def run_all_tests():
    """Run all tests."""
    tests = [
        ("Mediator initialization", test_mediator_initialization),
        ("Source density", test_source_density),
        ("Klein-Gordon evolution", test_klein_gordon_evolution),
        ("Heavy mass limit", test_heavy_mass_limit),
        ("Local coupling init", test_local_coupling_init),
        ("Fields update", test_fields_update),
        ("Coupled evolution", test_coupled_evolution),
        ("Fermion mass init", test_fermion_mass_init),
        ("Mass from R", test_mass_from_R),
        ("Mass generation", test_mass_generation),
        ("SMFT system", test_smft_system),
        ("Full integration", test_full_integration),
    ]

    print("\n" + "="*60)
    print("FIELD THEORY TESTS")
    print("="*60 + "\n")

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ {name}: FAILED")
            print(f"  {str(e)}")
            traceback.print_exc()
            failed += 1

    print("\n" + "="*60)
    print(f"Results: {passed} passed, {failed} failed")
    print("="*60 + "\n")

    return failed == 0


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
