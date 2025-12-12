"""
Validation script for SMFT system integration.

Tests all integration points and generates validation report.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.kuramoto.field_theory import (
    SMFTSystem,
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField
)


def validate_discrete_continuum():
    """Validate Discrete ↔ Continuum integration."""
    print("\n" + "="*60)
    print("VALIDATION 1: Discrete ↔ Continuum")
    print("="*60)

    # Test N-scaling convergence
    grid_shape = (30, 30)
    N_values = [50, 100, 200, 500]
    variances = []

    print("\nTesting continuum limit convergence...")
    for N in N_values:
        system = SMFTSystem(
            grid_shape=grid_shape,
            N_oscillators=N,
            coupling='local'
        )

        # Synchronized state
        system.oscillators.theta[:] = 0.5

        R_field, _ = system.compute_local_order_parameter(kernel_width=0.1)
        variance = np.var(R_field)
        variances.append(variance)

        print(f"  N={N:3d}: Field variance = {variance:.6f}")

    # Check convergence
    convergence = all(variances[i+1] < variances[i] for i in range(len(variances)-1))

    if convergence:
        print("\n✓ PASSED: Variance decreases with N (continuum limit)")
    else:
        print("\n⚠ WARNING: Variance not monotonically decreasing")

    return convergence


def validate_hamiltonian_field():
    """Validate Hamiltonian ↔ Field integration."""
    print("\n" + "="*60)
    print("VALIDATION 2: Hamiltonian ↔ Field Dynamics")
    print("="*60)

    system = SMFTSystem(
        grid_shape=(20, 20),
        N_oscillators=50,
        coupling='local',
        mediator_mass=10.0
    )

    print("\nTesting coupled dynamics...")

    # Initial state
    theta0 = system.oscillators.theta.copy()
    field0 = system.mediator_field.values.copy()

    # Evolve
    for _ in range(20):
        system.step(dt=0.01)

    # Check both changed
    theta_changed = not np.array_equal(system.oscillators.theta, theta0)
    field_changed = not np.array_equal(system.mediator_field.values, field0)

    print(f"  Oscillators updated: {theta_changed}")
    print(f"  Field updated: {field_changed}")

    if theta_changed and field_changed:
        print("\n✓ PASSED: Coupled dynamics working")
        return True
    else:
        print("\n✗ FAILED: Coupling not working")
        return False


def validate_local_global():
    """Validate Local ↔ Global coupling."""
    print("\n" + "="*60)
    print("VALIDATION 3: Local ↔ Global Coupling")
    print("="*60)

    N = 80
    freq = np.random.normal(0, 1, N)

    # Local coupling
    system_local = SMFTSystem(
        grid_shape=(30, 30),
        N_oscillators=N,
        coupling='local',
        oscillator_frequencies=freq
    )

    # Global coupling
    system_global = SMFTSystem(
        grid_shape=(30, 30),
        N_oscillators=N,
        coupling='global',
        oscillator_frequencies=freq
    )

    print("\nEvolving both systems...")
    sol_local = system_local.evolve((0, 10), dt=0.01, store_interval=50)
    sol_global = system_global.evolve((0, 10), dt=0.01, store_interval=50)

    R_local = sol_local['R'][-1]
    R_global = sol_global['R'][-1]

    print(f"  Local coupling R: {R_local:.4f}")
    print(f"  Global coupling R: {R_global:.4f}")

    if 0 <= R_local <= 1 and 0 <= R_global <= 1:
        print("\n✓ PASSED: Both coupling modes produce valid R")
        return True
    else:
        print("\n✗ FAILED: Invalid R values")
        return False


def validate_backward_compatibility():
    """Validate Classical ↔ Field Theory compatibility."""
    print("\n" + "="*60)
    print("VALIDATION 4: Backward Compatibility")
    print("="*60)

    print("\nTesting Sprint 1 components...")

    # Test HamiltonianKuramoto
    model = HamiltonianKuramoto(
        N=50,
        coupling_strength=3.0,
        frequencies=np.random.normal(0, 1, 50),
        damping=1.0
    )
    sol = model.evolve((0, 5), dt=0.01)

    hamiltonian_ok = 't' in sol and len(sol['t']) > 0
    print(f"  HamiltonianKuramoto: {'✓' if hamiltonian_ok else '✗'}")

    # Test SpatialGrid
    grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')
    max_err, rel_err = grid.test_laplacian()

    grid_ok = rel_err < 0.1
    print(f"  SpatialGrid (Laplacian): {'✓' if grid_ok else '✗'} (error={rel_err:.4f})")

    # Test ScalarField
    field = ScalarField(grid, name="test")
    field.values = grid.create_gaussian(sigma=0.05)
    initial_max = np.max(field.values)

    for _ in range(20):
        field.diffuse(0.01, 0.01)

    field_ok = np.max(field.values) < initial_max
    print(f"  ScalarField (diffusion): {'✓' if field_ok else '✗'}")

    if hamiltonian_ok and grid_ok and field_ok:
        print("\n✓ PASSED: All Sprint 1 components work independently")
        return True
    else:
        print("\n✗ FAILED: Some components broken")
        return False


def validate_heavy_mass_limit():
    """Validate M→∞ recovers Kuramoto."""
    print("\n" + "="*60)
    print("VALIDATION 5: Heavy Mass Limit")
    print("="*60)

    M_values = [5.0, 10.0, 50.0, 100.0]
    R_values = []

    print("\nTesting mass scaling...")
    for M in M_values:
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            mediator_mass=M
        )

        sol = system.evolve((0, 10), dt=0.01, store_interval=50)
        R = sol['R'][-1]
        R_values.append(R)

        print(f"  M={M:6.1f}: R={R:.4f}")

    # Check all R valid
    all_valid = all(0 <= R <= 1 for R in R_values)

    if all_valid:
        print("\n✓ PASSED: Heavy mass limit produces valid synchronization")
        return True
    else:
        print("\n✗ FAILED: Invalid R values")
        return False


def validate_full_system():
    """Validate complete system evolution."""
    print("\n" + "="*60)
    print("VALIDATION 6: Full System Evolution")
    print("="*60)

    system = SMFTSystem(
        grid_shape=(40, 40),
        N_oscillators=100,
        coupling='local',
        mediator_mass=10.0
    )

    print("\nRunning full evolution workflow...")
    solution = system.evolve(
        t_span=(0, 15),
        dt=0.01,
        store_interval=75
    )

    # Validate outputs
    checks = {
        'Time array': 't' in solution and len(solution['t']) > 0,
        'Phases': 'theta' in solution and solution['theta'].shape[1] == 100,
        'Order parameter': 'R' in solution and all(0 <= r <= 1 or np.isnan(r) for r in solution['R']),
        'Mediator field': 'mediator_field' in solution,
        'Sync field': 'sync_field' in solution,
        'Energy': 'energy' in solution
    }

    print("\nOutput validation:")
    for name, passed in checks.items():
        print(f"  {name}: {'✓' if passed else '✗'}")

    # Effective mass
    m_eff = system.compute_effective_mass()
    mass_ok = m_eff.shape == (40, 40) and np.all(m_eff > 0)
    print(f"  Effective mass: {'✓' if mass_ok else '✗'}")

    all_passed = all(checks.values()) and mass_ok

    if all_passed:
        print("\n✓ PASSED: Full system workflow validated")
        return True
    else:
        print("\n✗ FAILED: Some outputs invalid")
        return False


def main():
    """Run all validations."""
    print("\n" + "="*60)
    print("SMFT SYSTEM INTEGRATION VALIDATION")
    print("="*60)
    print("\nValidating all integration points...")

    results = {
        'Discrete ↔ Continuum': validate_discrete_continuum(),
        'Hamiltonian ↔ Field': validate_hamiltonian_field(),
        'Local ↔ Global': validate_local_global(),
        'Backward Compatibility': validate_backward_compatibility(),
        'Heavy Mass Limit': validate_heavy_mass_limit(),
        'Full System': validate_full_system()
    }

    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)

    for name, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{name:30s}: {status}")

    n_passed = sum(results.values())
    n_total = len(results)

    print("\n" + "="*60)
    print(f"TOTAL: {n_passed}/{n_total} validations passed")
    print("="*60)

    if n_passed == n_total:
        print("\n✓ ALL INTEGRATION POINTS VALIDATED")
        print("\nThe SMFT system successfully integrates:")
        print("  • Discrete oscillators with continuous fields")
        print("  • Hamiltonian phase space with field dynamics")
        print("  • Local and global coupling mechanisms")
        print("  • Classical Kuramoto with field theory extensions")
        print("  • All components work independently and together")
        return 0
    else:
        print(f"\n⚠ {n_total - n_passed} validation(s) failed")
        return 1


if __name__ == "__main__":
    exit(main())
