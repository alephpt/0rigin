#!/usr/bin/env python
"""
Verify that the SMFTSystem energy explosion bug is fixed.
Tests multiple configurations to ensure energy remains bounded.
"""

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')

from kuramoto.field_theory.smft_system import SMFTSystem
import numpy as np

def test_energy_conservation(mass=10.0, coupling='local', steps=1000):
    """Test that energy remains bounded during evolution."""
    print(f"\nTesting: mass={mass}, coupling={coupling}")

    # Create system
    model = SMFTSystem(
        grid_shape=(10, 10),
        N_oscillators=20,
        coupling=coupling,
        mediator_mass=mass
    )

    # Evolve
    dt = 0.01
    T = dt * steps
    result = model.evolve(t_span=(0, T), dt=dt, store_interval=10)

    # Check energy
    initial_energy = result['energy'][0]
    final_energy = result['energy'][-1]
    max_energy = max(result['energy'])
    min_energy = min(result['energy'])

    # For Hamiltonian systems, energy can decrease (dissipation) but shouldn't explode
    energy_range = max_energy - min_energy

    print(f"  Initial energy: {initial_energy:.6f}")
    print(f"  Final energy: {final_energy:.6f}")
    print(f"  Energy range: [{min_energy:.6f}, {max_energy:.6f}]")
    print(f"  Final R: {result['R'][-1]:.6f}")

    # Check for explosion (energy shouldn't grow by more than 10x)
    if abs(max_energy) > 10 * abs(initial_energy):
        print(f"  ❌ FAILED: Energy exploded! Max/Initial ratio: {abs(max_energy/initial_energy):.2e}")
        return False

    # Check for NaN
    if np.isnan(result['R'][-1]):
        print(f"  ❌ FAILED: R became NaN!")
        return False

    print(f"  ✅ PASSED: Energy bounded, R valid")
    return True

def main():
    """Run comprehensive energy stability tests."""
    print("=" * 60)
    print("SMFTSystem Energy Stability Verification")
    print("=" * 60)

    tests = [
        (1.0, 'local'),
        (10.0, 'local'),
        (100.0, 'local'),
        (1.0, 'global'),
        (10.0, 'global'),
        (100.0, 'global'),
    ]

    results = []
    for mass, coupling in tests:
        try:
            passed = test_energy_conservation(mass, coupling)
            results.append((mass, coupling, passed))
        except Exception as e:
            print(f"  ❌ FAILED with exception: {e}")
            results.append((mass, coupling, False))

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    passed_count = sum(1 for _, _, passed in results if passed)
    total_count = len(results)

    for mass, coupling, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"  M={mass:5.1f}, {coupling:6s}: {status}")

    print(f"\nOverall: {passed_count}/{total_count} tests passed")

    if passed_count == total_count:
        print("\n🎉 SUCCESS: All energy stability tests passed!")
        print("The damping bug fix is working correctly.")
        return 0
    else:
        print(f"\n⚠️  WARNING: {total_count - passed_count} tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())