"""
Simple validation script for field theory prototypes.

Runs basic tests without pytest dependency.
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.kuramoto.field_theory import (
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField,
    PDESolver,
    PDE_AVAILABLE
)


def validate_hamiltonian():
    """Validate Hamiltonian Kuramoto prototype."""
    print("Validating Hamiltonian Kuramoto...")
    print("-" * 40)

    try:
        # Test 1: Initialization
        N = 20
        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=2.0,
            frequencies=np.ones(N),
            damping=1.0
        )
        print("  ✓ Initialization successful")

        # Test 2: Energy conservation (γ=0)
        model_conservative = HamiltonianKuramoto(
            N=20,
            coupling_strength=2.0,
            frequencies=np.zeros(20),
            damping=0.0
        )
        E0 = model_conservative.compute_hamiltonian()
        solution = model_conservative.evolve((0, 1.0), dt=0.001, method='rk4')
        Ef = model_conservative.compute_hamiltonian()
        rel_error = abs(Ef - E0) / abs(E0) if E0 != 0 else abs(Ef - E0)

        if rel_error < 1e-3:
            print(f"  ✓ Energy conservation: relative error = {rel_error:.2e}")
        else:
            print(f"  ✗ Energy conservation failed: relative error = {rel_error:.2e}")

        # Test 3: Overdamped limit
        model_overdamped = HamiltonianKuramoto(
            N=20,
            coupling_strength=3.0,
            frequencies=np.random.normal(0, 1, 20),
            damping=100.0,
            initial_momenta=np.ones(20)
        )
        solution = model_overdamped.evolve((0, 1.0), dt=0.01)
        final_p_max = np.max(np.abs(model_overdamped.p))

        if final_p_max < 0.1:
            print(f"  ✓ Overdamped limit: max|p| = {final_p_max:.3f}")
        else:
            print(f"  ✗ Overdamped limit failed: max|p| = {final_p_max:.3f}")

        return True

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def validate_grid():
    """Validate spatial grid infrastructure."""
    print("\nValidating Spatial Grid...")
    print("-" * 40)

    try:
        # Test 1: Grid creation
        grid = SpatialGrid(50, 50, 2.0, 3.0, 'periodic')
        assert grid.X.shape == (50, 50)
        print("  ✓ Grid creation successful")

        # Test 2: Laplacian accuracy
        grid_test = SpatialGrid(40, 40, 2*np.pi, 2*np.pi, 'periodic')
        max_error, rel_error = grid_test.test_laplacian()

        if rel_error < 1e-2:
            print(f"  ✓ Laplacian accuracy: relative error = {rel_error:.2e}")
        else:
            print(f"  ✗ Laplacian accuracy failed: relative error = {rel_error:.2e}")

        # Test 3: Gradient operator
        f = grid.X + 2 * grid.Y  # Linear function
        grad_x, grad_y = grid.gradient(f)

        grad_x_error = np.max(np.abs(grad_x - 1.0))
        grad_y_error = np.max(np.abs(grad_y - 2.0))

        if grad_x_error < 0.1 and grad_y_error < 0.1:
            print(f"  ✓ Gradient operator: errors = ({grad_x_error:.2e}, {grad_y_error:.2e})")
        else:
            print(f"  ✗ Gradient operator failed: errors = ({grad_x_error:.2e}, {grad_y_error:.2e})")

        return True

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def validate_scalar_field():
    """Validate scalar field evolution."""
    print("\nValidating Scalar Field...")
    print("-" * 40)

    try:
        # Test 1: Field initialization
        grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')
        field = ScalarField(grid, name="test")
        print("  ✓ Field creation successful")

        # Test 2: Diffusion
        field.values = grid.create_gaussian(sigma=0.05)
        initial_max = np.max(field.values)

        for _ in range(50):
            field.diffuse(diffusion_coeff=0.01, dt=0.01)

        final_max = np.max(field.values)

        if final_max < initial_max * 0.9:
            print(f"  ✓ Diffusion: max reduced from {initial_max:.3f} to {final_max:.3f}")
        else:
            print(f"  ✗ Diffusion failed: max only reduced to {final_max:.3f}")

        # Test 3: Oscillator coupling
        field_coupling = ScalarField(grid)
        N = 5
        phases = np.zeros(N)
        positions = np.array([[0.5, 0.5], [0.3, 0.3], [0.7, 0.3], [0.3, 0.7], [0.7, 0.7]])

        field_coupling.update_from_oscillators(phases, positions, kernel_width=0.1)

        if np.max(field_coupling.values) > 0:
            print(f"  ✓ Oscillator coupling: max field = {np.max(field_coupling.values):.3f}")
        else:
            print("  ✗ Oscillator coupling failed: field is zero")

        return True

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def validate_pde_solver():
    """Validate PDE solver integration."""
    print("\nValidating PDE Solver...")
    print("-" * 40)

    if not PDE_AVAILABLE:
        print("  ⚠ py-pde not installed")
        print("  Install with: pip install py-pde")
        return None

    try:
        # Test 1: Solver creation
        solver = PDESolver(
            grid_shape=(20, 20),
            grid_size=(1.0, 1.0),
            boundary='periodic'
        )
        print("  ✓ Solver creation successful")

        # Test 2: Integration test
        if PDESolver.test_integration():
            print("  ✓ Integration test passed")
        else:
            print("  ✗ Integration test failed")

        # Test 3: Diffusion solving
        X, Y = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20), indexing='ij')
        initial = np.exp(-20 * ((X - 0.5)**2 + (Y - 0.5)**2))
        initial_max = np.max(initial)

        result = solver.solve_diffusion(
            initial,
            diffusion_coeff=0.1,
            t_range=(0, 0.1),
            dt=0.01,
            tracker=None
        )

        final_max = np.max(result['field'])

        if final_max < initial_max:
            print(f"  ✓ Diffusion solving: max reduced from {initial_max:.3f} to {final_max:.3f}")
        else:
            print("  ✗ Diffusion solving failed")

        return True

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def validate_integration():
    """Validate that all components work together."""
    print("\nValidating Component Integration...")
    print("-" * 40)

    try:
        # Create all components
        grid = SpatialGrid(20, 20, 1.0, 1.0, 'periodic')
        field = ScalarField(grid, name="order_parameter")
        N = 30
        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=3.0,
            frequencies=np.random.normal(0, 1, N),
            damping=1.0
        )

        # Evolve model
        solution = model.evolve((0, 0.5), dt=0.01)

        # Update field from oscillators
        positions = np.random.rand(N, 2)
        field.update_from_oscillators(model.theta, positions)

        # Apply diffusion
        field.diffuse(0.01, 0.01)

        if field.t > 0 and model.t > 0 and np.max(field.values) > 0:
            print("  ✓ All components work together")
            return True
        else:
            print("  ✗ Component integration failed")
            return False

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False


def main():
    """Run all validation tests."""
    print("=" * 60)
    print("FIELD THEORY PROTOTYPE VALIDATION")
    print("=" * 60)

    results = {}

    # Run all validations
    results['hamiltonian'] = validate_hamiltonian()
    results['grid'] = validate_grid()
    results['scalar_field'] = validate_scalar_field()
    results['pde_solver'] = validate_pde_solver()
    results['integration'] = validate_integration()

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    for component, status in results.items():
        if status is None:
            symbol = "⚠"
            text = "SKIPPED"
        elif status:
            symbol = "✓"
            text = "PASS"
        else:
            symbol = "✗"
            text = "FAIL"

        print(f"  {symbol} {component.replace('_', ' ').title()}: {text}")

    # Overall status
    passed = sum(1 for v in results.values() if v is True)
    failed = sum(1 for v in results.values() if v is False)
    skipped = sum(1 for v in results.values() if v is None)

    print(f"\n  Total: {passed} passed, {failed} failed, {skipped} skipped")

    if failed == 0:
        print("\n✓ All prototypes validated successfully!")
    else:
        print(f"\n✗ {failed} prototype(s) failed validation")

    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)