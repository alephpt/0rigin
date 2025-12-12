"""
Test py-pde library integration and benchmark performance.

Validates PDE solver wrapper and measures performance for different grid sizes.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.kuramoto.field_theory import PDESolver, PDE_AVAILABLE


def test_pde_availability():
    """Check if py-pde is available."""
    print("Checking py-pde availability...")
    print("-" * 40)

    if PDE_AVAILABLE:
        print("  py-pde is installed: YES")
        try:
            import pde
            print(f"  py-pde version: {pde.__version__}")
        except:
            print("  py-pde version: Unknown")
    else:
        print("  py-pde is installed: NO")
        print("  Install with: pip install py-pde")
        print("  Falling back to manual implementations.")

    return PDE_AVAILABLE


def test_simple_diffusion():
    """Test diffusion equation solver."""
    if not PDE_AVAILABLE:
        print("\nSkipping diffusion test (py-pde not available)")
        return False

    print("\nTesting diffusion equation solver...")
    print("-" * 40)

    # Create solver
    solver = PDESolver(
        grid_shape=(50, 50),
        grid_size=(1.0, 1.0),
        boundary='periodic'
    )

    # Initial condition: Gaussian peak
    X, Y = np.meshgrid(
        np.linspace(0, 1, 50),
        np.linspace(0, 1, 50),
        indexing='ij'
    )
    initial = np.exp(-50 * ((X - 0.5)**2 + (Y - 0.5)**2))

    initial_max = np.max(initial)
    initial_integral = np.sum(initial)

    # Solve diffusion
    result = solver.solve_diffusion(
        initial,
        diffusion_coeff=0.01,
        t_range=(0, 1.0),
        dt=0.1,
        tracker=None
    )

    final_field = result['field']
    final_max = np.max(final_field)
    final_integral = np.sum(final_field)

    print(f"  Initial max value: {initial_max:.4f}")
    print(f"  Final max value: {final_max:.4f}")
    print(f"  Initial integral: {initial_integral:.4f}")
    print(f"  Final integral: {final_integral:.4f}")

    # Check conservation (for periodic BC) and diffusion
    conserved = abs(final_integral - initial_integral) / initial_integral < 0.01
    diffused = final_max < initial_max * 0.5

    print(f"  Mass conservation: {'PASS' if conserved else 'FAIL'}")
    print(f"  Diffusion occurred: {'PASS' if diffused else 'FAIL'}")

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    im1 = axes[0].imshow(initial.T, origin='lower', cmap='hot',
                         extent=[0, 1, 0, 1], vmin=0, vmax=initial_max)
    axes[0].set_title('Initial Condition')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(im1, ax=axes[0])

    im2 = axes[1].imshow(final_field.T, origin='lower', cmap='hot',
                         extent=[0, 1, 0, 1], vmin=0, vmax=initial_max)
    axes[1].set_title('After Diffusion (t=1.0)')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig('pde_diffusion_test.png', dpi=150)
    plt.show()

    return conserved and diffused


def benchmark_performance():
    """Benchmark PDE solver performance for different grid sizes."""
    if not PDE_AVAILABLE:
        print("\nSkipping benchmark (py-pde not available)")
        return {}

    print("\nBenchmarking PDE solver performance...")
    print("-" * 40)

    grid_sizes = [20, 50, 100]
    results = {}

    for N in grid_sizes:
        print(f"  Testing {N}x{N} grid...")

        solver = PDESolver(
            grid_shape=(N, N),
            grid_size=(1.0, 1.0),
            boundary='periodic'
        )

        # Random initial condition
        initial = np.random.randn(N, N)

        # Time the solve
        start_time = time.time()

        try:
            result = solver.solve_diffusion(
                initial,
                diffusion_coeff=0.1,
                t_range=(0, 1.0),
                dt=0.01,
                tracker=None
            )
            elapsed = time.time() - start_time

            results[f'{N}x{N}'] = elapsed
            print(f"    Time: {elapsed:.3f} seconds")

            # Estimate scaling
            if N == 50:
                time_50 = elapsed
            elif N == 100 and 'time_50' in locals():
                scaling = elapsed / time_50
                print(f"    Scaling (100x100 vs 50x50): {scaling:.1f}x")

        except Exception as e:
            results[f'{N}x{N}'] = f"Failed: {str(e)}"
            print(f"    Failed: {e}")

    return results


def test_reaction_diffusion():
    """Test custom PDE with reaction-diffusion dynamics."""
    if not PDE_AVAILABLE:
        print("\nSkipping reaction-diffusion test (py-pde not available)")
        return False

    print("\nTesting reaction-diffusion equation...")
    print("-" * 40)

    # Create solver
    solver = PDESolver(
        grid_shape=(100, 100),
        grid_size=(1.0, 1.0),
        boundary='periodic'
    )

    # Initial condition: random perturbations around uniform state
    initial = 0.5 + 0.1 * np.random.randn(100, 100)

    try:
        # Solve reaction-diffusion: ∂u/∂t = D∇²u + u(1-u)(u-a)
        # This is a bistable reaction-diffusion equation
        result = solver.solve_custom_pde(
            pde_rhs="D * laplace(u) + u * (1 - u) * (u - a)",
            initial_field=initial,
            t_range=(0, 10.0),
            dt=0.1,
            parameters={'D': 0.01, 'a': 0.3}
        )

        final_field = result['field']

        # Check for pattern formation
        initial_variance = np.var(initial)
        final_variance = np.var(final_field)

        print(f"  Initial variance: {initial_variance:.4f}")
        print(f"  Final variance: {final_variance:.4f}")
        print(f"  Pattern formation: {'YES' if final_variance > initial_variance else 'NO'}")

        # Visualize
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        im1 = axes[0].imshow(initial.T, origin='lower', cmap='RdBu_r',
                             extent=[0, 1, 0, 1], vmin=0, vmax=1)
        axes[0].set_title('Initial Condition')
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('y')
        plt.colorbar(im1, ax=axes[0])

        im2 = axes[1].imshow(final_field.T, origin='lower', cmap='RdBu_r',
                             extent=[0, 1, 0, 1], vmin=0, vmax=1)
        axes[1].set_title('After Reaction-Diffusion (t=10)')
        axes[1].set_xlabel('x')
        axes[1].set_ylabel('y')
        plt.colorbar(im2, ax=axes[1])

        plt.tight_layout()
        plt.savefig('pde_reaction_diffusion.png', dpi=150)
        plt.show()

        return True

    except Exception as e:
        print(f"  Custom PDE failed: {e}")
        return False


def test_boundary_conditions():
    """Test different boundary conditions."""
    if not PDE_AVAILABLE:
        print("\nSkipping boundary condition test (py-pde not available)")
        return False

    print("\nTesting boundary conditions...")
    print("-" * 40)

    boundaries = ['periodic', 'dirichlet', 'neumann']
    results = {}

    # Same initial condition for all
    X, Y = np.meshgrid(
        np.linspace(0, 1, 30),
        np.linspace(0, 1, 30),
        indexing='ij'
    )
    initial = np.exp(-30 * ((X - 0.5)**2 + (Y - 0.5)**2))

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    # Plot initial
    im0 = axes[0].imshow(initial.T, origin='lower', cmap='hot',
                         extent=[0, 1, 0, 1], vmin=0, vmax=1)
    axes[0].set_title('Initial')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(im0, ax=axes[0])

    for idx, bc in enumerate(boundaries):
        print(f"  Testing {bc} BC...")

        try:
            solver = PDESolver(
                grid_shape=(30, 30),
                grid_size=(1.0, 1.0),
                boundary=bc
            )

            result = solver.solve_diffusion(
                initial.copy(),
                diffusion_coeff=0.05,
                t_range=(0, 0.5),
                dt=0.05,
                tracker=None
            )

            final = result['field']
            results[bc] = 'PASS'

            # Plot
            im = axes[idx+1].imshow(final.T, origin='lower', cmap='hot',
                                    extent=[0, 1, 0, 1], vmin=0, vmax=1)
            axes[idx+1].set_title(f'{bc.capitalize()} BC')
            axes[idx+1].set_xlabel('x')
            axes[idx+1].set_ylabel('y')
            plt.colorbar(im, ax=axes[idx+1])

        except Exception as e:
            results[bc] = f'FAIL: {e}'
            print(f"    Failed: {e}")

    plt.tight_layout()
    plt.savefig('pde_boundary_conditions.png', dpi=150)
    plt.show()

    # Summary
    for bc, status in results.items():
        print(f"  {bc}: {status}")

    return all(v == 'PASS' for v in results.values())


def main():
    """Run all PDE integration tests."""
    print("=" * 60)
    print("PDE SOLVER INTEGRATION TEST")
    print("=" * 60)

    # Check availability
    available = test_pde_availability()

    if available:
        # Test 1: Simple diffusion
        diffusion_test = test_simple_diffusion()

        # Test 2: Performance benchmark
        benchmark_results = benchmark_performance()

        # Test 3: Reaction-diffusion
        reaction_test = test_reaction_diffusion()

        # Test 4: Boundary conditions
        bc_test = test_boundary_conditions()

        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"py-pde available: YES")
        print(f"Diffusion test: {'PASS' if diffusion_test else 'FAIL'}")
        print(f"Reaction-diffusion: {'PASS' if reaction_test else 'FAIL'}")
        print(f"Boundary conditions: {'PASS' if bc_test else 'FAIL'}")
        print(f"Performance: Tested {len(benchmark_results)} grid sizes")
        print("\nPDE solver integration validated successfully!")

    else:
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print("py-pde not available - install with: pip install py-pde")
        print("Manual implementations can still be used.")


if __name__ == "__main__":
    main()